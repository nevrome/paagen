{-# LANGUAGE OverloadedStrings #-}

import           Paths_paagen                   (version)

import           Control.Applicative            ((<|>))
import           Control.Exception              (catch, throwIO, Exception)
import           Control.Monad                  (forM, guard)
import           Control.Monad.Random           (fromList, RandomGen, evalRand, newStdGen, getStdGen)
import           Data.Char                      (isSpace)
import           Data.Function                  (on)
import           Data.List                      (nub, tails, sortBy, intersect, maximumBy, group, sort, intercalate, elemIndex, elemIndices, unfoldr)
import           Data.Maybe                     (isJust, fromMaybe, catMaybes, fromJust)
import qualified Data.Vector                    as V
import           Data.Version                   (showVersion)
import qualified Options.Applicative            as OP
import           Pipes
import qualified Pipes.Prelude                  as P
import           Pipes.Safe                     (SafeT (..), runSafeT, throwM)
import           Poseidon.GenotypeData
import           Poseidon.Janno
import           Poseidon.Package
import           SequenceFormats.Eigenstrat     (EigenstratIndEntry (..),
                                                EigenstratSnpEntry (..), GenoLine,
                                                writeEigenstrat, GenoEntry (..))
import           SequenceFormats.Plink          (writePlink)
import           System.Console.ANSI            (hClearLine, hSetCursorColumn)
import           System.Exit                    (exitFailure)
import           System.FilePath                ((<.>), (</>))
import           System.IO                      (hPutStrLn, stderr, hPutStr)
import qualified Text.Parsec                    as P
import qualified Text.Parsec.String             as P
import qualified Text.Parsec.Number             as P
import Control.Exception (SomeException(SomeException))

-- data types

data IndWithPosition = IndWithPosition {
      ind :: String
    , unit :: String
    , pos :: SpatialTemporalPosition
} deriving (Show)

data SpatialTemporalPosition = SpatialTemporalPosition {
      time :: Int
    , lat :: Latitude 
    , lon :: Longitude
} deriving (Show)

instance Ord GenoEntry where
    compare HomRef Het     = GT
    compare Het HomAlt     = GT
    compare HomAlt Missing = GT
    compare _ _            = EQ

data PaagenException =
    PaagenCLIParsingException String
    deriving (Show)

instance Exception PaagenException

data SpaceTimeOptions = SpaceTimeOptions {   
      _optBaseDirs :: [FilePath]
    , _optIndWithPosition :: [IndWithPosition]
    , _optIndWithPositionFile :: Maybe FilePath
    , _optNumberOfNearestNeighbors :: Int
    , _optOutFormat :: GenotypeFormatSpec
    , _optOutPath :: FilePath
    }

data Options = CmdSpaceTime SpaceTimeOptions

-- CLI interface configuration

main :: IO ()
main = do
    cmdOpts <- OP.customExecParser p optParserInfo
    runCmd cmdOpts
    where
        p = OP.prefs OP.showHelpOnEmpty

runCmd :: Options -> IO ()
runCmd o = case o of
    CmdSpaceTime opts -> runSpaceTime opts

optParserInfo :: OP.ParserInfo Options
optParserInfo = OP.info (OP.helper <*> versionOption <*> optParser) (
    OP.briefDesc <>
    OP.progDesc "paagen generates artificial genotype profiles \
                \based on input data in the Poseidon format."
    )

versionOption :: OP.Parser (a -> a)
versionOption = OP.infoOption (showVersion version) (OP.long "version" <> OP.help "Show version")

optParser :: OP.Parser Options
optParser = OP.subparser (
        OP.command "spacetime" spaceTimeOptInfo
    )
  where
    spaceTimeOptInfo = OP.info (OP.helper <*> (CmdSpaceTime <$> spaceTimeOptParser))
        (OP.progDesc "Genotype profile generation based on spatiotemporal position")

spaceTimeOptParser :: OP.Parser SpaceTimeOptions
spaceTimeOptParser = SpaceTimeOptions <$> parseBasePaths
                                      <*> parseIndWithPositionDirect
                                      <*> parseIndWithPositionFromFile
                                      <*> parseNumberOfNearestNeighbors
                                      <*> parseOutGenotypeFormat
                                      <*> parseOutPath

parseBasePaths :: OP.Parser [FilePath]
parseBasePaths = OP.some (OP.strOption (
    OP.long "baseDir" <>
    OP.short 'd' <>
    OP.metavar "DIR" <>
    OP.help "A base directory to search for Poseidon Packages (could be a Poseidon repository)"
    ))

parseIndWithPositionDirect :: OP.Parser [IndWithPosition]
parseIndWithPositionDirect = OP.option (OP.eitherReader readIndWithPositionString) (
    OP.long "positionString" <>
    OP.short 'p' <>
    OP.value [] <>
    OP.help "Spatiotemporal positions of interest: each position is a string of the form \
            \\"[id:group](time,latitude,longitude)\". Multiple positions can be listed separated by ;. \
            \id and group are simple strings, latitude and longitude must be given in decimal degrees \
            \and time in calBC/AD, so 3245BC would be -3245 and 1148AD just 1148"
    )

parseIndWithPositionFromFile :: OP.Parser (Maybe FilePath)
parseIndWithPositionFromFile = OP.option (Just <$> OP.str) (OP.long "positionFile" <>
    OP.value Nothing <>
    OP.help "A file with a list of spatiotemporal positions. \
            \Works just as -p, but multiple values can also be separated by newline, not just by ;. \
            \-p and --positionFile can be combined."
    )

readIndWithPositionString :: String -> Either String [IndWithPosition]
readIndWithPositionString s = case P.runParser indWithPositionParser () "" s of
    Left p  -> Left (show p)
    Right x -> Right x

readIndWithPositionFromFile :: FilePath -> IO [IndWithPosition]
readIndWithPositionFromFile positionFile = do
    let multiPositionParser = indWithPositionParser `P.sepBy1` (P.newline *> P.spaces)
    eitherParseResult <- P.parseFromFile (P.spaces *> multiPositionParser <* P.spaces) positionFile
    case eitherParseResult of
        Left err -> throwIO $ PaagenCLIParsingException (show err)
        Right r -> return (concat r)

indWithPositionParser :: P.Parser [IndWithPosition]
indWithPositionParser = P.try (P.sepBy parseIndWithPosition (P.char ';' <* P.spaces))

parseIndWithPosition :: P.Parser IndWithPosition
parseIndWithPosition = do
    _ <- P.oneOf "["
    ind <- P.manyTill P.anyChar (P.string ":")
    unit <- P.manyTill P.anyChar (P.string "]")
    _ <- P.oneOf "("
    spatpos <- parseSpatialTemporalPosition
    _ <- P.oneOf ")"
    return (IndWithPosition ind unit spatpos)

parseSpatialTemporalPosition :: P.Parser SpatialTemporalPosition
parseSpatialTemporalPosition = do
    time <- pInt
    _ <- P.oneOf ","
    lat <- pLat
    _ <- P.oneOf ","
    lon <- pLon
    return (SpatialTemporalPosition time lat lon)

pInt :: P.Parser Int
pInt = read <$> P.many1 P.digit

pLat :: P.Parser Latitude 
pLat = do
    lat <- P.sign <*> P.floating2 True
    guard (lat >= -90 && lat <= 90) P.<?> "valid latitude (-90 to 90)"
    return (Latitude lat)

pLon :: P.Parser Longitude
pLon = do
    lon <- P.sign <*> P.floating2 True
    guard (lon >= -180 && lon <= 180) P.<?> "valid longitude (-180 to 180)"
    return (Longitude lon)

parseNumberOfNearestNeighbors :: OP.Parser Int
parseNumberOfNearestNeighbors = OP.option OP.auto (
    OP.long "neighbors" <> 
    OP.help "Number of nearest neighbors to consider for the calculation" <>
    OP.value 50 <>
    OP.showDefault
    )

parseOutGenotypeFormat :: OP.Parser GenotypeFormatSpec
parseOutGenotypeFormat = OP.option (OP.eitherReader readGenotypeFormat) (
    OP.long "outFormat" <>
    OP.help "The format of the output genotype data: EIGENSTRAT or PLINK" <>
    OP.value GenotypeFormatPlink
    )
    where
    readGenotypeFormat :: String -> Either String GenotypeFormatSpec
    readGenotypeFormat s = case s of
        "EIGENSTRAT" -> Right GenotypeFormatEigenstrat
        "PLINK"      -> Right GenotypeFormatPlink
        _            -> Left "must be EIGENSTRAT or PLINK"

parseOutPath :: OP.Parser FilePath
parseOutPath = OP.strOption (
    OP.long "outPath" <>
    OP.short 'o' <>
    OP.help "The output directory path"
    )

-- Actual program code

runSpaceTime :: SpaceTimeOptions -> IO ()
runSpaceTime (SpaceTimeOptions baseDirs poisDirect poisFile numNeighbors outFormat outDir) = do
    -- compile pois
    poisFromFile <- case poisFile of
        Nothing -> return []
        Just f -> readIndWithPositionFromFile f
    let pois = poisDirect ++ poisFromFile --this nub could also be relevant for forge
    -- load Poseidon packages
    allPackages <- readPoseidonPackageCollection True True False baseDirs
    -- load janno tables
    let jannos = concatMap posPacJanno allPackages
    -- transform to spatiotemporal positions
        stInds = jannoToSpaceTimePos jannos
    -- calculate distances
        distancesToPois = map (`distanceOneToAll` stInds) pois 
    -- get X closest inds
        closest = map (getClosestInds numNeighbors) distancesToPois
        closestIndividuals = map (map fst) closest
        closestDistances = map (map snd) closest
        closestWeights = map distToWeight closestDistances
    -- putStrLn $ "Closest individuals: " ++ intercalate ", " closestIndividuals
    -- determine relevant packages
    relevantPackages <- filterPackagesByInds (nub $ concat closestIndividuals) allPackages
    closestIndices <- mapM (`extractIndIndices` relevantPackages) closestIndividuals
    -- merge poi indices and weights
    let infoForIndividualPOIs = zip closestIndices closestWeights
    -- compile genotype data structure
    let [outInd, outSnp, outGeno] = case outFormat of 
            GenotypeFormatEigenstrat -> ["poi.ind", "poi.snp", "poi.geno"]
            GenotypeFormatPlink -> ["poi.fam", "poi.bim", "poi.bed"]
    -- compile genotype data
    runSafeT $ do
        (eigenstratIndEntries, eigenstratProd) <- getJointGenotypeData False False relevantPackages
        let [outG, outS, outI] = map (outDir </>) [outGeno, outSnp, outInd]
            newIndEntries = map (\x -> EigenstratIndEntry (ind x) Unknown (unit x)) pois
        let outConsumer = case outFormat of
                GenotypeFormatEigenstrat -> writeEigenstrat outG outS outI newIndEntries 
                GenotypeFormatPlink      -> writePlink      outG outS outI newIndEntries
        runEffect $ eigenstratProd >-> printSNPCopyProgress >-> P.mapM (sampleGenoForMultiplePOIs infoForIndividualPOIs) >-> outConsumer
        liftIO $ hClearLine stderr
        liftIO $ hSetCursorColumn stderr 0
        liftIO $ hPutStrLn stderr "SNPs processed: All done"

sampleGenoForMultiplePOIs :: [([Int], [Rational])] -> (EigenstratSnpEntry, GenoLine) -> SafeT IO (EigenstratSnpEntry, GenoLine)
sampleGenoForMultiplePOIs infoForIndividualPOIs (snpEntry, genoLine) = do
    entries <- mapM (\(x,y) -> sampleGenoForOnePOI x y genoLine) infoForIndividualPOIs
    return (snpEntry, V.fromList entries)

sampleGenoForOnePOI :: [Int] -> [Rational] -> GenoLine -> SafeT IO GenoEntry
sampleGenoForOnePOI individualIndices weights genoLine = do
    let relevantGenoEntries = [genoLine V.! i | i <- individualIndices]
        -- count occurrence of GenoEntries
        genoEntryIndices = getGenoIndices relevantGenoEntries
        -- sum distance-based weight for each GenoEntry
        weightsPerGenoEntry = sumWeights genoEntryIndices weights
    -- sample GenoEntry based on weight
    gen <- liftIO getStdGen
    -- liftIO $ hPutStrLn stderr (show gen)
    let selectedGenoEntry = sampleWeightedList gen weightsPerGenoEntry
    liftIO newStdGen
    -- return 
    return selectedGenoEntry

sampleWeightedList :: RandomGen g => g -> [(a, Rational)] -> a
sampleWeightedList gen weights = head $ evalRand m gen
    where m = sequence . repeat . fromList $ weights

getGenoIndices :: Eq a => [a] -> [(a, [Int])]
getGenoIndices xs = 
    let unique = nub xs
        indices = map (\v -> elemIndices v xs) unique
    in  zip unique indices

sumWeights :: Num b => [(a, [Int])] -> [b] -> [(a, b)]
sumWeights xs weights = map (\(x, ys) -> (x, sum $ subset ys weights)) xs
    where
        subset :: [Int] -> [a] -> [a]
        subset indices xs = [xs !! i | i <- indices]

jannoToSpaceTimePos :: [JannoRow] -> [IndWithPosition]
jannoToSpaceTimePos jannos =
    let filterPos x = (isJust (jDateBCADMedian x) || jDateType x == Just Modern) && isJust (jLatitude x) && isJust (jLongitude x)
        transformTo x = IndWithPosition (jIndividualID x) (head $ jGroupName x) $
            SpatialTemporalPosition 
                (fromMaybe 2000 (jDateBCADMedian x)) -- If empty then modern (after filter), which is approx. 2000AD
                (fromMaybe (Latitude 0) (jLatitude x))
                (fromMaybe (Longitude 0) (jLongitude x))
    in  map transformTo $ filter filterPos jannos

distanceOneToAll :: IndWithPosition -> [IndWithPosition] -> [(String, String, Double)]
distanceOneToAll poi = 
    map (distanceOneToOne poi) 
    where distanceOneToOne i1 i2 = (ind i1, ind i2, spatioTemporalDistance 1 (pos i1) (pos i2))

spatioTemporalDistance :: Double -> SpatialTemporalPosition -> SpatialTemporalPosition -> Double 
spatioTemporalDistance 
    scaling
    (SpatialTemporalPosition t1 (Latitude lat1) (Longitude lon1)) 
    (SpatialTemporalPosition t2 (Latitude lat2) (Longitude lon2)) =
    let tDist = fromIntegral (abs (t1 - t2)) * scaling
        sDist = haversineDist (lat1, lon1) (lat2, lon2)
    in sqrt $ tDist^2 + sDist^2

haversineDist :: (Double, Double) -> (Double, Double) -> Double 
haversineDist (lat1, lon1) (lat2, lon2) =
    let r = 6371000  -- radius of Earth in meters
        toRadians n = n * pi / 180
        square x = x * x
        cosr = cos . toRadians
        dlat = toRadians (lat1 - lat2) / 2
        dlon = toRadians (lon1 - lon2) / 2
        a = square (sin dlat) + cosr lat1 * cosr lat2 * square (sin dlon)
        c = 2 * atan2 (sqrt a) (sqrt (1 - a))
    in (r * c) / 1000

getClosestInds :: Int -> [(String, String, Double)] -> [(String, Double)]
getClosestInds n dists = map (\(_,x,y) -> (x,y)) $ take n $ sortBy (\(_,_,x) (_,_,y) -> compare x y) dists

distToWeight :: [Double] -> [Rational]
distToWeight distances = 
    let closeness = map (1/) distances
    in rescale 0 100 closeness

rescale :: Double -> Double -> [Double] -> [Rational]
rescale minNew maxNew xs = 
    let minOld = minimum xs
        maxOld = maximum xs
        a = (maxNew - minNew) / (maxOld - minOld)
        bs = map (\x -> x - maxOld) xs
        res = map (\x -> a * x + maxNew) bs
    in map toRational res

filterPackagesByInds :: [String] -> [PoseidonPackage] -> IO [PoseidonPackage]
filterPackagesByInds indNamesStats packages = do
    fmap catMaybes . forM packages $ \pac -> do
        inds <- getIndividuals pac
        let indNamesPac   = [ind   | EigenstratIndEntry ind _ _     <- inds]
        if  length (intersect indNamesPac indNamesStats) > 0
        then return (Just pac)
        else return Nothing

extractIndIndices :: [String] -> [PoseidonPackage] -> IO [Int]
extractIndIndices indNames relevantPackages = do
    let allPackageNames = map posPacTitle relevantPackages
    allIndEntries <- mapM getIndividuals relevantPackages
    let filterFunc (_ , pacName, EigenstratIndEntry ind _ group) = ind `elem` indNames
    return $ map extractFirst $ filter filterFunc (zipGroup allPackageNames allIndEntries)

extractFirst :: (a, b, c) -> a
extractFirst (a,_,_) = a

zipGroup :: [a] -> [[b]] -> [(Int,a,b)]
zipGroup list nestedList =
    let lenghtsNestedList = map length nestedList
        listWithlenghtsNestedList = zip lenghtsNestedList list
        longerA = map (uncurry replicate) listWithlenghtsNestedList
    in zip3 [0..] (concat longerA) (concat nestedList)
