{-# LANGUAGE OverloadedStrings #-}

import           Paths_paagen                   (version)

import           Control.Applicative            ((<|>))
import           Control.Exception              (catch)
import           Control.Monad                  (forM)
import           Data.Function                  (on)
import           Data.List                      (nub, tails, sortBy, intersect, maximumBy, group, sort, intercalate, elemIndex, elemIndices)
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
import           System.IO                      (hPutStrLn, stderr)

-- CLI interface configuration

data GenOptions = GenOptions {   
      _optBaseDirs :: [FilePath]
    , _optSpatioTemporalPosition :: SpatialTemporalPosition
    , _optNumberOfNearestNeighbors :: Int
    , _optOutFormat :: GenotypeFormatSpec
    , _optOutPath :: FilePath
    }

data Options = CmdGen GenOptions

main :: IO ()
main = do
    cmdOpts <- OP.customExecParser p optParserInfo
    runCmd cmdOpts
    where
        p = OP.prefs OP.showHelpOnEmpty

runCmd :: Options -> IO ()
runCmd o = case o of
    CmdGen opts -> runGen opts

optParserInfo :: OP.ParserInfo Options
optParserInfo = OP.info (OP.helper <*> versionOption <*> optParser) (
    OP.briefDesc <>
    OP.progDesc "paagen generates artificial genotype profiles for spatiotemporal positions \
                \based on input data in the Poseidon format."
    )

versionOption :: OP.Parser (a -> a)
versionOption = OP.infoOption (showVersion version) (OP.long "version" <> OP.help "Show version")

optParser :: OP.Parser Options
optParser = OP.subparser (
        OP.command "gen" genOptInfo
    )
  where
    genOptInfo = OP.info (OP.helper <*> (CmdGen <$> genOptParser))
        (OP.progDesc "a first generator algorithm")

genOptParser :: OP.Parser GenOptions
genOptParser = GenOptions <$> parseBasePaths
                          <*> parseSpatialTemporalPosition
                          <*> parseNumberOfNearestNeighbors
                          <*> parseOutGenotypeFormat
                          <*> parseOutPath

parseBasePaths :: OP.Parser [FilePath]
parseBasePaths = OP.some (OP.strOption (
    OP.long "baseDir" <>
    OP.short 'd' <>
    OP.metavar "DIR" <>
    OP.help "a base directory to search for Poseidon Packages (could be a Poseidon repository)"
    ))

parseSpatialTemporalPosition :: OP.Parser SpatialTemporalPosition
parseSpatialTemporalPosition = SpatialTemporalPosition <$> timeParser <*> latParser <*> lonParser
    where
        timeParser = OP.option OP.auto (
            OP.long "time" <>
            OP.help "temporal position of point of interest: time in calBC/AD, so 3245BC would be -3245 and 1148AD just 1148"
            )
        latParser = OP.option (OP.eitherReader readLatitude) (
            OP.long "latitude" <>
            OP.help "spatial position of point of interest: latitude in decimal degrees"
            )
        readLatitude :: String -> Either String Latitude
        readLatitude s = 
            let val = read s :: Double
            in  if val < -90 || val > 90
                then Left "Latitude not in -90 to 90 degrees"
                else Right $ Latitude val
        lonParser = OP.option (OP.eitherReader readLongitude) (
            OP.long "longitude" <>
            OP.help "spatial position of point of interest: longitude in decimal degrees"
            )
        readLongitude :: String -> Either String Longitude
        readLongitude s = 
            let val = read s :: Double
            in  if val < -180 || val > 180
                then Left "Longitude not in -180 to 180 degrees"
                else Right $ Longitude val

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
    OP.help "the format of the output genotype data: EIGENSTRAT or PLINK" <>
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
    OP.help "the output directory path"
    )

-- Actual program code

runGen :: GenOptions -> IO ()
runGen (GenOptions baseDirs poi numNeighbors outFormat outDir) = do
    -- load Poseidon packages -- 
    allPackages <- readPoseidonPackageCollection True True False baseDirs
    -- load janno tables
    let jannos = concatMap posPacJanno allPackages
    -- transform to spatiotemporal positions
        stInds = jannoToSpaceTimePos jannos
    -- calculate distances
        positionOfInterest = IndsWithPosition "poi" poi
        distancesToPoi = distanceOneToAll positionOfInterest stInds
    -- get X closest inds
        closest = getClosestInds numNeighbors distancesToPoi
        closestIndividuals = map fst closest
        closestDistances = map snd closest
        closestWeights = distToWeight closestDistances
    putStrLn $ "Closest individuals: " ++ intercalate ", " closestIndividuals
    -- determine relevant packages
    relevantPackages <- filterPackagesByInds closestIndividuals allPackages
    closestIndices <- extractIndIndices closestIndividuals relevantPackages
    -- compile genotype data structure
    let [outInd, outSnp, outGeno] = case outFormat of 
            GenotypeFormatEigenstrat -> ["poi.ind", "poi.snp", "poi.geno"]
            GenotypeFormatPlink -> ["poi.fam", "poi.bim", "poi.bed"]
    -- compile genotype data
    runSafeT $ do
        (eigenstratIndEntries, eigenstratProd) <- getJointGenotypeData False False relevantPackages
        let [outG, outS, outI] = map (outDir </>) [outGeno, outSnp, outInd]
        let outConsumer = case outFormat of
                GenotypeFormatEigenstrat -> writeEigenstrat outG outS outI [EigenstratIndEntry "poi" Unknown "group_of_poi"]
                GenotypeFormatPlink -> writePlink outG outS outI [EigenstratIndEntry "poi" Unknown "group_of_poi"]
        runEffect $ eigenstratProd >-> printSNPCopyProgress >-> P.map (mergeIndividuals closestIndices closestWeights) >-> outConsumer
        liftIO $ hClearLine stderr
        liftIO $ hSetCursorColumn stderr 0
        liftIO $ hPutStrLn stderr "SNPs processed: All done"

mergeIndividuals :: [Int] -> [Int] -> (EigenstratSnpEntry, GenoLine) -> (EigenstratSnpEntry, GenoLine)
mergeIndividuals individualIndices weights (snpEntry, genoLine) = 
    let relevantGenoEntries = [genoLine V.! i | i <- individualIndices]
        -- only keep Missing if there is no other option
        genoEntriesWithoutMissingIfOthersAvailable = 
            if nub relevantGenoEntries == [Missing]
            then relevantGenoEntries
            else filter (/= Missing) relevantGenoEntries
        -- count occurrence of GenoEntries
        genoEntryIndices = getGenoIndices genoEntriesWithoutMissingIfOthersAvailable
        -- sum distance-based weight for each GenoEntry
        weightsPerGenoEntry = sumWeights genoEntryIndices weights
        -- modeGenoEntry = mostCommon relevantGenoEntries
        -- select GenoEntry with highest weight sum
        selectedGenoEntry = fst $ maximumBy (\ (_, a) (_, b) -> compare a b) weightsPerGenoEntry
    in (snpEntry, V.fromList [selectedGenoEntry])

getGenoIndices :: Eq a => [a] -> [(a, [Int])]
getGenoIndices xs = 
    let unique = nub xs
        indices = map (\v -> elemIndices v xs) unique
    in  zip unique indices

sumWeights :: [(a, [Int])] -> [Int] -> [(a, Int)]
sumWeights xs weights = map (\(x, ys) -> (x, sum $ subset ys weights)) xs
    where
        subset :: [Int] -> [a] -> [a]
        subset indices xs = [xs !! i | i <- indices]

-- mostCommon :: Ord a => [a] -> a
-- mostCommon = head . maximumBy (compare `on` length) . group . sort

data IndsWithPosition = IndsWithPosition {
      ind :: String
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

jannoToSpaceTimePos :: [JannoRow] -> [IndsWithPosition]
jannoToSpaceTimePos jannos =
    let filterPos x = (isJust (jDateBCADMedian x) || jDateType x == Just Modern) && isJust (jLatitude x) && isJust (jLongitude x)
        transformTo x = IndsWithPosition (jIndividualID x) $
            SpatialTemporalPosition 
                (fromMaybe 2000 (jDateBCADMedian x)) -- If empty then modern (after filter), which is approx. 2000AD
                (fromMaybe (Latitude 0) (jLatitude x))
                (fromMaybe (Longitude 0) (jLongitude x))
    in  map transformTo $ filter filterPos jannos

distanceOneToAll :: IndsWithPosition -> [IndsWithPosition] -> [(String, String, Double)]
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

distToWeight :: [Double] -> [Int]
distToWeight distances = 
    let closeness = map (1/) distances
    in map round $ rescale 0 100 closeness

rescale :: Double -> Double -> [Double] -> [Double]
rescale minNew maxNew xs = 
    let minOld = minimum xs
        maxOld = maximum xs
        a = (maxNew - minNew) / (maxOld - minOld)
        bs = map (\x -> x - maxOld) xs
    in map (\x -> a * x + maxNew) bs

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
