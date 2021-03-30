{-# LANGUAGE OverloadedStrings #-}

import           Paths_paagen      (version)

import           Control.Applicative    ((<|>))
import           Control.Exception      (catch)
import           Data.Version           (showVersion)
import qualified Options.Applicative    as OP
import           Poseidon.Package       (PoseidonPackage (..),
                                        getIndividuals,
                                        getJointGenotypeData,
                                        readPoseidonPackageCollection)
import           System.Exit            (exitFailure)
import           System.IO              (hPutStrLn, stderr)

import Data.Maybe (isJust, fromMaybe, catMaybes, fromJust)
import Poseidon.Janno
import Data.List (nub, tails, sortBy, intersect, maximumBy, group, sort, intercalate, elemIndex, elemIndices)
import qualified Data.Vector                as V
import           SequenceFormats.Eigenstrat (EigenstratIndEntry (..),
                                             EigenstratSnpEntry (..), GenoLine,
                                             writeEigenstrat, GenoEntry (..))
import Poseidon.GenotypeData
import qualified Pipes.Prelude              as P
import           Pipes.Safe                 (SafeT (..), runSafeT, throwM)
import Control.Monad (forM)
import           System.FilePath            ((<.>), (</>))
import Pipes
import           System.Console.ANSI        (hClearLine, hSetCursorColumn)
import Data.Function (on)

data TestOptions = TestOptions
    { _inTest :: String
    }

data Options = CmdTest TestOptions

main :: IO ()
main = do
    cmdOpts <- OP.customExecParser p optParserInfo
    runCmd cmdOpts
    --catch (runCmd cmdOpts) handler
    where
        p = OP.prefs OP.showHelpOnEmpty
        --handler :: PoseidonException -> IO ()
        --handler e = do
        --    hPutStrLn stderr $ renderPoseidonException e
        --    exitFailure

runCmd :: Options -> IO ()
runCmd o = case o of
    CmdTest opts -> runTest opts

optParserInfo :: OP.ParserInfo Options
optParserInfo = OP.info (OP.helper <*> versionOption <*> optParser) (
    OP.briefDesc <>
    OP.progDesc "paagen"
    )

versionOption :: OP.Parser (a -> a)
versionOption = OP.infoOption (showVersion version) (OP.long "version" <> OP.help "Show version")

optParser :: OP.Parser Options
optParser = OP.subparser (
        OP.command "test" testOptInfo <>
        OP.commandGroup "test:"
    )
  where
    testOptInfo = OP.info (OP.helper <*> (CmdTest <$> testOptParser))
        (OP.progDesc "test")

testOptParser :: OP.Parser TestOptions
testOptParser = TestOptions <$> parseTest

parseTest :: OP.Parser String 
parseTest = OP.strOption (
    OP.long "test" <> 
    OP.help "test" <>
    OP.value "test" <>
    OP.showDefault
    )

-- Actual program code

runTest :: TestOptions -> IO ()
runTest (TestOptions test) = do
    -- load Poseidon packages -- 
    allPackages <- readPoseidonPackageCollection True True False ["/home/clemens/test/fetchtest/already_ready"]
    -- load janno tables
    let jannos = concatMap posPacJanno allPackages
    -- transform to spatiotemporal positions
        stInds = jannoToSpaceTimePos jannos
    -- calculate distances
        positionOfInterest = IndsWithPosition "poi" $ SpatialTemporalPosition 1000 (Latitude 47.82) (Longitude 47.82)
        distancesToPoi = distanceOneToAll positionOfInterest stInds
    -- get X closest inds
        closest = getClosestInds 50 distancesToPoi
        closestIndividuals = map fst closest
        closestDistances = map snd closest
        closestWeights = distToWeight closestDistances
    putStrLn $ "Closest individuals: " ++ intercalate ", " closestIndividuals
    -- determine relevant packages
    relevantPackages <- filterPackagesByInds closestIndividuals allPackages
    closestIndices <- extractIndIndices closestIndividuals relevantPackages
    -- compile genotype data
    runSafeT $ do
        (eigenstratIndEntries, eigenstratProd) <- getJointGenotypeData False False relevantPackages
        let [outG, outS, outI] = map ("/home/clemens/test/paagentest" </>) ["huhu.geno", "huhu.snp", "huhu.ind"]
        let outConsumer = writeEigenstrat outG outS outI [EigenstratIndEntry "poi" Unknown "group_of_poi"]
        runEffect $ eigenstratProd >-> printSNPCopyProgress >-> P.map (mergeIndividuals closestIndices closestWeights) >-> outConsumer
        liftIO $ hClearLine stderr
        liftIO $ hSetCursorColumn stderr 0
        liftIO $ hPutStrLn stderr "SNPs processed: All done"

mergeIndividuals :: [Int] -> [Int] -> (EigenstratSnpEntry, GenoLine) -> (EigenstratSnpEntry, GenoLine)
mergeIndividuals individualIndices weights (snpEntry, genoLine) = 
    let relevantGenoEntries = [genoLine V.! i | i <- individualIndices]
        genoEntryIndices = getGenoIndices relevantGenoEntries
        weightsPerGenoEntry = sumWeights genoEntryIndices weights
        -- modeGenoEntry = mostCommon relevantGenoEntries
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
    let filterPos x = isJust (jDateBCADMedian x) && isJust (jLatitude x) && isJust (jLongitude x)
        transformTo x = IndsWithPosition (jIndividualID x) $
            SpatialTemporalPosition 
                (fromMaybe 0 (jDateBCADMedian x)) 
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
