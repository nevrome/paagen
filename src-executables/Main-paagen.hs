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

import Data.Maybe (isJust, fromMaybe, catMaybes)
import Poseidon.Janno
import Data.List (nub, tails, sortBy, intersect, maximumBy, group, sort)
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
        stInds = jannosToSTInds jannos
    -- calculate distances
        positionOfInterest = IndsWithPosition "poi" $ SpatialTemporalPosition 1000 (Latitude 47.82) (Longitude 47.82)
        distancesToPoi = distanceOneToAll positionOfInterest stInds
    -- get X closest inds
        closestInds = getXClosestInds 50 distancesToPoi
    -- determine relevant packages
    relevantPackages <- filterPackagesByInds closestInds allPackages
    indices <- extractIndIndices closestInds relevantPackages
    print indices
    -- compile genotype data
    runSafeT $ do
        (eigenstratIndEntries, eigenstratProd) <- getJointGenotypeData True True relevantPackages
        let eigenstratIndEntriesV = V.fromList eigenstratIndEntries
        let newEigenstratIndEntries = [eigenstratIndEntriesV V.! i | i <- indices]
        let [outG, outS, outI] = map ("/home/clemens/test/paagentest" </>) ["huhu.geno", "huhu.snp", "huhu.ind"]
        let outConsumer = writeEigenstrat outG outS outI newEigenstratIndEntries
        runEffect $ eigenstratProd >-> printSNPCopyProgress >-> P.map (mergeIndividuals indices) >-> outConsumer
        liftIO $ hClearLine stderr
        liftIO $ hSetCursorColumn stderr 0
        liftIO $ hPutStrLn stderr "SNPs processed: All done"
    -- 
    print closestInds

selectIndices :: [Int] -> (EigenstratSnpEntry, GenoLine) -> (EigenstratSnpEntry, GenoLine)
selectIndices indices (snpEntry, genoLine) = (snpEntry, V.fromList [genoLine V.! i | i <- indices])

mergeIndividuals :: [Int] -> (EigenstratSnpEntry, GenoLine) -> (EigenstratSnpEntry, GenoLine)
mergeIndividuals indices (snpEntry, genoLine) = 
    let relevantGenoEntries = [genoLine V.! i | i <- indices]
        modeGenoEntry = mostCommon relevantGenoEntries
    in (snpEntry, V.fromList [modeGenoEntry])

mostCommon :: Ord a => [a] -> a
mostCommon = head . maximumBy (compare `on` length) . group . sort

instance Ord GenoEntry where
    compare HomRef Het     = GT
    compare Het HomAlt     = GT
    compare HomAlt Missing = GT
    compare _ _            = EQ

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

filterPackagesByInds :: [String] -> [PoseidonPackage] -> IO [PoseidonPackage]
filterPackagesByInds indNamesStats packages = do
    fmap catMaybes . forM packages $ \pac -> do
        inds <- getIndividuals pac
        let indNamesPac   = [ind   | EigenstratIndEntry ind _ _     <- inds]
        if  length (intersect indNamesPac indNamesStats) > 0
        then return (Just pac)
        else return Nothing

getXClosestInds :: Int -> [(String, String, Double)] -> [String]
getXClosestInds n dists = map (\(_,x,_) -> x) $ take n $ sortBy (\(_,_,x) (_,_,y) -> compare x y) dists

distanceOneToAll :: IndsWithPosition -> [IndsWithPosition] -> [(String, String, Double)]
distanceOneToAll poi = map (distanceOneToOne poi)

distanceOneToOne :: IndsWithPosition -> IndsWithPosition -> (String, String, Double)
distanceOneToOne i1 i2 = (ind i1, ind i2, spatioTemporalDistance 1 (pos i1) (pos i2))

-- distanceForPair :: (IndsWithPosition, IndsWithPosition) -> (String, String, Double)
-- distanceForPair (i1, i2) = (ind i1, ind i2, spatioTemporalDistance 1 (pos i1) (pos i2))

-- pairs :: [a] -> [(a, a)]
-- pairs l = [(x,y) | (x:ys) <- tails l, y <- ys]

spatioTemporalDistance :: Double -> SpatialTemporalPosition -> SpatialTemporalPosition -> Double 
spatioTemporalDistance 
    scaling
    (SpatialTemporalPosition t1 (Latitude lat1) (Longitude lon1)) 
    (SpatialTemporalPosition t2 (Latitude lat2) (Longitude lon2)) =
    let tDist = fromIntegral (abs (t1 - t2)) * scaling
        sDist = haversineDist (lat1, lon1) (lat2, lon2)
    in sqrt $ tDist^2 + sDist^2

jannosToSTInds :: [JannoRow] -> [IndsWithPosition]
jannosToSTInds jannos = 
    let filterPos x = isJust (jDateBCADMedian x) && isJust (jLatitude x) && isJust (jLongitude x)
        transformTo x = IndsWithPosition (jIndividualID x) $
            SpatialTemporalPosition 
                (fromMaybe 0 (jDateBCADMedian x)) 
                (fromMaybe (Latitude 0) (jLatitude x))
                (fromMaybe (Longitude 0) (jLongitude x))
    in  map transformTo $ filter filterPos jannos

data IndsWithPosition = IndsWithPosition {
      ind :: String
    , pos :: SpatialTemporalPosition
} deriving (Show)

data SpatialTemporalPosition = SpatialTemporalPosition {
      time :: Int
    , lat :: Latitude 
    , lon :: Longitude
} deriving (Show)

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