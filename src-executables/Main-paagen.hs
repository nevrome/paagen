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

import Data.Maybe (isJust, fromMaybe)
import Poseidon.Janno
import Data.List (nub, tails, sortBy)


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
        closestInds = getXClosestInds 5 distancesToPoi
    -- get 
    print closestInds

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