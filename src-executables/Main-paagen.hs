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
    -- 
    let startCoords = (51.510357, -0.116773)
    hPutStrLn stderr $ show (haversineDist (51.510357, -0.116773) (38.889931, -77.009003))

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