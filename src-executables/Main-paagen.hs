{-# LANGUAGE OverloadedStrings #-}

import           Paths_paagen      (version)

import           Control.Applicative    ((<|>))
import           Control.Exception      (catch)
import           Data.Version           (showVersion)
import qualified Options.Applicative    as OP
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

runTest :: TestOptions -> IO ()
runTest (TestOptions test) = do
    hPutStrLn stderr $ "test"

