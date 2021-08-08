{-# LANGUAGE OverloadedStrings #-}

import           Paths_paagen                   (version)
import           Paagen.Parsers
import           Paagen.Types
import           Paagen.CLI.SpaceTime           (runSpaceTime, 
                                                 SpaceTimeOptions (..))
import           Paagen.CLI.AdmixPops           (runAdmixPops,
                                                 AdmixPopsOptions (..))
import           Poseidon.GenotypeData

import           Control.Applicative            ((<|>))
import           Data.Char                      (isSpace)
import           Data.Function                  (on)
import           Data.Version                   (showVersion)
import qualified Options.Applicative            as OP
import           System.Exit                    (exitFailure)
import Control.Exception (SomeException(SomeException))

-- data types

data Options = 
      CmdSpaceTime SpaceTimeOptions
    | CmdAdmixPops AdmixPopsOptions

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
    CmdAdmixPops opts -> runAdmixPops opts

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
        OP.command "spacetime" spaceTimeOptInfo <>
        OP.command "admixpops" admixPopsOptInfo
    )
  where
    spaceTimeOptInfo = OP.info (OP.helper <*> (CmdSpaceTime <$> spaceTimeOptParser))
        (OP.progDesc "Genotype profile generation based on spatiotemporal position")
    admixPopsOptInfo = OP.info (OP.helper <*> (CmdAdmixPops <$> admixPopsOptParser))
        (OP.progDesc "Genotype profile generation based on admixture proportions")

spaceTimeOptParser :: OP.Parser SpaceTimeOptions
spaceTimeOptParser = SpaceTimeOptions <$> parseBasePaths
                                      <*> parseIndWithPositionDirect
                                      <*> parseIndWithPositionFromFile
                                      <*> parseNumberOfNearestNeighbors
                                      <*> parseTemporalDistanceScaling
                                      <*> parseOutGenotypeFormat
                                      <*> parseOutPath

admixPopsOptParser = AdmixPopsOptions <$> parseBasePaths
                                      <*> parseIndWithAdmixtureSetDirect
                                      <*> parseIndWithAdmixtureSetFromFile
                                      <*> parseOutGenotypeFormat
                                      <*> parseOutPath

parseIndWithAdmixtureSetDirect :: OP.Parser [IndWithAdmixtureSet]
parseIndWithAdmixtureSetDirect = OP.option (OP.eitherReader readIndWithAdmixtureSetString) (
    OP.long "admixString" <>
    OP.short 'a' <>
    OP.value [] <>
    OP.help "Population setup of interest: Each setup is a string of the form \
            \\"[id:group](population1=10+population2=30+...)\". Multiple setups can be listed separated by ;. \
            \id and group are simple strings. \
            \The population fractions must be simple integers and sum to 100."
    )

parseIndWithAdmixtureSetFromFile :: OP.Parser (Maybe FilePath)
parseIndWithAdmixtureSetFromFile = OP.option (Just <$> OP.str) (OP.long "admixFile" <>
    OP.value Nothing <>
    OP.help "A file with a list of spatiotemporal positions. \
            \Works just as -p, but multiple values can be given separated by newline. \
            \-a and --admixFile can be combined."
    )

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
    OP.help "Spatiotemporal positions of interest: Each position is a string of the form \
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

parseNumberOfNearestNeighbors :: OP.Parser Int
parseNumberOfNearestNeighbors = OP.option OP.auto (
    OP.long "neighbors" <> 
    OP.help "Number of nearest neighbors to consider for the calculation" <>
    OP.value 50 <>
    OP.showDefault
    )

parseTemporalDistanceScaling :: OP.Parser Double
parseTemporalDistanceScaling = OP.option OP.auto (
    OP.long "temporalDistanceScaling" <> 
    OP.help "Scaling factor for temporal distance in relation to spatial distance. \
            \Default is 1, so 1year = 1km." <>
    OP.value 1 <>
    OP.showDefault
    )

parseOutGenotypeFormat :: OP.Parser GenotypeFormatSpec
parseOutGenotypeFormat = OP.option (OP.eitherReader readGenotypeFormat) (
    OP.long "outFormat" <>
    OP.help "The format of the output genotype data: EIGENSTRAT or PLINK" <>
    OP.value GenotypeFormatEigenstrat
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
