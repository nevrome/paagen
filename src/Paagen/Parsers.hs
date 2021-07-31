module Paagen.Parsers where

import Paagen.Types
import Paagen.Utils

import           Poseidon.Janno

import           Control.Exception              (throwIO)
import           Control.Monad                  (guard)
import qualified Text.Parsec                    as P
import qualified Text.Parsec.String             as P
import qualified Text.Parsec.Number             as P

readPopulationWithFractionString :: String -> Either String [PopulationWithFraction]
readPopulationWithFractionString s = case P.runParser populationWithFractionParser () "" s of
    Left p  -> Left (show p)
    Right x -> Right x

populationWithFractionParser :: P.Parser [PopulationWithFraction]
populationWithFractionParser = P.try (P.sepBy parsePopulationWithFraction (P.char ';' <* P.spaces))

parsePopulationWithFraction :: P.Parser PopulationWithFraction
parsePopulationWithFraction = do
    pop <- P.many (P.noneOf ",")
    _ <- P.oneOf ","
    frac <- read <$> P.many1 P.digit
    return (PopulationWithFraction pop frac)

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