module Paagen.Types (
    IndWithAdmixtureSet (..),
    AdmixtureSet (..),
    PopulationWithFraction (..),
    IndWithPosition (..),
    SpatialTemporalPosition (..),
    GenoEntry (..)
) where

import           Poseidon.Janno
import           Poseidon.Package
import           Poseidon.GenotypeData

import           SequenceFormats.Eigenstrat     (GenoEntry (..))

data IndWithAdmixtureSet = IndWithAdmixtureSet {
      admixInd :: String
    , admixUnit :: String 
    , admixSet :: AdmixtureSet
} deriving (Show)

data AdmixtureSet = AdmixtureSet {
    popFracList :: [PopulationWithFraction]
} deriving (Show)

data PopulationWithFraction = PopulationWithFraction {
      pop :: String
    , frac :: Int
} deriving (Show)

data IndWithPosition = IndWithPosition {
      spatInd :: String
    , spatUnit :: String
    , spatPos :: SpatialTemporalPosition
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
