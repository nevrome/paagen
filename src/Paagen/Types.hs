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
      admixtureId :: String
    , admixUnit :: String 
    , admixtureSet :: AdmixtureSet
}

data AdmixtureSet = AdmixtureSet {
    popFracList :: [PopulationWithFraction]
} 

data PopulationWithFraction = PopulationWithFraction {
      pop :: String
    , frac :: Int
} deriving (Show)

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
