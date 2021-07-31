module Paagen.Utils (
    PaagenException (..),
    extractFirst,
    zipGroup
) where

import           Control.Exception (Exception)

data PaagenException =
    PaagenCLIParsingException String
    deriving (Show)

instance Exception PaagenException


extractFirst :: (a, b, c) -> a
extractFirst (a,_,_) = a

zipGroup :: [a] -> [[b]] -> [(Int,a,b)]
zipGroup list nestedList =
    let lenghtsNestedList = map length nestedList
        listWithlenghtsNestedList = zip lenghtsNestedList list
        longerA = map (uncurry replicate) listWithlenghtsNestedList
    in zip3 [0..] (concat longerA) (concat nestedList)