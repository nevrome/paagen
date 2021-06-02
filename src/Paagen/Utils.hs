module Paagen.Utils (
    PaagenException (..),
) where

import           Control.Exception (Exception)

data PaagenException =
    PaagenCLIParsingException String
    deriving (Show)

instance Exception PaagenException
