module Paagen.Utils where

import                  Paagen.Types 

import                  Control.Exception               (Exception)
import                  Control.Monad.Random            (fromList, RandomGen, evalRand, newStdGen, getStdGen)
import                  Data.List
import                  Data.Ratio                      ((%))   
import                  SequenceFormats.Eigenstrat
import                  Pipes
import                  Pipes.Safe
import qualified        Data.Vector as V

data PaagenException =
    PaagenCLIParsingException String
    deriving (Show)

instance Exception PaagenException

sampleGenoForMultipleIndWithAdmixtureSet :: [[([Int], Rational)]] -> (EigenstratSnpEntry, GenoLine) -> SafeT IO (EigenstratSnpEntry, GenoLine)
sampleGenoForMultipleIndWithAdmixtureSet infoForIndividualInd (snpEntry, genoLine) = do
    entries <- mapM (`sampleGenoForOneIndWithAdmixtureSet` genoLine) infoForIndividualInd
    return (snpEntry, V.fromList entries)

sampleGenoForOneIndWithAdmixtureSet :: [([Int], Rational)] -> GenoLine -> SafeT IO GenoEntry
sampleGenoForOneIndWithAdmixtureSet xs genoLine = do
    gen <- liftIO getStdGen
    let sampledAllelesPerPop = map (\(x,y) -> (sampleWeightedList gen $ getAlleleFrequencyInPopulation x genoLine, y)) xs
        sampledAlleleAcrossPops = sampleWeightedList gen sampledAllelesPerPop
    liftIO newStdGen
    return sampledAlleleAcrossPops

getAlleleFrequencyInPopulation :: [Int] -> GenoLine -> [(GenoEntry, Rational)]
getAlleleFrequencyInPopulation individualIndices genoLine =
    let relevantGenoEntries = [genoLine V.! i | i <- individualIndices]
    in  if all (Missing ==) relevantGenoEntries
        then [(Missing, 1)]
        else calcFractions $ filter (Missing /=) relevantGenoEntries

calcFractions :: Ord a => [a] -> [(a, Rational)]
calcFractions xs =
    let ls = toInteger $ length xs
    in map (\x -> (head x, toInteger (length x) % ls)) $ group $ sort xs

extractFirst :: (a, b, c) -> a
extractFirst (a,_,_) = a

zipGroup :: [a] -> [[b]] -> [(Int,a,b)]
zipGroup list nestedList =
    let lenghtsNestedList = map length nestedList
        listWithlenghtsNestedList = zip lenghtsNestedList list
        longerA = map (uncurry replicate) listWithlenghtsNestedList
    in zip3 [0..] (concat longerA) (concat nestedList)

sampleGenoForMultiplePOIs :: [([Int], [Rational])] -> (EigenstratSnpEntry, GenoLine) -> SafeT IO (EigenstratSnpEntry, GenoLine)
sampleGenoForMultiplePOIs infoForIndividualPOIs (snpEntry, genoLine) = do
    entries <- mapM (\(x,y) -> sampleGenoForOnePOI x y genoLine) infoForIndividualPOIs
    return (snpEntry, V.fromList entries)

sampleGenoForOnePOI :: [Int] -> [Rational] -> GenoLine -> SafeT IO GenoEntry
sampleGenoForOnePOI individualIndices weights genoLine = do
    let relevantGenoEntries = [genoLine V.! i | i <- individualIndices]
        -- count occurrence of GenoEntries
        genoEntryIndices = getGenoIndices relevantGenoEntries
        -- sum distance-based weight for each GenoEntry
        weightsPerGenoEntry = sumWeights genoEntryIndices weights
    -- sample GenoEntry based on weight
    gen <- liftIO getStdGen
    -- liftIO $ hPutStrLn stderr (show gen)
    let selectedGenoEntry = sampleWeightedList gen weightsPerGenoEntry
    liftIO newStdGen
    -- return 
    return selectedGenoEntry

sampleWeightedList :: RandomGen g => g -> [(a, Rational)] -> a
sampleWeightedList gen weights = head $ evalRand m gen
    where m = sequence . repeat . fromList $ weights

getGenoIndices :: Eq a => [a] -> [(a, [Int])]
getGenoIndices xs = 
    let unique = nub xs
        indices = map (\v -> elemIndices v xs) unique
    in  zip unique indices

sumWeights :: Num b => [(a, [Int])] -> [b] -> [(a, b)]
sumWeights xs weights = map (\(x, ys) -> (x, sum $ subset ys weights)) xs
    where
        subset :: [Int] -> [a] -> [a]
        subset indices xs = [xs !! i | i <- indices]