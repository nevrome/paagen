module Paagen.CLI.AdmixPops where

import           Paagen.Parsers
import           Paagen.Types                   (PopulationWithFraction (..), GenoEntry (..))
import           Paagen.Utils

import           Control.Exception              (catch, throwIO, Exception)
import           Control.Monad                  (forM, guard, when)
import           Control.Monad.Random           (fromList, RandomGen, evalRand, newStdGen, getStdGen)
import           Data.List                      
import           Data.Maybe                     
import           Data.Ratio                     ((%))
import qualified Data.Vector                    as V
import           Pipes
import qualified Pipes.Prelude                  as P
import           Pipes.Safe                     (SafeT (..), runSafeT, throwM)
import           Poseidon.GenotypeData
import           Poseidon.Janno
import           Poseidon.Package
import           SequenceFormats.Eigenstrat     (EigenstratIndEntry (..),
                                                EigenstratSnpEntry (..), GenoLine,
                                                writeEigenstrat)
import           SequenceFormats.Plink          (writePlink)
import           System.Console.ANSI            (hClearLine, hSetCursorColumn)
import           System.Directory               (createDirectoryIfMissing)
import           System.FilePath                ((<.>), (</>))
import           System.IO                      (hPutStrLn, stderr, hPutStr)

data AdmixPopsOptions = AdmixPopsOptions {   
      _optBaseDirs :: [FilePath]
    , _optPopulationsWithFractions :: [[PopulationWithFraction]]
    , _optPopulationsWithFractionsFile :: Maybe FilePath
    , _optOutFormat :: GenotypeFormatSpec
    , _optOutPath :: FilePath
    }

pacReadOpts :: PackageReadOptions
pacReadOpts = defaultPackageReadOptions {
      _readOptVerbose          = True
    , _readOptStopOnDuplicates = True
    , _readOptIgnoreChecksums  = True --False
    , _readOptIgnoreGeno       = True --False
    , _readOptGenoCheck        = True
    }

runAdmixPops :: AdmixPopsOptions -> IO ()
runAdmixPops (AdmixPopsOptions baseDirs popsWithFracsDirect popsWithFracsFile outFormat outDir) = do
    -- check input
    popsWithFracsFromFile <- case popsWithFracsFile of
        Nothing -> return []
        Just f -> readPopulationWithFractionFromFile f
    let popsWithFracs = popsWithFracsDirect ++ popsWithFracsFromFile 
        pops = map (map pop) popsWithFracs
        fracs = map (map frac) popsWithFracs
    -- when (sum fracs /= 100) $ do
    --     throwIO $ PaagenCLIParsingException "Fractions have to sum to 100%"
    -- load Poseidon packages
    allPackages <- readPoseidonPackageCollection pacReadOpts baseDirs
    -- determine relevant packages and indices
    relevantPackages <- filterPackagesByPops (concat pops) allPackages
    indicesPerPop <- mapM (mapM (`extractIndsPerPop` relevantPackages)) pops
    let fracsInds = zipWith zip fracs indicesPerPop
        weights = map weightsPerInd fracsInds
    -- prepare weights per individual
    -- compile genotype data structure
    let [outInd, outSnp, outGeno] = case outFormat of 
            GenotypeFormatEigenstrat -> ["poi.ind", "poi.snp", "poi.geno"]
            GenotypeFormatPlink -> ["poi.fam", "poi.bim", "poi.bed"]
    putStrLn $ show $ zip3 pops fracs indicesPerPop
    -- create output directory
    createDirectoryIfMissing True outDir
    -- compile genotype data
    runSafeT $ do
        (eigenstratIndEntries, eigenstratProd) <- getJointGenotypeData False False relevantPackages
        let [outG, outS, outI] = map (outDir </>) [outGeno, outSnp, outInd]
            newIndEntry = EigenstratIndEntry "testInd" Unknown "testGroup"
        let outConsumer = case outFormat of
                GenotypeFormatEigenstrat -> writeEigenstrat outG outS outI [newIndEntry]
                GenotypeFormatPlink      -> writePlink      outG outS outI [newIndEntry]
        runEffect $ eigenstratProd >-> printSNPCopyProgress >-> P.mapM (sampleGenoForMultiplePOIs weights) >-> outConsumer
        liftIO $ hPutStrLn stderr "Done"

weightsPerInd :: [(Int, [Int])] -> ([Int], [Rational])
weightsPerInd fracAndInds = unzip $ concatMap
            (\(x,y) -> zipWith3 (\a b c -> (a, b % c)) 
                y 
                (replicate (length y) (toInteger x)) 
                (replicate (length y) (toInteger (length y)))
            ) fracAndInds

filterPackagesByPops :: [String] -> [PoseidonPackage] -> IO [PoseidonPackage]
filterPackagesByPops pops packages = do
    fmap catMaybes . forM packages $ \pac -> do
        inds <- getIndividuals pac
        let groupNamesPac = [group | EigenstratIndEntry _   _ group <- inds]
        if   not (null (groupNamesPac `intersect` pops))
        then return (Just pac)
        else return Nothing

extractIndsPerPop :: String -> [PoseidonPackage] -> IO [Int]
extractIndsPerPop pop relevantPackages = do
    let allPackageNames = map posPacTitle relevantPackages
    allIndEntries <- mapM getIndividuals relevantPackages
    let filterFunc (_ , pacName, EigenstratIndEntry ind _ group) = group == pop
    return $ map extractFirst $ filter filterFunc (zipGroup allPackageNames allIndEntries)
