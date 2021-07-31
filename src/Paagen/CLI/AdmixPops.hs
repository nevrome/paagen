module Paagen.CLI.AdmixPops where

import           Paagen.Parsers
import           Paagen.Types                   (PopulationWithFraction (..), GenoEntry (..))
import           Paagen.Utils

import           Control.Exception              (catch, throwIO, Exception)
import           Control.Monad                  (forM, guard)
import           Control.Monad.Random           (fromList, RandomGen, evalRand, newStdGen, getStdGen)
import           Data.List                      
import           Data.Maybe                     
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
import           System.FilePath                ((<.>), (</>))
import           System.IO                      (hPutStrLn, stderr, hPutStr)

data AdmixPopsOptions = AdmixPopsOptions {   
      _optBaseDirs :: [FilePath]
    , _optPopulationsWithFractions :: [PopulationWithFraction]
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
runAdmixPops (AdmixPopsOptions baseDirs popsWithFracs outFormat outDir) = do
    -- load Poseidon packages
    allPackages <- readPoseidonPackageCollection pacReadOpts baseDirs
    -- load janno tables
    let jannos = concatMap posPacJanno allPackages
    -- determine relevant packages and indices
    let pops = map pop popsWithFracs
    relevantPackages <- filterPackagesByPops pops allPackages
    indicesPerPop <- mapM (`extractIndsPerPop` relevantPackages) pops
    let popsFracsInds = zip3 pops (map frac popsWithFracs) indicesPerPop
    -- compile genotype data structure
    let [outInd, outSnp, outGeno] = case outFormat of 
            GenotypeFormatEigenstrat -> ["poi.ind", "poi.snp", "poi.geno"]
            GenotypeFormatPlink -> ["poi.fam", "poi.bim", "poi.bed"]
    putStrLn $ show popsFracsInds
    return ()
    -- compile genotype data
    -- runSafeT $ do
    --     (eigenstratIndEntries, eigenstratProd) <- getJointGenotypeData False False relevantPackages
    --     let [outG, outS, outI] = map (outDir </>) [outGeno, outSnp, outInd]
    --         newIndEntries = map (\x -> EigenstratIndEntry (ind x) Unknown (unit x)) pois
    --     let outConsumer = case outFormat of
    --             GenotypeFormatEigenstrat -> writeEigenstrat outG outS outI newIndEntries 
    --             GenotypeFormatPlink      -> writePlink      outG outS outI newIndEntries
    --     runEffect $ eigenstratProd >-> printSNPCopyProgress >-> P.mapM (sampleGenoForMultiplePOIs infoForIndividualPOIs) >-> outConsumer
    --     liftIO $ hPutStrLn stderr "Done"

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
