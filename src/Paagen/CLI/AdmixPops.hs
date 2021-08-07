module Paagen.CLI.AdmixPops where

import           Paagen.Parsers
import           Paagen.SampleGeno
import           Paagen.Types
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
import           Poseidon.BibFile
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
      _admixBaseDirs :: [FilePath]
    , _admixIndWithAdmixtureSet :: [IndWithAdmixtureSet]
    , _admixIndWithAdmixtureSetFile :: Maybe FilePath
    , _admixOutFormat :: GenotypeFormatSpec
    , _admixOutPath :: FilePath
    }

pacReadOpts :: PackageReadOptions
pacReadOpts = defaultPackageReadOptions {
      _readOptVerbose          = True
    , _readOptStopOnDuplicates = True
    , _readOptIgnoreChecksums  = False
    , _readOptIgnoreGeno       = False
    , _readOptGenoCheck        = True
    }

runAdmixPops :: AdmixPopsOptions -> IO ()
runAdmixPops (AdmixPopsOptions baseDirs popsWithFracsDirect popsWithFracsFile outFormat outDir) = do
    -- compile individuals
    popsWithFracsFromFile <- case popsWithFracsFile of
        Nothing -> return []
        Just f -> readIndWithAdmixtureSetFromFile f
    let requestedInds = popsWithFracsDirect ++ popsWithFracsFromFile 
        popsWithFracs = map (popFracList . admixSet) requestedInds
        pops = map (map pop) popsWithFracs
    print requestedInds
    -- when (sum fracs /= 100) $ do
    --     throwIO $ PaagenCLIParsingException "Fractions have to sum to 100%"
    -- load Poseidon packages
    allPackages <- readPoseidonPackageCollection pacReadOpts baseDirs
    -- determine relevant packages and indices
    relevantPackages <- filterPackagesByPops (concat pops) allPackages
    popsFracsInds <- mapM (mapM (`extractIndsPerPop` relevantPackages)) popsWithFracs
    -- compile genotype data structure
    let [outInd, outSnp, outGeno] = case outFormat of
            GenotypeFormatEigenstrat -> ["admixpops_package.ind", "admixpops_package.snp", "admixpops_package.geno"]
            GenotypeFormatPlink -> ["admixpops_package.fam", "admixpops_package.bim", "admixpops_package.bed"]
    -- create output directory
    createDirectoryIfMissing True outDir
    -- compile genotype data
    hPutStrLn stderr "Compiling chimeras"
    runSafeT $ do
        (eigenstratIndEntries, eigenstratProd) <- getJointGenotypeData False False relevantPackages
        let [outG, outS, outI] = map (outDir </>) [outGeno, outSnp, outInd]
            newIndEntries = map (\x -> EigenstratIndEntry (admixInd x) Unknown (admixUnit x)) requestedInds
        let outConsumer = case outFormat of
                GenotypeFormatEigenstrat -> writeEigenstrat outG outS outI newIndEntries
                GenotypeFormatPlink      -> writePlink      outG outS outI newIndEntries
        runEffect $ eigenstratProd >-> printSNPCopyProgress >-> P.mapM (sampleGenoForMultipleIndWithAdmixtureSet popsFracsInds) >-> outConsumer
        liftIO $ hPutStrLn stderr "Done"
    -- complete poseidon package
    hPutStrLn stderr "Wrapping generated genotype data in a Poseidon package"
    let genotypeData = GenotypeDataSpec outFormat outGeno Nothing outSnp Nothing outInd Nothing Nothing
    inds <- loadIndividuals outDir genotypeData
    pac <- newPackageTemplate outDir "admixpops_package" genotypeData (Just inds) Nothing Nothing
    writePoseidonPackage pac
    writeJannoFile (outDir </> "admixpops_package" <.> "janno") $ posPacJanno pac
    writeBibTeXFile (outDir </> "admixpops_package" <.> "bib") $ posPacBib pac

filterPackagesByPops :: [String] -> [PoseidonPackage] -> IO [PoseidonPackage]
filterPackagesByPops pops packages = do
    fmap catMaybes . forM packages $ \pac -> do
        inds <- getIndividuals pac
        let groupNamesPac = [group | EigenstratIndEntry _   _ group <- inds]
        if   not (null (groupNamesPac `intersect` pops))
        then return (Just pac)
        else return Nothing

extractIndsPerPop :: PopulationWithFraction -> [PoseidonPackage] -> IO ([Int], Rational)
extractIndsPerPop (PopulationWithFraction pop frac) relevantPackages = do
    let allPackageNames = map posPacTitle relevantPackages
    allIndEntries <- mapM getIndividuals relevantPackages
    let filterFunc (_ , pacName, EigenstratIndEntry ind _ group) = group == pop
    return (map extractFirst $ filter filterFunc (zipGroup allPackageNames allIndEntries), toInteger frac % 100)
