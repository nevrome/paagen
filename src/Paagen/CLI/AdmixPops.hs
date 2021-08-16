module Paagen.CLI.AdmixPops where

import           Paagen.Parsers
import           Paagen.SampleGeno
import           Paagen.Types
import           Paagen.Utils

import           Control.Monad                  (forM, when, unless)
import           Data.List                      
import           Data.Maybe
import           Data.Ratio                     ((%))                     
import           Pipes
import qualified Pipes.Prelude                  as P
import           Pipes.Safe                     (runSafeT)
import           Poseidon.GenotypeData
import           Poseidon.Janno
import           Poseidon.Package
import           SequenceFormats.Eigenstrat     (EigenstratIndEntry (..), writeEigenstrat)
import           SequenceFormats.Plink          (writePlink)
import           System.Directory               (createDirectoryIfMissing)
import           System.FilePath                ((</>))
import           System.IO                      (hPutStrLn, stderr, hPrint)
import Control.Exception (throwIO)

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
    hPutStrLn stderr "Checking chimeras"
    mapM_ (hPrint stderr) (take 5 requestedInds)
    when (length requestedInds > 5) $ do
        hPutStrLn stderr "..."
    -- validating input
    let individualsGrouped = filter (\x -> length x > 1) $ group $ sort $ map admixInd requestedInds
    unless (null individualsGrouped) $ do
                throwIO $ PaagenCLIParsingException $
                    "Duplicate individual names: " ++ intercalate "," (nub $ concat individualsGrouped)
    mapM_ checkIndsWithAdmixtureSets requestedInds
    -- load Poseidon packages
    allPackages <- readPoseidonPackageCollection pacReadOpts baseDirs
    -- determine relevant packages and indices
    relevantPackages <- filterPackagesByPops (concat pops) allPackages
    popsFracsInds <- mapM (mapM (`extractIndsPerPop` relevantPackages)) popsWithFracs
    -- compile genotype data structure
    let [outInd, outSnp, outGeno] = case outFormat of
            GenotypeFormatEigenstrat -> ["admixpops_package.ind", "admixpops_package.snp", "admixpops_package.geno"]
            GenotypeFormatPlink -> ["admixpops_package.fam", "admixpops_package.bim", "admixpops_package.bed"]
    -- create output poseidon package
    hPutStrLn stderr "Creating output Poseidon package"
    createDirectoryIfMissing True outDir
    let genotypeData = GenotypeDataSpec outFormat outGeno Nothing outSnp Nothing outInd Nothing Nothing
        pac = newMinimalPackageTemplate outDir "admixpops_package" genotypeData
    writePoseidonPackage pac
    -- compile genotype data
    hPutStrLn stderr "Compiling chimeras"
    runSafeT $ do
        (_, eigenstratProd) <- getJointGenotypeData False False relevantPackages
        let [outG, outS, outI] = map (outDir </>) [outGeno, outSnp, outInd]
            newIndEntries = map (\x -> EigenstratIndEntry (admixInd x) Unknown (admixUnit x)) requestedInds
        let outConsumer = case outFormat of
                GenotypeFormatEigenstrat -> writeEigenstrat outG outS outI newIndEntries
                GenotypeFormatPlink      -> writePlink      outG outS outI newIndEntries
        runEffect $ eigenstratProd >-> printSNPCopyProgress >-> P.mapM (sampleGenoForMultipleIndWithAdmixtureSet popsFracsInds) >-> outConsumer
        liftIO $ hPutStrLn stderr "Done"

checkIndsWithAdmixtureSets :: IndWithAdmixtureSet -> IO ()
checkIndsWithAdmixtureSets cur@(IndWithAdmixtureSet _ _ (AdmixtureSet _popFracList)) = do
    checkPopFracList _popFracList
    where
        checkPopFracList :: [PopulationWithFraction] -> IO ()
        checkPopFracList xs = do
            let fracs = map frac xs
            when (sum fracs /= 100) $ do
                throwIO $ PaagenCLIParsingException $
                    "Fractions in " ++ show cur ++ " do not to sum to 100%"

filterPackagesByPops :: [String] -> [PoseidonPackage] -> IO [PoseidonPackage]
filterPackagesByPops pops packages = do
    fmap catMaybes . forM packages $ \pac -> do
        inds <- getIndividuals pac
        let groupNamesPac = [groupName | EigenstratIndEntry _ _ groupName <- inds]
        if   not (null (groupNamesPac `intersect` pops))
        then return (Just pac)
        else return Nothing

extractIndsPerPop :: PopulationWithFraction -> [PoseidonPackage] -> IO ([Int], Rational)
extractIndsPerPop (PopulationWithFraction _pop _frac) relevantPackages = do
    let allPackageNames = map posPacTitle relevantPackages
    allIndEntries <- mapM getIndividuals relevantPackages
    let filterFunc (_,_,EigenstratIndEntry _ _ _group) = _group == _pop
    return (map extractFirst $ filter filterFunc (zipGroup allPackageNames allIndEntries), toInteger _frac % 100)
