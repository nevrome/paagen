module Paagen.CLI.SpaceTime where

import           Paagen.Parsers
import           Paagen.Types
import           Paagen.Utils

import           Control.Exception              (catch, throwIO, Exception)
import           Control.Monad                  (forM, guard)
import           Control.Monad.Random           (newStdGen, getStdGen)
import           Data.List                      (nub, tails, sortBy, intersect, maximumBy, group, sort, intercalate, elemIndex, elemIndices, unfoldr)
import           Data.Maybe                     (isJust, fromMaybe, catMaybes, fromJust)
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

data SpaceTimeOptions = SpaceTimeOptions {   
      _spatBaseDirs :: [FilePath]
    , _spatIndWithPosition :: [IndWithPosition]
    , _spatIndWithPositionFile :: Maybe FilePath
    , _spatNumberOfNearestNeighbors :: Int
    , _spatTemporalDistanceScaling :: Double
    , _spatOutFormat :: GenotypeFormatSpec
    , _spatOutPath :: FilePath
    }

pacReadOpts :: PackageReadOptions
pacReadOpts = defaultPackageReadOptions {
      _readOptVerbose          = True
    , _readOptStopOnDuplicates = True
    , _readOptIgnoreChecksums  = False
    , _readOptIgnoreGeno       = False
    , _readOptGenoCheck        = True
    }

runSpaceTime :: SpaceTimeOptions -> IO ()
runSpaceTime (SpaceTimeOptions baseDirs poisDirect poisFile numNeighbors temporalDistanceScaling outFormat outDir) = do
    -- compile pois
    poisFromFile <- case poisFile of
        Nothing -> return []
        Just f -> readIndWithPositionFromFile f
    let pois = poisDirect ++ poisFromFile
    -- load Poseidon packages
    allPackages <- readPoseidonPackageCollection pacReadOpts baseDirs
    -- load janno tables
    let jannos = concatMap posPacJanno allPackages
    -- transform to spatiotemporal positions
        stInds = jannoToSpaceTimePos jannos
    -- calculate distances
        distancesToPois = map (\x -> distanceOneToAll temporalDistanceScaling x stInds) pois 
    -- get X closest inds
        closest = map (getClosestInds numNeighbors) distancesToPois
        closestIndividuals = map (map fst) closest
        closestDistances = map (map snd) closest
        closestWeights = map distToWeight closestDistances
    -- putStrLn $ "Closest individuals: " ++ intercalate ", " closestIndividuals
    -- determine relevant packages
    relevantPackages <- filterPackagesByInds (nub $ concat closestIndividuals) allPackages
    closestIndices <- mapM (`extractIndIndices` relevantPackages) closestIndividuals
    -- merge poi indices and weights
    let infoForIndividualPOIs = zip closestIndices closestWeights
    -- compile genotype data structure
    let [outInd, outSnp, outGeno] = case outFormat of
            GenotypeFormatEigenstrat -> ["res.ind", "res.snp", "res.geno"]
            GenotypeFormatPlink -> ["res.fam", "res.bim", "res.bed"]
    -- compile genotype data
    runSafeT $ do
        (eigenstratIndEntries, eigenstratProd) <- getJointGenotypeData False False relevantPackages
        let [outG, outS, outI] = map (outDir </>) [outGeno, outSnp, outInd]
            newIndEntries = map (\x -> EigenstratIndEntry (spatInd x) Unknown (spatUnit x)) pois
        let outConsumer = case outFormat of
                GenotypeFormatEigenstrat -> writeEigenstrat outG outS outI newIndEntries 
                GenotypeFormatPlink      -> writePlink      outG outS outI newIndEntries
        runEffect $ eigenstratProd >-> printSNPCopyProgress >-> P.mapM (sampleGenoForMultiplePOIs infoForIndividualPOIs) >-> outConsumer
        liftIO $ hPutStrLn stderr "Done"

jannoToSpaceTimePos :: [JannoRow] -> [IndWithPosition]
jannoToSpaceTimePos jannos =
    let filterPos x = (isJust (jDateBCADMedian x) || jDateType x == Just Modern) && isJust (jLatitude x) && isJust (jLongitude x)
        transformTo x = IndWithPosition (jIndividualID x) (head $ getJannoList $ jGroupName x) $
            SpatialTemporalPosition 
                (fromMaybe 2000 (jDateBCADMedian x)) -- If empty then modern (after filter), which is approx. 2000AD
                (fromMaybe (Latitude 0) (jLatitude x))
                (fromMaybe (Longitude 0) (jLongitude x))
    in  map transformTo $ filter filterPos jannos

distanceOneToAll :: Double -> IndWithPosition -> [IndWithPosition] -> [(String, String, Double)]
distanceOneToAll temporalDistanceScaling poi = 
    map (distanceOneToOne poi) 
    where distanceOneToOne i1 i2 = (spatInd i1, spatInd i2, spatioTemporalDistance temporalDistanceScaling (spatPos i1) (spatPos i2))

spatioTemporalDistance :: Double -> SpatialTemporalPosition -> SpatialTemporalPosition -> Double 
spatioTemporalDistance 
    scaling
    (SpatialTemporalPosition t1 (Latitude lat1) (Longitude lon1)) 
    (SpatialTemporalPosition t2 (Latitude lat2) (Longitude lon2)) =
    let tDist = fromIntegral (abs (t1 - t2)) * scaling
        sDist = haversineDist (lat1, lon1) (lat2, lon2)
    in sqrt $ tDist^2 + sDist^2

haversineDist :: (Double, Double) -> (Double, Double) -> Double 
haversineDist (lat1, lon1) (lat2, lon2) =
    let r = 6371000  -- radius of Earth in meters
        toRadians n = n * pi / 180
        square x = x * x
        cosr = cos . toRadians
        dlat = toRadians (lat1 - lat2) / 2
        dlon = toRadians (lon1 - lon2) / 2
        a = square (sin dlat) + cosr lat1 * cosr lat2 * square (sin dlon)
        c = 2 * atan2 (sqrt a) (sqrt (1 - a))
    in (r * c) / 1000

getClosestInds :: Int -> [(String, String, Double)] -> [(String, Double)]
getClosestInds n dists = map (\(_,x,y) -> (x,y)) $ take n $ sortBy (\(_,_,x) (_,_,y) -> compare x y) dists

distToWeight :: [Double] -> [Rational]
distToWeight distances = 
    let closeness = map (1/) distances
    in rescale 0 100 closeness

rescale :: Double -> Double -> [Double] -> [Rational]
rescale minNew maxNew xs = 
    let minOld = minimum xs
        maxOld = maximum xs
        a = (maxNew - minNew) / (maxOld - minOld)
        bs = map (\x -> x - maxOld) xs
        res = map (\x -> a * x + maxNew) bs
    in map toRational res

filterPackagesByInds :: [String] -> [PoseidonPackage] -> IO [PoseidonPackage]
filterPackagesByInds indNamesStats packages = do
    fmap catMaybes . forM packages $ \pac -> do
        inds <- getIndividuals pac
        let indNamesPac   = [ind   | EigenstratIndEntry ind _ _     <- inds]
        if  length (intersect indNamesPac indNamesStats) > 0
        then return (Just pac)
        else return Nothing

extractIndIndices :: [String] -> [PoseidonPackage] -> IO [Int]
extractIndIndices indNames relevantPackages = do
    let allPackageNames = map posPacTitle relevantPackages
    allIndEntries <- mapM getIndividuals relevantPackages
    let filterFunc (_ , pacName, EigenstratIndEntry ind _ group) = ind `elem` indNames
    return $ map extractFirst $ filter filterFunc (zipGroup allPackageNames allIndEntries)
