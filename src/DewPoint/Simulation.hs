module DewPoint.Simulation (integrate, paramAdvectionDiff) where

import Data.Time.Clock
import Data.Array.Repa ( (!), (+^), (-^) )
import qualified Data.Array.Repa as R
import qualified DewPoint.Grid as G
import qualified DewPoint.Physics as Ph
--import Control.Exception.Base
import DewPoint.Geo(Lat)

import Debug.Trace
import Data.Array.Repa.Repr.Unboxed(Unbox)
import Data.Array.Repa.Eval(Elt)
import Data.Array.Repa.Algorithms.Randomish(randomishIntArray)
import Data.Time.Clock.POSIX(utcTimeToPOSIXSeconds)

-- | temperature increase due to solar influx at given cell at given time interval in seconds
tempDiffSolarInflux :: (Ord a, Floating a, Unbox a) =>
                       G.Grid
                    -> R.Array R.U G.Ix a -- ^ cloud cover
                    -> UTCTime
                    -> a         -- ^ time interval in seconds
                    -> G.Ix
                    -> a
tempDiffSolarInflux g gridCloudCover t dt ix =
    dt * Ph.temperatureDiffPerSecond ( (G.gridCellGeoCoord g) ! ix ) t ( gridCloudCover ! ix )


-- | temperature change due to emission, always negative
tempDiffEmission:: (Ord a, Floating a, Elt a, Unbox a) =>
                   R.Array R.U G.Ix a -- ^ air temperature
                -> R.Array R.U G.Ix a -- ^ temperature increase due to solar influx and advection
                -> R.Array R.D G.Ix a -- ^ decrease of temperature due to emission
tempDiffEmission gridAirTemp tempAbsorbed =
    let minTemp = Ph.toKelvin (-70)
        adjustTemp te' = ( max (te' - minTemp) 0 ) ** 4
        tempE = R.map adjustTemp gridAirTemp
        tempEsum = R.sumAllS tempE
        tempCellEmissionCoeff = R.map (\v -> -v / tempEsum) tempE
    in R.map ( * R.sumAllS tempAbsorbed ) tempCellEmissionCoeff


-- | parameter advection due to winds
--
--     -----> U
--
--       dx
--  V  -------
--     |_____|_______  _
--  ^  |     |       |  |
--  |  | cXY |  cX   |  |  dy
--  |  |_____|_______| _|
--     |     |       |
--     | cY  |  c0   |
--     |     |       |
--     ---------------
--
--
paramAdvectionDiff :: (Floating a, Ord a) =>
                      a -- ^ U-wind along X (m/s)
                   -> a -- ^ V-wind along Y (m/s)
                   -> a -- ^ x-source cell quantity
                   -> a -- ^ y-source cell quantity
                   -> a -- ^ length of cell side (assuming cell is always a square)
                   -> a -- ^ time interval (sec)
                   -> a -- ^ initial quantity
                   -> (a, a, a) -- ^ change of quantity in the cell + influx of quantity  (adX, adY)
paramAdvectionDiff u v aXs aYs l dt a0 =
  let dx = min (abs(u) * dt) l
      dy = min (abs(v) * dt) l
      totalArea = l**2
      c0 = (l - dy) * (l - dx) / totalArea
      cX = (l - dy) * dx / totalArea
      cY = (l - dx) * dy / totalArea
      cXY = 1 - c0 - cX - cY
      cX' = if (dx + dy > 0) then dx / (dx + dy) else 0
      cY' = if (dx + dy > 0) then 1 - cX' else 0
      adX = aXs * ( cX + cX' * cXY ) -- ^ loss quantity at source cell X
      adY = aYs * ( cY + cY' * cXY ) -- ^ loss quantity at source cell Y
      a1 = c0 * a0 + adX + adY      -- ^ new temperature in the cell
  in (a1 - a0, adX, adY)


tempDiffAdvection :: (Floating a, Ord a, Unbox a) =>
                     G.Grid
                  -> R.Array R.U G.Ix a -- air temp grid
                  -> R.Array R.U G.Ix a -- u-wind grid
                  -> R.Array R.U G.Ix a -- v-wind grid
                  -> a
                  -> G.Ix
                  -> a
tempDiffAdvection g gridAirTemp gridUWind gridVWind dt ix
    = let u = gridUWind ! ix
          v = gridVWind ! ix
          t0 = gridAirTemp ! ix
          ixU = G.gridSrcCellX g ix u
          ixV = G.gridSrcCellY g ix v
          tXs = gridAirTemp ! ixU
          tYs = gridAirTemp ! ixV
          l = G.gridCellLengthMeters g
          (dT, _, _) = paramAdvectionDiff u v tXs tYs l dt t0
      in dT

diffTimeToSec :: Floating a => DiffTime -> a
diffTimeToSec = fromRational . toRational

-- | second order central differences
centralDifferences :: Floating a => a -> a -> a -> a
centralDifferences p n step = (n - p) / (2 * step)

-- | forward difference for the first element
forwardDifferences :: Floating a => a -> a -> a -> a -> a
forwardDifferences e0 e1 e2 step = -(3.0*e0 - 4.0*e1 + e2) / (2.0*step)

-- | backward difference for the last element
backwardDifferences :: Floating a => a -> a -> a -> a -> a
backwardDifferences l1 l2 l3 = negate . forwardDifferences l1 l2 l3


-- | calculate gradient value by X coordinate at given grid cell
calculateGrad :: (Ord a, Floating a, Unbox a) =>
                 (G.Ix -> a -> G.Ix) -- ^ function to get index of src cell given vector value
              -> R.Array R.U G.Ix a -- ^ 2D array to calculate gradient for
              -> a -- ^ gradient step
              -> G.Ix   -- ^ cell index
              -> a -- ^ gradient value
calculateGrad getSrc ar s ix =
    let pix1 = getSrc ix 1
        nix1 = getSrc ix (-1)
        val = R.index ar
        adaptiveDiff a_1 a a1
            | a_1 == a = forwardDifferences (val a) (val a1) (val ( getSrc a1 (-1) ) ) s
            | a1  == a = backwardDifferences (val a) (val a_1) (val ( getSrc a_1 1 ) )  s
            | otherwise  = centralDifferences (val a_1) (val a1) s
    in adaptiveDiff pix1 ix nix1

-- | calculate wind component
calculateWind :: (Ord a, Floating a, Unbox a) =>
                 G.Grid
              -> R.Array R.U G.Ix a     -- ^ 500Mb heights
              -> (Lat -> a -> a)      -- ^ calculate wind at given latitude and gradient value
              -> (G.Ix -> a -> G.Ix)  -- ^ function to get index src cell given vector value
              -> G.Ix                 -- ^ cell index
              -> a
calculateWind g ar calcWindFn getSrcFn ix =
    let (lat,_) = (G.gridCellGeoCoord g) ! ix
    in calcWindFn lat $ (calculateGrad getSrcFn) ar (G.gridCellLengthMeters g) ix

-- | calculate U-wind component (zonal wind)
calculateUWind :: (Ord a, Floating a, Unbox a) => G.Grid -> R.Array R.U G.Ix a -> G.Ix-> a
calculateUWind g ar ix = calculateWind g ar Ph.geostrophicZonalWind (G.gridSrcCellY g) ix

-- | calculate V-wind component (meridional wind)
calculateVWind :: (Ord a, Floating a, Unbox a) => G.Grid -> R.Array R.U G.Ix a -> G.Ix -> a
calculateVWind g ar ix = calculateWind g ar Ph.geostrophicMeridionalWind (G.gridSrcCellX g) ix

-- | calculate pressure at 500Mb heights
calculate500MbHeights :: (Floating a, Unbox a) =>
                         R.Array R.U G.Ix a  -- ^ temperatures
                      -> G.Ix                -- ^ cell index
                      -> a
calculate500MbHeights tempAr = Ph.calculatePressureHeight 0.5 . R.index tempAr

-- | diff in water density
calculateWaterVaporDiff :: (Ord a, Floating a, Unbox a) =>
                           R.Array R.U G.Ix a
                        -> R.Array R.U G.Ix a
                        -> a
                        -> G.Ix
                        -> a
calculateWaterVaporDiff gridAirTemp gridEvapCoeff dts ix =
    dts * Ph.surfaceEvaporation (gridAirTemp ! ix) (gridEvapCoeff ! ix)


-- | calculate water vapor condensation
calculateWaterCondDiff :: (Ord a, Floating a, Unbox a) =>
                           R.Array R.U G.Ix a -- ^ water vapor density
                        -> R.Array R.U G.Ix a -- ^ air temperature
                        -> R.Array R.U G.Ix a -- ^ surface height
                        -> a                  -- ^ time interval in sec
                        -> G.Ix
                        -> a
calculateWaterCondDiff gridWaterVapor gridAirTemp gridSurfaceHeight dt ix  =
    if temp <= dewPointT
    then min water (water * Ph.condensationRate * dt)
    else 0
        where rh = Ph.relativeHumidity water temp ps
              ps = Ph.altitudePressure (gridSurfaceHeight ! ix) temp
              dewPointT = Ph.dewpointTemperature temp rh
              water = gridWaterVapor ! ix
              temp = gridAirTemp ! ix

-- | perform one iteration step given grid, start_time, time_step
--   returns updated grid
integrate :: (Monad m) => G.Grid -> UTCTime -> DiffTime -> m (G.Grid, UTCTime)
integrate grid t dt = do
  let sh = G.gridShape grid
      dts = diffTimeToSec  dt -- seconds
      t' = addUTCTime ( (fromRational . toRational) dt )  t -- new time
      gridAirTemp = G.gridAirTemp grid
      gridUWind = G.gridUWind grid
      gridVWind = G.gridVWind grid
      gridSurfaceHeight = G.gridSurfaceHeight grid
      gridVaporWater = G.gridVaporWater grid
      gridCondensedWater = G.gridCondensedWater grid
      gridEvapCoeff = G.gridEvapCoeff grid
      gridWaterDropletsSz = G.gridWaterDropletsSz grid
      gridPrecipationTriggered = G.gridPrecipationTriggered grid
      satAddMin min' a  b = max min' ( a + b ) -- ^ addition with saturation ny min' value
  -- temperature increase due to solar radiation influx
  let tempDiffS = R.fromFunction sh $
                  tempDiffSolarInflux grid gridCondensedWater t dts
  -- changes of temperature due to advection caused by wind
  let tempDiffA = R.fromFunction sh $
                  tempDiffAdvection grid gridAirTemp gridUWind gridVWind dts
  -- total change of temperature due to solar radiation and advection
  tempAbsorbed <- R.computeUnboxedP $ tempDiffS +^ tempDiffA
  -- we need to emit same amount of temperature as absorbed to keep planet temperature in balance
  -- each cell will emit T proportional to it's T^4, so cell with lower temp will emit much less
  -- than cells with higher T
  tempEmitted <- R.computeUnboxedP $ tempDiffEmission gridAirTemp tempAbsorbed
  tempTotal <- R.computeUnboxedP $ gridAirTemp +^ tempEmitted +^ tempAbsorbed
  -- 500Mb geopotential height
  h500Mb <- R.computeUnboxedP $ R.fromFunction sh (calculate500MbHeights gridAirTemp )
  -- V-winds
  vWind <- R.computeUnboxedP $ R.fromFunction sh (calculateVWind grid h500Mb)
  -- U-winds
  uWind <- R.computeUnboxedP $ R.fromFunction sh (calculateUWind grid h500Mb)
  -- surface water evaporation
  let waterVaporDiff = R.fromFunction sh
                  $ calculateWaterVaporDiff gridAirTemp gridEvapCoeff dts

  -- water condensation (cloud formation) and precipation (water loss due to rain)
  let waterCondensation = R.fromFunction sh
                         $ calculateWaterCondDiff gridVaporWater gridAirTemp gridSurfaceHeight dts
      waterPrecipation = R.map ( (*) dts . Ph.precipationRate )
                               $ R.zipWith (\a b -> if a then b else 0)
                                 gridPrecipationTriggered gridWaterDropletsSz

  waterCondensationDiff <- R.computeUnboxedP $ waterCondensation -^ waterPrecipation
  -- ! water vapor advection
  -- new water vapor density: + vapor - condensation
  waterVapor <- R.computeUnboxedP
                $ R.zipWith (satAddMin 0) gridVaporWater $ waterVaporDiff -^ waterCondensation
  -- cloud formation from condensed water
  waterCondensed <- R.computeUnboxedP $ R.zipWith (satAddMin 0) gridCondensedWater waterCondensationDiff
  -- condensed water must eventually rain
  let dropletSizeIncrease = R.map ( * (dts * Ph.rainDropletsFormRate) ) gridCondensedWater
      dropletScaleCoeff 0 _ = 0
      dropletScaleCoeff a b =  ( a / (a + b) )
      dropletsScaled = R.zipWith (*) gridWaterDropletsSz
                       $ R.zipWith dropletScaleCoeff gridCondensedWater waterCondensationDiff
  waterDropletsSz <- R.computeUnboxedP $ R.zipWith (satAddMin 0) dropletsScaled dropletSizeIncrease
  -- precipation
  let precipationProb = R.zipWith precipationChance waterDropletsSz gridPrecipationTriggered
      precipationChance dsz prec = round $ (if prec then 200 else 100) * Ph.precipationProbability dsz
      randSeed = (round . toRational . utcTimeToPOSIXSeconds) t
      randIntAr = randomishIntArray sh 10 100 randSeed
  precipationTriggered <- R.computeUnboxedP $ R.zipWith (<) randIntAr precipationProb

  return $! (grid {
               G.gridAirTemp = tempTotal
             , G.gridVWind = vWind
             , G.gridUWind = uWind
             , G.gridVaporWater = waterVapor
             , G.gridCondensedWater = waterCondensed
             , G.gridWaterDropletsSz = waterDropletsSz
             , G.gridPrecipationTriggered = precipationTriggered
             } , t')
