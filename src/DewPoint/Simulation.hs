{-# LANGUAGE BangPatterns #-}
module DewPoint.Simulation (integrate, paramAdvectionDiff) where

import Data.Time.Clock
import Data.Array.Repa ( (!), (+^) )
import qualified Data.Array.Repa as R
import qualified DewPoint.Grid as G
import DewPoint.Physics
import Control.Exception.Base

-- | temperature increase due to solar influx at given cell at given time interval in seconds
tempDiffSolarInflux :: G.Grid
                    -> UTCTime
                    -> Double
                    -> G.Ix
                    -> Double
tempDiffSolarInflux g t dt ix = let tdps = temperatureDiffPerSecond
                                           ( (G.gridCellGeoCoord g) ! ix )
                                           t
                                           ( (G.gridCloudCover g) ! ix )
                                in tdps * dt

-- | parameter advection due to winds
--
--     -----> U
--
--       dx
--     -------
--     |_____|_______  _
--  |  |     |       |  |
--  |  | cXY |  cX   |  |  dy
--  V  |_____|_______| _|
--     |     |       |
--     | cY  |  c0   |
--     |     |       |
--     ---------------
--
--
paramAdvectionDiff :: Double -- ^ U-wind along X (m/s)
                   -> Double -- ^ V-wind along Y (m/s)
                   -> Double -- ^ x-source cell quantity
                   -> Double -- ^ y-source cell quantity
                   -> Double -- ^ length of cell side (assuming cell is always a square)
                   -> Double -- ^ time interval (sec)
                   -> Double -- ^ initial quantity
                   -> (Double, Double, Double) -- ^ change of quantity in the cell + influx of quantity  (adX, adY)
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


tempDiffAdvection :: G.Grid -> Double -> G.Ix -> Double
tempDiffAdvection !g !dt !ix = let u = (G.gridUWind g) ! ix
                                   v = (G.gridUWind g) ! ix
                                   t0 = (G.gridAirTemp g) ! ix
                                   ixU = G.gridSrcCellX g ix u
                                   ixV = G.gridSrcCellY g ix v
                                   tXs = (G.gridAirTemp g) ! ixU
                                   tYs = (G.gridAirTemp g) ! ixV
                                   l = G.gridCellLengthMeters g
                                   (dT, _, _) = paramAdvectionDiff u v tXs tYs l dt t0
                               in assert ( (not . isNaN) dT ) dT

diffTimeToSec :: DiffTime -> Double
diffTimeToSec = fromRational . toRational

-- | perform one iteration step given grid, start_time, time_step
--   returns updated grid
integrate :: (Monad m) => G.Grid -> UTCTime -> DiffTime -> m (G.Grid, UTCTime)
integrate grid t dt = do
  let sh = G.gridShape grid
      dts = diffTimeToSec  dt -- seconds
      tempDiffFunc = tempDiffSolarInflux grid t dts
      t' = addUTCTime ( (fromRational . toRational) dt )  t -- new time
  -- temperature increase due to solar radiation influx
  let tempDiffS = R.fromFunction sh tempDiffFunc
  -- changes of temperature due to advection caused by wind
  let tempDiffA = R.fromFunction sh (tempDiffAdvection grid dts)
  -- total increase of temperature
  tempAbsorbed <- R.computeUnboxedP $ tempDiffS +^ tempDiffA
  let tempAbsorbedSum = R.sumAllS tempAbsorbed
  -- we need to emit same amount of temperature as absorbed to keep planet temperature in balance
  -- each cell will emit T proportional to it's T^8, so cell with lower temp will emit much less
  -- than cells with higher T
  let minTemp = toKelvin (-50)
      minBoundT te' = if te' < minTemp then 0 else te' - minTemp
      adjustTemp te' = (minBoundT te') ** 4
  tempE <- R.computeUnboxedP $ R.map adjustTemp (G.gridAirTemp grid)
  let tempEsum = R.sumAllS tempE
  tempCellEmissionCoeff <- R.computeUnboxedP $ R.map (\v -> -v / tempEsum) tempE
  tempEmitted <- R.computeUnboxedP $ R.map ( * tempAbsorbedSum ) tempCellEmissionCoeff
  tempTotal <- R.computeUnboxedP $ (G.gridAirTemp grid) +^ tempEmitted +^ tempAbsorbed
  -- total changes
--  tempNew'   <- R.computeUnboxedS  tempNew
  return $! (grid {
               G.gridAirTemp = tempTotal
             } , t')
