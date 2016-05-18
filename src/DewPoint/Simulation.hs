module DewPoint.Simulation (integrate, paramAdvectionDiff) where

import Data.Time.Clock
import Data.Array.Repa ( (!), (+^) )
import qualified Data.Array.Repa as R
import qualified DewPoint.Grid as G
import DewPoint.Physics

-- | temperature increase due to solar influx at given cell at given time interval in seconds
tempDiffSolarInflux :: G.Grid -> UTCTime -> Double -> G.Ix -> Double
tempDiffSolarInflux g t dt ix = temperatureDiffPerSecond ( (G.gridCellGeoCoord g) ! ix ) t dt

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
      totalArea = l * l
      c0 = (l - dy) * (l - dx) / totalArea
      cX = (l - dy) * dx / totalArea
      cY = (l - dx) * dy / totalArea
      cXY = 1 - c0 - cX - cY
      cX' = dx / (dx + dy)
      cY' = 1 - cX'
      adX = aXs * ( cX + cX' * cXY ) -- ^ loss quantity at source cell X
      adY = aYs * ( cY + cY' * cXY ) -- ^ loss quantity at source cell Y
      a1 = c0 * a0 + adX + adY      -- ^ new temperature in the cell
  in (a1 - a0, adX, adY)

{-
tempDiffAdvection :: G.Grid -> Double -> G.Ix -> Double
tempDiffAdvection g dt ix = let u = (G.gridUWind g) ! ix
                                v = (G.gridUWind g) ! ix
                                t0 = (G.gridAirTemp g) ! ix
                                ixU = G.gridSrcCellX g ix u
                                ixV = G.gridSrcCellY g ix v
                                tXs = (G.gridAirTemp g) ! ixU
                                tYs = (G.gridAirTemp g) ! ixV
                                l = G.gridCellLengthMeters g
                                (dT, _, _) = paramAdvectionDiff u v tXs tYs l dt t0
                            in dT
-}

-- | perform one iteration step given grid, start_time, time_step
--   returns updated grid
integrate :: (Monad m) => G.Grid -> UTCTime -> DiffTime -> m (G.Grid, UTCTime)
integrate grid t dt = do let sh = G.gridShape grid
                             dts = (fromRational . toRational)  dt -- seconds
                             tempDiffFunc = tempDiffSolarInflux grid t dts
                             t' = addUTCTime ( (fromRational . toRational) dt )  t -- new time
                         -- temperature increase due to solar radiation influx
                         let tempDiffS = R.fromFunction sh tempDiffFunc
                         -- changes of temperature due to advection caused by wind
--                         let tempDiffA = R.fromFunction sh (tempDiffAdvection grid dts)
                         -- total changes
                         let tempNew =  (G.gridAirTemp grid) +^ tempDiffS -- +^ tempDiffS
--                         tempNew'   <- R.computeUnboxedS  tempNew
                         return $! (grid {
                           G.gridAirTemp = R.computeUnboxedS  tempNew
                         } , t')
