module DewPoint.Simulation (integrate) where

import Data.Time.Clock
import Data.Array.Repa ( (!) )
import qualified Data.Array.Repa as R
import qualified DewPoint.Grid as G
import DewPoint.Physics

--type DUArray = R.Array R.U R.DIM2 Double

-- | temperature increase due to solar influx at given cell at given time interval in seconds
tempDiffSolarInflux :: G.Grid -> UTCTime -> Int -> G.Ix -> Double
tempDiffSolarInflux g t dt ix = temperatureDiffPerSecond ( (G.gridCellGeoCoord g) ! ix ) t (fromIntegral dt)


-- | perform one iteration step given grid, start_time, time_step
--   returns updated grid
integrate :: (Monad m) => G.Grid -> UTCTime -> DiffTime -> m (G.Grid, UTCTime)
integrate grid t dt = do let sh = G.gridShape grid
                             (dts, _) = properFraction  dt -- seconds
                             tempDiffFunc = tempDiffSolarInflux grid t dts
                             t' = addUTCTime ( (fromRational . toRational) dt )  t -- new time
                         -- temperature increase due to solar radiation influx
                         tempDiff <- R.computeUnboxedP $ R.fromFunction sh tempDiffFunc
                         return $! (grid {
                           G.gridAirTemp = R.computeS $ R.zipWith (+) (G.gridAirTemp grid) tempDiff
                         } , t')
