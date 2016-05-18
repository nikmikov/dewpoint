module Main where

import Paths_dewpoint (version)
import Data.Version (showVersion)
import Data.Time.Clock(UTCTime(..), DiffTime,  secondsToDiffTime)
import qualified CommandLineOptions as Opt
import Options.Applicative((<>), execParser, info, helper, fullDesc, progDesc, header)

import Data.List(intersperse)

import qualified DewPoint.Grid as G
import qualified DewPoint.Physics as Ph
import qualified DewPoint.Simulation as Sim
import qualified Data.Array.Repa as R

-- | Return version string
versionInfo :: String
versionInfo = "Weather simulation v." ++ showVersion version

data MeteoStation = MeteoStation
                    {
                      iataCode :: String
                    , latLon   :: (Double, Double)
                    , gridCell :: G.Ix
                    }

-- | output station data in format:
--   SYD|-33.86,151.21,39|2015-12-23T05:02:12Z|Rain|+12.5|1004.3|97
stationString :: G.Grid -> UTCTime -> MeteoStation -> String
stationString g t st = concat $ intersperse "|"
  [
    iataCode st
  , (show lat ++ "," ++ show lon)
  , show t
  , "Rain"
  , show $ Ph.toCelsius temp
  , show ps
  , show rh
  ]
  where ix = gridCell st
        (lat, lon) = latLon st
        temp = (flip R.index ix . G.gridAirTemp ) g
        alt = (flip R.index ix . G.gridSurfaceHeight ) g
        ps  = Ph.altitudePressure (fromIntegral alt) temp
        wd  = (flip R.index ix . G.gridWaterDensity ) g
        rh  = Ph.relativeHumidity wd temp ps


validateOpt :: Monad m => Opt.ProgramOpt -> m Opt.ProgramOpt
validateOpt o = if (Opt.outputType o) == Opt.Text && null (Opt.stationListIATACodes o)
                then fail "Must specify at least one station when --output set to \"text\""
                else return o

integrate' :: Opt.ProgramOpt -> G.Grid -> UTCTime -> DiffTime -> IO (G.Grid, UTCTime)
integrate' _ g t dt = do (g', t') <- Sim.integrate g t dt
                         let st = MeteoStation "SYD" (-33.86,151.21) (G.xyToIx 20 20)
                         print $ stationString g t st
                         return $! (g', t')

run :: Opt.ProgramOpt -> IO()
run opt = do
  let grid   = G.gridCreateFixed 100 200
      tstart = UTCTime (Opt.startDate opt) 0
      timeInterval = secondsToDiffTime (Opt.interval  opt) * 60
  loop timeInterval (grid, tstart)
    where loop dt (g, t) = integrate' opt g t dt >>=  loop dt

main :: IO ()
main = execParser opts >>= validateOpt >>= run
  where
    opts = info (helper <*> Opt.parseOptions)
      ( fullDesc
     <> progDesc "Run weather simulation"
     <> header versionInfo )
