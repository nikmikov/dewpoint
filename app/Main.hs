module Main where

import Paths_dewpoint (version)
import Data.Version (showVersion)
import Data.Time.Clock(UTCTime(..), DiffTime, secondsToDiffTime, diffUTCTime)
import Data.Time.Format(formatTime, defaultTimeLocale)
import qualified CommandLineOptions as Opt
import Options.Applicative((<>), execParser, info, helper, fullDesc, progDesc, header)

import Data.List(intersperse)
import Data.List.Split(splitOn)
import qualified Data.Set as S

import qualified DewPoint.Grid as G
import qualified DewPoint.Geo as Geo
import qualified DewPoint.Physics as Ph
import qualified DewPoint.Simulation as Sim
import qualified Data.Array.Repa as R
import System.FilePath
import Text.Printf
import Control.Monad(when)

-- | Return version string
versionInfo :: String
versionInfo = "Weather simulation v." ++ showVersion version

data MeteoStation = MeteoStation
                    {
                      iataCode :: String
                    , country  :: String
                    , city     :: String
                    , stationName :: String
                    , latLon   :: (Double, Double)
                    , gridCell :: G.Ix
                    }

-- | output station data in format:
--   SYD|-33.86,151.21,39|2015-12-23T05:02:12Z|Rain|+12.5|1004.3|97
stationString :: G.Grid -> UTCTime -> MeteoStation -> String
stationString g t st = printf "%3s|%9.5f,%9.5f|%s|%s|%+5.1f|%f|%d|%10.4f"
                       (iataCode st)
                       lat lon
                       (formatTime defaultTimeLocale "%FT%TZ" t)
                       "Rain"
                       (Ph.toCelsius temp)
                       ps
                       rh
                       se
  where ix = gridCell st
        (lat, lon) = latLon st
        temp = (flip R.index ix . G.gridAirTemp ) g
        alt = (flip R.index ix . G.gridSurfaceHeight ) g
        ps  = Ph.altitudePressure (fromIntegral alt) temp
        wd  = (flip R.index ix . G.gridWaterDensity ) g
        rh  = toInteger $ round $ Ph.relativeHumidity wd temp ps
        se = Ph.solarEnergyInflux (Geo.toLatLon lat lon) t


loadStations :: String -> S.Set String -> G.Grid -> IO [MeteoStation]
loadStations dataDir st g = do
  let f = map (toStation) . filter (flip S.member st . head ) . map (splitOn ";") . lines
  content <- readFile $ dataDir </> "stations.csv"
  let stations = f content
  if (length stations == length st)
  then return stations
  else fail $ "Stations not found: " ++ toDiffStr stations
    where toStation (iata':country':city':name':lat':lon':[])
              = MeteoStation iata' country' city' name' (lat, lon)
                $ G.gridGetLocationIx g (Geo.toLatLon lat lon)
                    where (lat,lon) = (read lat', read lon')
          toStation _ = error "Malformed station line in input CSV"
          toDiffStr = concat . intersperse ", " . S.toList . S.difference st . S.fromList . map iataCode

validateOpt :: Monad m => Opt.ProgramOpt -> m Opt.ProgramOpt
validateOpt o = if (Opt.outputType o) == Opt.Text && null (Opt.stationListIATACodes o)
                then fail "Must specify at least one station when --output set to \"text\""
                else return o


printOutput :: [MeteoStation]
            -> UTCTime           -- ^ start time
            -> Integer           -- ^ print interval in seconds
            -> (G.Grid, UTCTime) -- ^ (grid, currentTime)
            -> IO (G.Grid, UTCTime)
printOutput st t0 dt (g, t) =
    do let secSinceStart = round $ realToFrac $ diffUTCTime t t0
       when (secSinceStart `mod` dt == 0) $
            mapM_ (putStrLn . stationString g t) st
       return (g, t)

initGrid :: G.Grid
initGrid = let g = G.gridCreateFixed 100 200
               sh = G.gridShape g
           in g {
                    G.gridAirTemp = (R.computeS . R.fromFunction sh . const) Ph.avgSurfaceTemp
                }

run :: Opt.ProgramOpt -> IO()
run opt = do
  let grid   = initGrid
      tstart = UTCTime (Opt.startDate opt) 0
      timeInterval = secondsToDiffTime (Opt.interval  opt) * 60
      printInterval = (Opt.outputInterval  opt) * 60
      st = S.fromList (Opt.stationListIATACodes opt)
  stations <- loadStations (Opt.dataDir opt) st grid
  loop stations tstart timeInterval printInterval (grid, tstart)
    where loop st t0 dt dt' (g, t) = Sim.integrate g t dt
                                     >>= printOutput st t0 dt'
                                     >>= loop st t0 dt dt'

main :: IO ()
main = execParser opts >>= validateOpt >>= run
  where
    opts = info (helper <*> Opt.parseOptions)
      ( fullDesc
     <> progDesc "Run weather simulation"
     <> header versionInfo )
