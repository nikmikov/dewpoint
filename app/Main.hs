--{-# LANGUAGE OverloadedStrings #-}
module Main where

import Paths_dewpoint (version)
import Data.Version (showVersion)
import Data.Time.Clock(UTCTime(..), DiffTime, secondsToDiffTime, diffUTCTime)
import Data.Time.Format(formatTime, defaultTimeLocale)
import qualified CommandLineOptions as Opt
import qualified LoadGrid as LG
import qualified Plot as Plot
import Options.Applicative((<>), execParser, info, helper, fullDesc, progDesc, header)

import qualified Data.Set as S

import qualified DewPoint.Grid as G
import qualified DewPoint.Geo as Geo
import qualified DewPoint.Physics as Ph
import qualified DewPoint.Simulation as Sim
import qualified Data.Array.Repa as R
import System.FilePath
import Text.Printf
import Control.Monad(when)
import qualified Data.ByteString as BS
import qualified Data.ByteString.Char8 as BS8
import qualified Data.Text as T
import qualified Data.Text.Encoding as T
import Data.Text(snoc)
import qualified Data.Text.IO as T
import Data.Foldable(toList)

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
                    } deriving (Show)

-- | output station data in format:
--   SYD|-33.86,151.21,39|2015-12-23T05:02:12Z|Rain|+12.5|1004.3|97
stationString :: G.Grid -> UTCTime -> MeteoStation -> String
stationString g t st =
    printf
    "%3s|%9.5f,%9.5f|%s|%5s|%+5.1f|%f|%3d|%10.4f|%10.4f|%10.4f|%10.10f|%12.10f|%2.2f|%12.10f"
    (iataCode st)
    lat lon
    (formatTime defaultTimeLocale "%FT%TZ" t)
    rainTrig
    (Ph.toCelsius temp)
    ps
    rh
    se
    u
    v
    wvd
    wcd
    cloud
    rainDroplets
  where ix = gridCell st
        (lat, lon) = latLon st
        temp = (flip R.index ix . G.gridAirTemp ) g
        alt = (flip R.index ix . G.gridSurfaceHeight ) g
        ps  = Ph.altitudePressure alt temp
        wvd  = (flip R.index ix . G.gridVaporWater ) g
        wcd  = (flip R.index ix . G.gridCondensedWater ) g
        rh  = toInteger $ round $ 100 * Ph.relativeHumidity wvd temp ps
        se :: Double
        se  = Ph.solarEnergyInflux (Geo.toLatLon lat lon) t
        v   = (G.gridVWind g) R.! ix
        u   = (G.gridUWind g) R.! ix
        cloud = Ph.cloudCover wcd
        rainDroplets = (G.gridWaterDropletsSz g) R.! ix
        rainTrig = show $ (G.gridPrecipationTriggered g) R.! ix


loadStations :: String -> S.Set T.Text -> G.Grid -> IO [MeteoStation]
loadStations dataDir st g = do
  let f = map (toStation . map T.unpack) . filter (flip S.member st . head ) . map split' . T.lines
      split' = T.splitOn (T.pack ";")
  content <- T.readFile $ dataDir </> "stations.csv"
  let stations = f content
  if (length stations == length st)
  then return stations
  else fail $ "Stations not found: " ++ toDiffStr stations
    where toStation (iata':country':city':name':lat':lon':[])
              = MeteoStation iata' country' city' name' (lat, lon)
                $ G.gridGetLocationIx g (Geo.toLatLon lat lon)
                    where (lat,lon) = (read lat', read lon')
          toStation _ = error "Malformed station line in input CSV"
          setToStr = T.unpack . T.concat . map( `snoc` ',' ) .  S.toList
          toDiffStr = setToStr . S.difference st . S.fromList . map (T.pack . iataCode)

validateOpt :: Monad m => Opt.ProgramOpt -> m Opt.ProgramOpt
validateOpt o = if (Opt.outputType o) == Opt.Text && null (Opt.stationListIATACodes o)
                then fail "Must specify at least one station when --output set to \"text\""
                else return o


printOutput :: ( (G.Grid, UTCTime) -> BS.ByteString) -- ^ print function
            -> UTCTime           -- ^ start time
            -> Integer           -- ^ print interval in seconds
            -> (G.Grid, UTCTime) -- ^ (grid, currentTime)
            -> IO (G.Grid, UTCTime)
printOutput pFn t0 dt (g, t) =
    do let secSinceStart = (round . toRational) $ diffUTCTime t t0
       when (secSinceStart `mod` dt == 0) $
            BS8.putStrLn $ pFn (g, t)
       return (g, t)

printStationdData :: [MeteoStation] -> (G.Grid, UTCTime) -> BS.ByteString
printStationdData = flip f'
    where f' (g, t) = T.encodeUtf8 . T.intercalate (T.pack "\n") . fmap (T.pack . stationString g t)

plotHeader :: G.Grid -> BS.ByteString
plotHeader g = Plot.encodeHeader "Temperature and U,V winds plot"
               (toList $ Geo.fromLat <$> G.gridLats g)
               (toList $ Geo.fromLon <$> G.gridLons g)

plotFrame :: (G.Grid, UTCTime) -> BS.ByteString
plotFrame (g, t) = Plot.encodeFrame t
                   ( (map Ph.toCelsius .  R.toList . R.transpose . G.gridAirTemp)  g)
                   (R.toList $ R.transpose $ G.gridUWind g)
                   (R.toList $ R.transpose $ G.gridVWind g)

initGrid :: G.Grid
initGrid = g {
             G.gridAirTemp = (R.computeS . R.fromFunction sh . const) Ph.avgSurfaceTemp
           , G.gridEvapCoeff = (R.computeS . R.fromFunction sh . const) 1
           }
    where g = G.gridCreateFixed 100 200
          sh = G.gridShape g

run :: Opt.ProgramOpt -> IO()
run opt = do
  let
      tstart = UTCTime (Opt.startDate opt) 0
      timeInterval = secondsToDiffTime (Opt.interval  opt) * 60
      printInterval = (Opt.outputInterval  opt) * 60
      st = S.fromList (Opt.stationListIATACodes opt)
      dataDir = Opt.dataDir opt
      gridFile = dataDir </> "grid.json"
  maybeGrid <- LG.loadGrid <$> BS.readFile gridFile
  grid <- case maybeGrid of
    (Just a) -> return a
    _        -> fail $ "Unable to load grid from " ++ show gridFile
                       ++ " run scripts/import_data.py to import data"

  stations <- loadStations  dataDir st grid
  let isPlot = Opt.Plot  == Opt.outputType opt
      printFn = if isPlot then plotFrame else printStationdData stations
  when isPlot $  (BS8.putStrLn . plotHeader) grid
  loop printFn tstart timeInterval printInterval (grid, tstart)
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
