module MeteoStation where


import Data.Time.Clock(UTCTime(..))
import Data.Time.Format(formatTime, defaultTimeLocale)
import qualified Data.Set as S

import qualified DewPoint.Grid as G
import qualified DewPoint.Geo as Geo
import qualified DewPoint.Physics as Ph
import qualified Data.Array.Repa as R
import System.FilePath
import Text.Printf
import qualified Data.Text as T
import Data.Text(snoc)
import qualified Data.Text.IO as T

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
stationString :: Bool -> G.Grid -> UTCTime -> MeteoStation -> String
stationString verbose g t st =
    if verbose
    then
        printf
        "%3s|%9.5f,%9.5f|%s|%13s|%+5.1f|%5.1f|%3d|%10.4f|%10.4f|%10.4f|%10.10f|%12.10f|%2.2f|%12.10f"
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
    else
        printf
        "%3s|%9.5f,%9.5f|%s|%13s|%+5.1f|%5.1f|%3d"
        (iataCode st)
        lat lon
        (formatTime defaultTimeLocale "%FT%TZ" t)
        rainTrig
        (Ph.toCelsius temp)
        ps
        rh
  where ix = gridCell st
        (lat, lon) = latLon st
        temp = (flip R.index ix . G.gridAirTemp ) g
        alt = (flip R.index ix . G.gridSurfaceHeight ) g
        ps  = Ph.altitudePressure alt temp / 100
        wvd  = (flip R.index ix . G.gridVaporWater ) g
        wcd  = (flip R.index ix . G.gridCondensedWater ) g
        rh  = (round::Float->Int) $ 100 * Ph.relativeHumidity wvd temp ps
        se :: Double
        se  = Ph.solarEnergyInflux (Geo.toLatLon lat lon) t
        v   = (G.gridVWind g) R.! ix
        u   = (G.gridUWind g) R.! ix
        cloud = Ph.cloudCover wcd
        rainDroplets = (G.gridWaterDropletsSz g) R.! ix
        rainTrig = stateString  ((G.gridPrecipationTriggered g) R.! ix) temp cloud


stateString :: (Floating a, Ord a) => Bool -> a -> a -> String
stateString precTriggered temp cloudCover
    | precTriggered && temp < (-5)         = "Snow"
    | precTriggered                        = "Rain"
    | cloudCover > 0.2 && cloudCover < 0.6 = "Partly Cloudy"
    | cloudCover >= 0.6                    = "Cloudy"
    | otherwise                            = "Sunny"

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
