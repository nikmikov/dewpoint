module DewPoint.Physics where

import Data.Time.Clock
import Data.Time.Calendar(toGregorian, fromGregorian, diffDays)
import DewPoint.Geo

import Debug.Trace

-- | aboslute zero in K
tZero :: Double
tZero = -273.15

-- | convert Kelvin to Celsius
toCelsius :: Double -> Double
toCelsius = (+) tZero

-- | Celsius to Kelvin
toKelvin :: Double -> Double
toKelvin v = v - tZero

-- | Average surface temperature
avgSurfaceTemp :: Double
avgSurfaceTemp = toKelvin 14

-- | degrees to radians
toRad :: Double -> Double
toRad a = pi * a / 180.0

toDegrees :: Double -> Double
toDegrees r = 180.0 * r / pi

solarConst :: Double
solarConst = 1367 -- W/m^2

secondsInADay :: Double
secondsInADay = 24*60*60

-- | sea level pressure
pressureSeaLevel :: Double
pressureSeaLevel = 101325

-- | gravitational const
gravConst :: Double
gravConst = 9.807

-- | molar mass of dry air
molarMassAir :: Double
molarMassAir = 0.0289644

-- | molar mass of dry air
molarMassWater :: Double
molarMassWater = 0.01801528

-- | universal gas constant
gasConst :: Double
gasConst = 8.3143

-- | return pressure at given altitude in Pa
altitudePressure :: Double -> Double -> Double
altitudePressure _ 0 = 0
altitudePressure alt temp = pressureSeaLevel * exp ( -gravConst * molarMassAir * alt / (gasConst * temp) )

-- | pressure of water component in atmosphere
partialAtmosphericWaterPressure :: Double -- ^ water denisty
                                -> Double -- ^ air temperature
                                -> Double
partialAtmosphericWaterPressure wd temp = wd * gasConst * temp / molarMassWater

-- | Using approximation Buck formula
equilibriumWaterPressure  :: Double -- ^ surface pressure
                          -> Double -- ^ air temperature
                          -> Double
equilibriumWaterPressure ps temp = let pb = ps / 100 -- convert to millibars
                                       tC = toCelsius temp
                                       tcoeff = 17.502 * tC / (240.97 + tC)
                                       e = exp 1
                                   in 1.0007 + 3.46e-6 * pb * 6.1121 * e**tcoeff

-- | relative humdity 0..1
relativeHumidity :: Double -- ^ water density
                 -> Double -- ^ temperature
                 -> Double -- ^ surface pressure
                 -> Double
relativeHumidity wd temp ps = partialAtmosphericWaterPressure wd temp / equilibriumWaterPressure ps temp


-- | extract seconds since midnight from UTCTime and convert it to Double
secondsSinceMidnight :: UTCTime -> Double
secondsSinceMidnight = fromRational . toRational . utctDayTime

-- | return declination of sun angle for given day of year in radians
declinationOfSunAngle :: Integer -> Double
declinationOfSunAngle n = let n' = fromIntegral n
                              val = 0.98565 * (n' + 10.0) + 1.914 * sin ( toRad(0.98565 * (n' - 2)) )
                          in -asin( 0.39779 * cos (toRad val) )

-- | calculate solar hour angle in radians - angle from solar noon
solarHourAngle :: Lon -> UTCTime -> Double
solarHourAngle (Lon lon) t =
    let secondsToDegrees = 180.0 - (secondsSinceMidnight t) * 360.0 / secondsInADay
    in toRad (lon - secondsToDegrees)

-- | calculate solar zenith angle in radians
solarZenithAngle :: Lat -> Integer -> Double -> Double
solarZenithAngle (Lat lat) dayOfYear hourAngle =
  let da = declinationOfSunAngle dayOfYear
      latr = toRad lat
      a = sin(latr) * sin (da) + cos (latr) * cos (da) * cos (hourAngle)
  in acos(a)

-- | return solar energy influx (W/m^2) at the given coordinate and time in UTC
solarEnergyInflux :: (Lat, Lon) -> UTCTime -> Double
solarEnergyInflux (lat, lon) t = let (year, _, _) = toGregorian $ utctDay t
                                     janFst =  fromGregorian year 1 1
                                     dayOfYear = diffDays (utctDay t) janFst
                                     ha = solarHourAngle lon t
                                     sza = solarZenithAngle lat dayOfYear ha
                                 in if cos (sza) > 0 then solarConst * cos (sza) else 0

-- | calculate temperature increase (K/s) at given coordinate
--   cloudCover: cloud cover coefficient (0..1)
temperatureDiffPerSecond :: (Lat, Lon) -> UTCTime -> Double -> Double
temperatureDiffPerSecond coord t cloudCover =
    let rawSolarInflux = solarEnergyInflux coord t
        albedo = 0.3 + (0.6 * cloudCover) -- heuristic formula
        absorbedSolarInflux = (1.0 - albedo) * rawSolarInflux
    in absorbedSolarInflux / 1.2e6 -- experimental formula
