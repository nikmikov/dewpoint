module DewPoint.Physics where

import Data.Time.Clock
import Data.Time.Calendar(toGregorian, fromGregorian, diffDays)
import DewPoint.Geo

import Debug.Trace

-- | aboslute zero in K
tZero :: Floating a => a
tZero = -273.15

-- | convert Kelvin to Celsius
toCelsius :: Floating a => a -> a
toCelsius = (+) tZero

-- | Celsius to Kelvin
toKelvin :: Floating a => a -> a
toKelvin v = v - tZero

-- | Average surface temperature
avgSurfaceTemp :: Floating a => a
avgSurfaceTemp = toKelvin 14

-- | degrees to radians
toRad :: Floating a => a -> a
toRad a = pi * a / 180.0

toDegrees :: Floating a => a -> a
toDegrees r = 180.0 * r / pi

solarConst :: Floating a => a
solarConst = 1367 -- W/m^2

secondsInADay :: Floating a => a
secondsInADay = 24*60*60

-- | sea level pressure
pressureSeaLevel :: Floating a => a
pressureSeaLevel = 101325

-- | gravitational const
gravConst :: Floating a => a
gravConst = 9.807

-- | molar mass of dry air
molarMassAir :: Floating a => a
molarMassAir = 0.0289644

-- | molar mass of dry air
molarMassWater :: Floating a => a
molarMassWater = 0.01801528

-- | universal gas constant
gasConst :: Floating a => a
gasConst = 8.3143

-- | given pressure in atmospheres calculate geopotential height
calculatePressureHeight :: Floating a => a -> a -> a
calculatePressureHeight ps temp = -(log ps) * gasConst * temp / (molarMassAir * gravConst)

-- | return pressure at given altitude in Pa
altitudePressure :: (Floating a, Eq a) => a -> a -> a
altitudePressure _ 0 = 0
altitudePressure alt temp = pressureSeaLevel *
                            exp ( -gravConst * molarMassAir * alt / (gasConst * temp) )

-- | pressure of water component in atmosphere
partialAtmosphericWaterPressure :: Floating a => a -- ^ water denisty
                                -> a -- ^ air temperature
                                -> a
partialAtmosphericWaterPressure wd temp = wd * gasConst * temp / molarMassWater

-- | Using approximation Buck formula
equilibriumWaterPressure  :: Floating a => a -- ^ surface pressure
                          -> a -- ^ air temperature
                          -> a
equilibriumWaterPressure ps temp = let pb = ps / 100 -- convert to millibars
                                       tC = toCelsius temp
                                       tcoeff = 17.502 * tC / (240.97 + tC)
                                       e = exp 1
                                   in 1.0007 + 3.46e-6 * pb * 6.1121 * e**tcoeff

-- | relative humdity 0..1
relativeHumidity :: Floating a => a -- ^ water density
                 -> a -- ^ temperature
                 -> a -- ^ surface pressure
                 -> a
relativeHumidity wd temp ps = partialAtmosphericWaterPressure wd temp / equilibriumWaterPressure ps temp


-- | Calculate coriolis force at the given latitude
coriolisForce :: Floating a => Lat -> a
coriolisForce lat =  -2 * 7.292e-5 * sin(toRad $ fromRational $ toRational $ fromLat lat)

-- | Return coriolis force at given latitude
--   but not less than 4e-5 in northern hemisphere and -4e-5 southern hemisphere
coriolisForce' :: (Floating a, Ord a) => Lat -> a
coriolisForce' = clamp' 4e-5 . coriolisForce
    where clamp' m' 0 = m'
          clamp' m' v = (max m' (abs v) )  * (signum v)

-- | Calculates geostrophic zonal wind (west-east) ( U-wind )
--   given latitude,
--         gradient of pressure at 500Mb heights by Y axis(Pa/meters)
--         gradient step (meters)
--   return U-wind in m/sec
geostrophicZonalWind :: (Floating a, Ord a) => Lat -> a -> a
geostrophicZonalWind lat gradY = gradY * gravConst / (coriolisForce' lat)

-- | Calculates geostrophic meridional wind (north-south) ( V-wind )
--   given latitude,
--         gradient of pressure at 500Mb heights by X axis(Pa/meters)
--   return V-wind in m/sec
geostrophicMeridionalWind :: (Floating a, Ord a) => Lat -> a -> a
geostrophicMeridionalWind lat gradX = gradX * gravConst / (coriolisForce' lat)


-- | extract seconds since midnight from UTCTime and convert it to Double
secondsSinceMidnight :: Floating a => UTCTime -> a
secondsSinceMidnight = fromRational . toRational . utctDayTime

-- | return declination of sun angle for given day of year in radians
declinationOfSunAngle :: Floating a => Integer -> a
declinationOfSunAngle n = let n' = fromIntegral n
                              val = 0.98565 * (n' + 10.0) + 1.914 * sin ( toRad(0.98565 * (n' - 2)) )
                          in -asin( 0.39779 * cos (toRad val) )

-- | calculate solar hour angle in radians - angle from solar noon
solarHourAngle :: Floating a  => Lon -> UTCTime -> a
solarHourAngle lon t =
    let secondsToDegrees = 180.0 - (secondsSinceMidnight t) * 360.0 / secondsInADay
    in toRad (fromLon lon - secondsToDegrees)

-- | calculate solar zenith angle in radians
solarZenithAngle :: Floating a => Lat -> Integer -> a -> a
solarZenithAngle lat dayOfYear hourAngle =
  let da = declinationOfSunAngle dayOfYear
      latr = toRad ( fromLat lat)
      a = sin(latr) * sin (da) + cos (latr) * cos (da) * cos (hourAngle)
  in acos(a)

-- | return solar energy influx (W/m^2) at the given coordinate and time in UTC
solarEnergyInflux :: (Ord a, Floating a) => (Lat, Lon) -> UTCTime -> a
solarEnergyInflux (lat, lon) t = let (year, _, _) = toGregorian $ utctDay t
                                     janFst =  fromGregorian year 1 1
                                     dayOfYear = diffDays (utctDay t) janFst
                                     ha = solarHourAngle lon t
                                     sza = solarZenithAngle lat dayOfYear ha
                                 in if cos (sza) > 0 then solarConst * cos (sza) else 0

-- | calculate temperature increase (K/s) at given coordinate
--   cloudCover: cloud cover coefficient (0..1)
temperatureDiffPerSecond :: (Ord a, Floating a) => (Lat, Lon) -> UTCTime -> a -> a
temperatureDiffPerSecond coord t cloudCover =
    let rawSolarInflux = solarEnergyInflux coord t
        albedo = 0.3 + (0.6 * cloudCover) -- heuristic formula
        absorbedSolarInflux = (1.0 - albedo) * rawSolarInflux
    in absorbedSolarInflux / 1.8e6 -- experimental formula

-- |  Very primitive calculation of evaporation of water (kilogramms/m^3 per second)
-- with temp <= 0 Celsium no evaporation
-- with temp > 0 increase relative humidity is linear relative to air temp and max with t=50C
-- second parameter is evaporation coefficent [0..1]
-- return evaporation of water from surface (kg/ (m^3 * second) )
surfaceEvaporation :: (Ord a, Floating a) => a -> a -> a
surfaceEvaporation ta ec = let t = toCelsius ta
                               tc = (min t 50) / 50
                               pconst = 1e-9 -- kg/ (m^3 * sec)
                           in if t <= 0 then 0
                              else ec * tc * pconst
