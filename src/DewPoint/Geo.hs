module DewPoint.Geo(Lat(..), Lon(..), toLat, toLon, fromLat, fromLon, toLatLon,
                    earthRadius, earthSurfaceArea) where

data Lat = Lat Double deriving (Ord, Eq, Show)
data Lon = Lon Double deriving (Ord, Eq, Show)

instance Bounded Lat where
  minBound = Lat (-90)
  maxBound = Lat 90

instance Bounded Lon where
  minBound = Lon 0
  maxBound = Lon 360

clamp :: (Ord a, Bounded a) => a -> a
clamp val = if val > minBound
            then if val < maxBound then val else maxBound
            else minBound

toLat :: Double -> Lat
toLat val = clamp (Lat val)

fromLat :: Floating a => Lat -> a
fromLat (Lat v) = (fromRational . toRational ) v

fromLon :: Floating a => Lon -> a
fromLon (Lon v) = (fromRational . toRational ) v

toLon :: Double -> Lon
toLon val = if val < 0 then toLon (val + 360.0)  else clamp (Lon val)

toLatLon :: Double -> Double -> (Lat, Lon)
toLatLon lat lon = (toLat lat, toLon lon)

-- | radius of earth (m)
earthRadius :: Floating a => a
earthRadius = 6378100.0

earthSurfaceArea :: Floating a => a
earthSurfaceArea = 4 * pi * earthRadius**2
