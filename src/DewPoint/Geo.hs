module DewPoint.Geo(Lat(..), Lon(..), toLat, toLon, fromLat, fromLon, toLatLon) where

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

fromLat :: Lat -> Double
fromLat (Lat v) = v

fromLon :: Lon -> Double
fromLon (Lon v) = v

toLon :: Double -> Lon
toLon val = clamp (Lon val)

toLatLon :: Double -> Double -> (Lat, Lon)
toLatLon lat lon = (toLat lat, toLon lon)
