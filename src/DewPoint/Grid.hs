module DewPoint.Grid(Grid(..), CellBox(..), Ix,
                     gridCreate, gridCreateFixed,
                     gridSizeXY,
                     gridGetLocationIx,
                     gridGetCellBBox,
                     gridCellLengthMeters,
                     gridSrcCellX, gridSrcCellY,
                     ixToXY, xyToIx
) where

import DewPoint.Geo
import Data.Array.Repa ( (:.)(..) )
import qualified Data.Array.Repa as R
import Data.Sequence(Seq)
import qualified Data.Sequence as S
import qualified Data.Array.Repa.Repr.Vector as R

{-
Grid over earth surface represenation

Y
|______
|_|_|_|
|_|_|_|
|_|_|_|_ X

X axis directed from west to east
Y axis directed from south to north

for each cell C(x,y) with cell box (lat1, lon1, lat2, lon2) and point p(lat, lon)
p is inside the cell if p.lat >= lat1 && p.lat < lat2 && p.lon >= lon1 && p.lon < lon2

-}

data CellBox  = CellBox {
  cellBoxTopLeft :: (Lat, Lon)
  , cellBoxBottomRight :: (Lat, Lon)
  }

type Ix = R.DIM2

type DUArray = R.Array R.U Ix Double

data Grid = Grid {
  gridLats :: Seq Lat
  , gridLons :: Seq Lon
  , gridCellGeoCoord :: R.Array R.V Ix (Lat, Lon)
  , gridSurfaceHeight :: R.Array R.U Ix Int
  , gridShape :: Ix
  , gridAirTemp :: !DUArray
  , gridUWind :: !DUArray
  , gridVWind :: !DUArray
  , gridCloudCover :: !DUArray
  , gridWaterDensity :: !DUArray -- ^ density of water vapor in atmosphere
  } deriving (Show)


ixToXY :: Ix -> (Int, Int)
ixToXY (_ :. x :. y) = (x, y)

-- | create grid from given list of lats,lons
gridCreate :: [Lat] -> [Lon] -> Grid
gridCreate lats lons =
  Grid
  {
    gridLats = lats'
  , gridLons = lons'
  , gridShape = sh
  , gridCellGeoCoord = R.computeVectorS $ R.fromFunction sh (computeCellCenter lats' lons')
  , gridSurfaceHeight = izeros
  , gridAirTemp       = dzeros
  , gridUWind         = dzeros
  , gridVWind         = dzeros
  , gridCloudCover    = dzeros
  , gridWaterDensity  = dzeros
  }
  where lats' = S.sort $ S.fromList lats
        lons' = S.sort $ S.fromList lons
        numX = (S.length lons' - 1)
        numY = (S.length lats' - 1)
        sh = R.Z :. numX :. numY
        dzeros = (R.computeS . R.fromFunction sh . const) 0
        izeros = (R.computeS . R.fromFunction sh . const) 0

-- |create grid given number of cell by X and Y axis
gridCreateFixed :: Int -> Int -> Grid
gridCreateFixed sizeX sizeY = let latLen = fromLat maxBound - fromLat minBound
                                  lonLen = fromLon maxBound - fromLon minBound
                                  stepX = lonLen / fromIntegral (sizeX + 1)
                                  stepY = latLen / fromIntegral (sizeY + 1)
                                  lons = toLon <$> [ fromIntegral x * stepX | x <- [0..sizeX-1]  ]
                                  minLat = fromLat minBound
                                  lats = toLat <$> [ minLat + fromIntegral y * stepY| y <- [0..sizeY-1]  ]
                              in gridCreate (lats ++ maxBound:[]) (lons ++ maxBound:[])

-- | number of cells by X, Y
gridSizeXY :: Grid -> (Int, Int)
gridSizeXY grid = let (_ :. x :. y) = gridShape grid
                  in (x, y)

-- | return cell index on the grid with given geo coordinates
gridGetLocationIx :: Grid -> (Lat, Lon) -> Ix
gridGetLocationIx grid (lat, lon) = R.Z :. findLonX grid lon :. findLatY grid lat


-- | return bounding box of cell
gridGetCellBBox :: Grid -> Ix -> CellBox
gridGetCellBBox grid (_ :. x :. y) =
  CellBox
  {
    cellBoxTopLeft = (gridGetLatAt grid y, gridGetLonAt grid x)
  , cellBoxBottomRight = ( gridGetLatAt grid (succ y), gridGetLonAt grid (succ x))
  }

-- | aproximate length of cell side in meters (assuming cell is a square)
gridCellLengthMeters :: Grid -> Double
gridCellLengthMeters g = let (numX, numY) = gridSizeXY g
                             cellArea = earthSurfaceArea / fromIntegral (numX * numY)
                         in sqrt(cellArea)

xyToIx :: Int -> Int -> Ix
xyToIx x y = (R.Z :. x :. y)

gridMaxX :: Grid -> Int
gridMaxX = (-) 1 . fst . gridSizeXY

gridMaxY :: Grid -> Int
gridMaxY = (-) 1 . snd . gridSizeXY

gridSrcCellX :: Grid -> Ix -> Double -> Ix
gridSrcCellX g (R.Z :. x :. y) u
  | u > 0  = if x > 0 then xyToIx (x - 1) y else xyToIx maxX y
  | u < 0  = if x >= maxX then xyToIx 0 y else xyToIx (x + 1) y
  where maxX = gridMaxX g
gridSrcCellX _ ix _ = ix

gridSrcCellY :: Grid -> Ix -> Double -> Ix
gridSrcCellY g ix@(R.Z :. x :. y) v
  | v > 0  = if y > 0 then  xyToIx x (y - 1) else ix
  | v < 0  = if y < maxY then  xyToIx x (y + 1) else ix
  where maxY = gridMaxY g
gridSrcCellY _ ix _ = ix


----------------------------------------------------------------------------------------------------------------
-- non exported functions

computeCellCenter :: Seq Lat -> Seq Lon -> Ix -> (Lat, Lon)
computeCellCenter lats lons ix = let (x, y) = ixToXY ix
                                     (Lat lat1, Lat lat2) = ( S.index lats y, S.index lats (y + 1) )
                                     (Lon lon1, Lon lon2) = ( S.index lons x, S.index lons (x + 1) )
                                     middle v1 v2 = v1 + (v2 - v1) / 2
                                 in toLatLon (middle lat1 lat2) (middle lon1 lon2)


gridGetLonAt :: Grid -> Int -> Lon
gridGetLonAt grid = S.index (gridLons grid)

gridGetLatAt :: Grid -> Int -> Lat
gridGetLatAt = S.index . gridLats


-- | return index i such as xs[i] <= val < xs[i + 1]
findIdx :: Ord a => a -> Seq a -> Int
findIdx val xs = case  S.findIndexL (> val) xs of
                   (Just i) -> i - 1
                   _        -> length xs - 1

-- | get y coordinate of the point with given longitude
--   !!! uses very inefficient seqscan
findLatY :: Grid -> Lat -> Int
findLatY grid lat = findIdx lat ( gridLats grid )

-- | get x coordinate of the point with given latitude
--   !!! uses very inefficient seqscan
findLonX :: Grid -> Lon -> Int
findLonX grid lon = findIdx lon ( gridLons grid )
