module DewPoint.Grid(Grid(..), CellBox(..), Ix,
                     gridCreate, gridCreateFixed,
                     gridSizeXY,
                     gridGetLocationIx,
                     gridGetCellBBox,
                     gridCellLengthMeters,
                     gridSrcCellX, gridSrcCellY,
                     ixToXY,
                     gridMaxX, gridMaxY,
                     gridCellRight, gridCellLeft, gridCellUp, gridCellDown
) where

import DewPoint.Geo
import Data.Array.Repa ( (:.)(..) )
import qualified Data.Array.Repa as R
import Data.Sequence(Seq)
import qualified Data.Sequence as S
import qualified Data.Array.Repa.Repr.Vector as R

--import Debug.Trace

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

type DUArray = R.Array R.U Ix Float

data Grid = Grid {
      gridLats :: Seq Lat
    , gridLons :: Seq Lon
    , gridCellGeoCoord :: R.Array R.V Ix (Lat, Lon)
    , gridShape :: Ix
    , gridSurfaceHeight :: !DUArray
    , gridAirTemp :: !DUArray
    , gridUWind :: !DUArray
    , gridVWind :: !DUArray
    , gridCondensedWater :: !DUArray  -- ^ density of condensed water in atmosphere (clouds)
    , gridVaporWater :: !DUArray      -- ^ density of water vapor in atmosphere
    , gridWaterDropletsSz :: !DUArray -- ^ avg size of water droplet in the gridCondensedWater
                                      --   when droplet become big enough it's going to rain
    , gridEvapCoeff :: !DUArray
    , gridPrecipationTriggered :: R.Array R.U Ix Bool
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
  , gridSurfaceHeight   = fill' 0
  , gridAirTemp         = fill' 0
  , gridUWind           = fill' 0
  , gridVWind           = fill' 0
  , gridCondensedWater  = fill' 0
  , gridVaporWater      = fill' 0
  , gridEvapCoeff       = fill' 0
  , gridWaterDropletsSz = fill' 0
  , gridPrecipationTriggered = (R.computeS . R.fromFunction sh . const) False
  }
  where lats' = S.sort $ S.fromList lats
        lons' = S.sort $ S.fromList lons
        numX = S.length lons' - 1
        numY = S.length lats' - 1
        sh = R.Z :. numX :. numY
        fill' = R.computeS . R.fromFunction sh . const

-- |create grid given number of cell by X and Y axis
gridCreateFixed :: Int -> Int -> Grid
gridCreateFixed sizeX sizeY = let latLen = fromLat maxBound - fromLat minBound
                                  lonLen = fromLon maxBound - fromLon minBound
                                  stepX = lonLen / fromIntegral (sizeX + 1)
                                  stepY = latLen / fromIntegral (sizeY + 1)
                                  lons = toLon <$> [ fromIntegral x * stepX | x <- [0..sizeX-1]  ]
                                  minLat = fromLat minBound
                                  lats = toLat <$> [ minLat + fromIntegral y * stepY |
                                                     y <- [0..sizeY-1]  ]
                              in gridCreate (lats ++ [maxBound]) (lons ++ [maxBound])

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
gridCellLengthMeters :: Floating a => Grid -> a
gridCellLengthMeters g = let (numX, numY) = gridSizeXY g
                             cellArea = earthSurfaceArea / fromIntegral (numX * numY)
                         in sqrt cellArea

gridMaxX :: Grid -> Int
gridMaxX = flip (-) 1 . fst . gridSizeXY

gridMaxY :: Grid -> Int
gridMaxY = flip (-) 1 . snd . gridSizeXY

gridSrcCellX :: (Ord a, Num a) => Grid -> Ix -> a -> Ix
gridSrcCellX g ix u
  | u > 0  = gridCellLeft g ix
  | u < 0  = gridCellRight g ix
gridSrcCellX _ ix _ = ix

gridSrcCellY :: (Ord a, Num a) => Grid -> Ix -> a -> Ix
gridSrcCellY g ix v
  | v > 0  = gridCellDown g ix
  | v < 0  = gridCellUp g ix
gridSrcCellY _ ix _ = ix

-- | get cell on the right from the given index
gridCellRight :: Grid -> Ix -> Ix
gridCellRight g (R.Z :. x :. y) = if x == gridMaxX g then R.ix2 0 y else R.ix2 (x + 1) y

-- | get cell on the left from the given index
gridCellLeft :: Grid -> Ix -> Ix
gridCellLeft g (R.Z :. x :. y) = if x == 0 then R.ix2 (gridMaxX g) y else R.ix2 (x - 1) y

-- | get cell on the top from the given index (positive direction Y axis)
gridCellUp :: Grid -> Ix -> Ix
gridCellUp g ix@(R.Z :. x :. y) = if y == gridMaxY g then ix else R.ix2 x (y + 1)

-- | get cell on the bottom from the given index (negative direction Y axis)
gridCellDown :: Grid -> Ix -> Ix
gridCellDown _ ix@(R.Z :. x :. y) = if y == 0 then ix else R.ix2 x (y - 1)

-------------------------------------------------------------------------------------------
-- non exported functions

computeCellCenter :: Seq Lat -> Seq Lon -> Ix -> (Lat, Lon)
computeCellCenter lats lons ix =
    let (x, y) = ixToXY ix
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
