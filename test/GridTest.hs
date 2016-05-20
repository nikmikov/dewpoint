module GridTest(
                prop_gridGetXY
               , prop_gridSrcCellX
               ) where

import Test.QuickCheck.Gen
import Test.QuickCheck.Arbitrary
import Test.QuickCheck.Modifiers
import Test.QuickCheck.Property

import DewPoint.Grid
import DewPoint.Geo

import Data.Array.Repa hiding ( (++) )

import Debug.Trace

instance Arbitrary Lat where
    arbitrary = toLat <$> choose (fromLat minBound, fromLat maxBound)

instance Arbitrary Lon where
    arbitrary = toLon <$> choose (fromLon minBound, fromLon maxBound)

instance Arbitrary Grid where
    arbitrary = gridCreateFixed <$> pos <*> pos
        where pos = getPositive <$> arbitrary

prop_gridGetXY :: Grid -> (Lat, Lon) -> Property
prop_gridGetXY grid (lat, lon) = True
                                 ==> let ix = gridGetLocationIx grid (lat, lon)
                                         cb = gridGetCellBBox grid ix
                                         (lat1, lon1) = cellBoxTopLeft cb
                                         (lat2, lon2) = cellBoxBottomRight cb
                                     in (lat >= lat1 && lat < lat2) && (lon >= lon1 && lon1 < lon2 )

prop_gridSrcCellX :: Grid -> (NonNegative Int) -> (NonNegative Int) -> Property
prop_gridSrcCellX g (NonNegative x) (NonNegative y) = x <= gridMaxX g && y <= gridMaxY g
                                                ==> let ix = ix2 x y
                                                        ix' = gridSrcCellX g ix 1
                                                    in  ix == gridSrcCellX g ix' (-1)

prop_gridSrcCellY :: Grid -> (NonNegative Int) -> (NonNegative Int) -> Property
prop_gridSrcCellY g (NonNegative x) (NonNegative y) = x <= gridMaxX g && y <= gridMaxY g
                                                      ==> let ix = ix2 x y
                                                              ix' = gridSrcCellY g ix (-1)
                                                          in  ix == gridSrcCellY g ix' 1