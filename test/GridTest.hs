module GridTest(gridTests) where

import Test.Tasty
import Test.Tasty.HUnit
import Test.Tasty.QuickCheck

import DewPoint.Grid
import DewPoint.Geo

import Data.Array.Repa hiding ( (++) )

instance Arbitrary Lat where
    arbitrary = toLat <$> choose (fromLat minBound, fromLat maxBound)

instance Arbitrary Lon where
    arbitrary = toLon <$> choose (fromLon minBound, fromLon maxBound)

instance Arbitrary Grid where
    arbitrary = gridCreateFixed <$> pos <*> pos
        where pos = getPositive <$> arbitrary

gridU :: Grid
gridU = gridCreateFixed 100 200

minBoundForGridTest :: Assertion
minBoundForGridTest = (0,0)
                      @=? (ixToXY $ gridGetLocationIx gridU (minBound, minBound))


gridSrcCellY_testLowBorder :: Assertion
gridSrcCellY_testLowBorder = let ix = ix2 0 0
                             in ix @=? gridSrcCellY gridU ix (1::Double)

gridSrcCellY_testUpperBorder :: Assertion
gridSrcCellY_testUpperBorder = let ix = ix2 0 (gridMaxY gridU)
                               in ix @=? gridSrcCellY gridU ix (-1::Double)

gridTests :: TestTree
gridTests = testGroup "Grid Tests" [
             testGroup "Unit tests"
                           [ testCase "gridMinBound" minBoundForGridTest
                           , testCase "test lower border gridSrcCellY" gridSrcCellY_testLowBorder
                           , testCase "test upper border gridSrcCellY" gridSrcCellY_testUpperBorder
                           ],

             testGroup "QuickCheck tests"
                           [ testProperty "gridGetXY test" prop_gridGetXY
                           , testProperty "gridSrcCellX test" prop_gridSrcCellX
                           , testProperty "gridSrcCellY test" prop_gridSrcCellY
                           ]
            ]

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
                                                        ix' = gridSrcCellX g ix (1::Double)
                                                    in  ix == gridSrcCellX g ix' (-1.0::Double)

prop_gridSrcCellY :: Grid -> (NonNegative Int) -> (Positive Int) -> Property
prop_gridSrcCellY g (NonNegative x) (Positive y) = x <= gridMaxX g && y < gridMaxY g && y > 0
                                                      ==> let ix = ix2 x y
                                                              ix' = gridSrcCellY g ix (-1::Double)
                                                          in  ix == gridSrcCellY g ix' (1::Double)
