import Test.Tasty
import Test.Tasty.HUnit
import Test.Tasty.QuickCheck
import qualified Test.Tasty.SmallCheck as SC
import Test.QuickCheck.Gen

import DewPoint.Grid
import DewPoint.Geo

instance Arbitrary Lat where
  arbitrary = toLat <$> choose (fromLat minBound, fromLat maxBound)

instance Arbitrary Lon where
  arbitrary = toLon <$> choose (fromLon minBound, fromLon maxBound)

instance Arbitrary Grid where
  arbitrary = gridCreateFixed <$> pos <*> pos
    where pos = getPositive <$> arbitrary

negation :: Integer -> Bool
negation x = abs (x^2) >= x


prop_gridGetXY :: Grid -> (Lat, Lon) -> Property
prop_gridGetXY grid (lat, lon) = True ==> let ix = gridGetLocationIx grid (lat, lon)
                                              cb = gridGetCellBBox grid ix
                                              (lat1, lon1) = cellBoxTopLeft cb
                                              (lat2, lon2) = cellBoxBottomRight cb
                                          in (lat >= lat1 && lat < lat2) && (lon >= lon1 && lon1 < lon2  )

suite :: TestTree
suite = testGroup "Test Suite" [
    testGroup "Units"
      [ testCase "gridGeXY"  $ (0,0) @=? (ixToXY $ gridGetLocationIx (gridCreateFixed 100 200) (minBound, minBound))
      ],

    testGroup "QuickCheck tests"
      [ testProperty "Quickcheck test" prop_gridGetXY
      ],

    testGroup "SmallCheck tests"
      [ SC.testProperty "Negation" negation
      ]
  ]




main :: IO ()
main = defaultMain suite
