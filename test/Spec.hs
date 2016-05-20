import Test.Tasty
import Test.Tasty.HUnit
import Test.Tasty.QuickCheck
import qualified Test.Tasty.SmallCheck as SC

import DewPoint.Grid
import DewPoint.Geo

import GridTest


negation :: Integer -> Bool
negation x = abs (x^2) >= x


suite :: TestTree
suite = testGroup "Test Suite" [
    testGroup "Units"
      [ testCase "gridGeXY"  $ (0,0) @=? (ixToXY $ gridGetLocationIx (gridCreateFixed 100 200) (minBound, minBound))
      ],

    testGroup "QuickCheck tests"
      [ testProperty "Quickcheck test 1" prop_gridGetXY
      , testProperty "Quickcheck test 2" prop_gridSrcCellX
      ],

    testGroup "SmallCheck tests"
      [ SC.testProperty "Negation" negation
      ]
  ]




main :: IO ()
main = defaultMain suite
