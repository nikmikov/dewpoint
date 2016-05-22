module AdvectionTest(advectionTests)
    where

import Test.Tasty
import Test.Tasty.HUnit
import Test.Tasty.QuickCheck

import DewPoint.Advection


advectionTestZeroes :: Assertion
advectionTestZeroes = (0,0) @=? advectionCoefficients 0 0 0 (0::Double)

advectionTestSimple :: Assertion
advectionTestSimple = (0.11, 0.57) @=? advectionCoefficients 2 6 10 (1::Double)

prop_advectionCoeffSumLessIsThan1 ::
    Double -> Double -> NonNegative Double -> NonNegative Double -> Property
prop_advectionCoeffSumLessIsThan1 u v (NonNegative l) (NonNegative dt)
    = True ==> let (x, y) = advectionCoefficients u v l dt
               in (x + y) <= 1

prop_advectionCoeffPositiveUZeroV ::
    NonZero Double -> NonNegative Double -> NonNegative Double -> Property
prop_advectionCoeffPositiveUZeroV (NonZero u) (NonNegative l) (NonNegative dt)
    = abs u > 1e-10 && l > 1 && dt > 1
      ==> let (x, y) = advectionCoefficients u 0 l dt
          in abs x > 0 && y == 0

prop_advectionCoeffZeroUPositiveV ::
    NonZero Double -> NonNegative Double -> NonNegative Double -> Property
prop_advectionCoeffZeroUPositiveV (NonZero v) (NonNegative l) (NonNegative dt)
    = abs v > 1e-10 && l > 1 && dt > 1
      ==> let (x, y) = advectionCoefficients 0 v l dt
          in abs y > 0 && x == 0

prop_advectionCoeffLargeUV ::
    Double -> Double -> Positive Double -> Positive Double -> Property
prop_advectionCoeffLargeUV u v (Positive l) (Positive dt)
    =  abs u * dt > l && abs v * dt > l
       ==> let (x, y) = advectionCoefficients u v l dt
               in (x + y) == 1

prop_advectionCoeffEqualUV ::
    NonZero Double -> Positive Double -> Positive Double -> Property
prop_advectionCoeffEqualUV  (NonZero uv)  (Positive l) (Positive dt) =
    True ==> let (x, y) = advectionCoefficients uv uv l dt
             in x == y


advectionTests :: TestTree
advectionTests = testGroup "Advection Tests"
                 [
                  testGroup "Unit tests"
                                [ testCase "Zero params" advectionTestZeroes
                                , testCase "Simple case " advectionTestSimple
                                ],
                  testGroup "QuickCheck Tests"
                             [ testProperty "Sum of advection coefficients (X,Y) is less than 1"
                                            prop_advectionCoeffSumLessIsThan1
                             , testProperty "Non-zero U-wind and Zero V-wind "
                                            prop_advectionCoeffPositiveUZeroV
                             , testProperty "Zero U-wind and Non-Zero V-wind "
                                            prop_advectionCoeffZeroUPositiveV
                             , testProperty "Large values of U-wind  V-wind "
                                            prop_advectionCoeffLargeUV
                             , testProperty "Equal values of U-wind and  V-wind "
                                            prop_advectionCoeffEqualUV

                             ]
                 ]
