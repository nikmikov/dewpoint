import Test.Tasty

import GridTest
import AdvectionTest

suite :: TestTree
suite = testGroup "Test Suite" [
         advectionTests,
         gridTests
  ]

main :: IO ()
main = defaultMain suite
