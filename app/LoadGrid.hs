{-# LANGUAGE OverloadedStrings #-}
module LoadGrid(loadGrid) where

import Data.Aeson
import DewPoint.Grid
import DewPoint.Geo
import Data.ByteString.Lazy(fromStrict)
import Data.ByteString
import qualified Data.Array.Repa as R
import qualified DewPoint.Physics as Ph
import Control.Monad(mzero)

data GridJson = GridJson [Double] [Double] [Double] [Int]

instance FromJSON GridJson where
    parseJSON (Object v) =  GridJson <$>
                            v .: "lats" <*>
                            v .: "lons" <*>
                            v .: "height" <*>
                            v .: "land_sea_mask"
    parseJSON _ = mzero


toGrid :: GridJson -> Grid
toGrid (GridJson lats lons height seamask) =
    g {
        gridAirTemp = (R.computeS . R.fromFunction sh . const) Ph.avgSurfaceTemp
      , gridEvapCoeff = evapCoeff
      , gridSurfaceHeight = surfaceHeight
    }
    where g = gridCreate (toLat <$> lats') (toLon <$> lons' )
          appendHead (x:xs) = (x:xs) ++ [x]
          appendHead _ = []
          (lats', lons') = (appendHead  lats, appendHead lons)
          sh = gridShape g
          ec v = if v == 1 then 1 else 0.3
          evapCoeff = (R.computeS . R.map ec . R.fromListUnboxed sh) seamask
          height' = fromRational . toRational . min 0
          surfaceHeight = (R.computeS . R.map height' . R.fromListUnboxed sh) height

loadGrid :: ByteString -> Maybe Grid
loadGrid = fmap toGrid . decode . fromStrict
