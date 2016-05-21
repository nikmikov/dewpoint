{-# LANGUAGE OverloadedStrings #-}
module Plot(encodeHeader, encodeFrame) where

import Data.Aeson
import Data.ByteString
import Data.ByteString.Lazy(toStrict)
import Data.Time.Clock

--  output plot data to stdout in .json
--  see scripts/plot_data.py for format description

encodeHeader :: String -> [Double] -> [Double] -> ByteString
encodeHeader title lats lons = (toStrict . encode . object)  ["title" .= title,
                                                              "lats" .= lats,
                                                              "lons" .= lons ]

encodeFrame :: (Floating a, ToJSON a) => UTCTime -> [a] -> [a] -> [a] -> ByteString
encodeFrame time temp uwind vwind = (toStrict . encode . object) ["frame" .= time,
                                                                  "temp" .= temp,
                                                                  "vwind" .= vwind,
                                                                  "uwind" .= uwind ]
