--{-# LANGUAGE OverloadedStrings #-}
module Main where

import Paths_dewpoint (version)
import Data.Version (showVersion)
import Data.Time.Clock(UTCTime(..), secondsToDiffTime, diffUTCTime)
import qualified CommandLineOptions as Opt
import qualified LoadGrid as LG
import MeteoStation
import qualified Plot as Plot
import Options.Applicative((<>), execParser, info, helper, fullDesc, progDesc, header)

import qualified Data.Set as S

import qualified DewPoint.Grid as G
import qualified DewPoint.Geo as Geo
import qualified DewPoint.Physics as Ph
import qualified DewPoint.Simulation as Sim
import qualified Data.Array.Repa as R
import System.FilePath
import Control.Monad(when)
import qualified Data.ByteString as BS
import qualified Data.ByteString.Char8 as BS8
import qualified Data.Text as T
import qualified Data.Text.Encoding as T
import Data.Foldable(toList)

-- | Return version string
versionInfo :: String
versionInfo = "Weather simulation v." ++ showVersion version

validateOpt :: Monad m => Opt.ProgramOpt -> m Opt.ProgramOpt
validateOpt o = if (Opt.outputType o) == Opt.Text && null (Opt.stationListIATACodes o)
                then fail "Must specify at least one station when --output set to \"text\""
                else return o


printOutput :: ( (G.Grid, UTCTime) -> BS.ByteString) -- ^ print function
            -> UTCTime           -- ^ start time
            -> Integer           -- ^ print interval in seconds
            -> (G.Grid, UTCTime) -- ^ (grid, currentTime)
            -> IO (G.Grid, UTCTime)
printOutput pFn t0 dt (g, t) =
    do let secSinceStart = (round . toRational) $ diffUTCTime t t0
       when (secSinceStart `mod` dt == 0) $
            BS8.putStrLn $ pFn (g, t)
       return (g, t)

printStationdData :: Bool -> [MeteoStation] -> (G.Grid, UTCTime) -> BS.ByteString
printStationdData verbose = flip f'
    where f' (g, t) = T.encodeUtf8 . T.intercalate (T.pack "\n")
                      . fmap (T.pack . stationString verbose g t)

plotHeader :: G.Grid -> BS.ByteString
plotHeader g = Plot.encodeHeader "Temperature and U,V winds plot"
               (toList $ Geo.fromLat <$> G.gridLats g)
               (toList $ Geo.fromLon <$> G.gridLons g)

plotFrame :: (G.Grid, UTCTime) -> BS.ByteString
plotFrame (g, t) = Plot.encodeFrame t
                   ( (map Ph.toCelsius .  R.toList . R.transpose . G.gridAirTemp)  g)
                   (R.toList $ R.transpose $ G.gridUWind g)
                   (R.toList $ R.transpose $ G.gridVWind g)


run :: Opt.ProgramOpt -> IO()
run opt = do
  let
      tstart = UTCTime (Opt.startDate opt) 0
      timeInterval = secondsToDiffTime (Opt.interval  opt) * 60
      printInterval = (Opt.outputInterval  opt) * 60
      st = S.fromList (Opt.stationListIATACodes opt)
      dataDir = Opt.dataDir opt
      gridFile = dataDir </> "grid.json"
      beVerbose = Opt.verbosity opt == Opt.Verbose
  maybeGrid <- LG.loadGrid <$> BS.readFile gridFile
  grid <- case maybeGrid of
    (Just a) -> return a
    _        -> fail $ "Unable to load grid from " ++ show gridFile
                       ++ " run scripts/import_data.py to import data"

  stations <- loadStations  dataDir st grid
  let isPlot = Opt.Plot  == Opt.outputType opt
      printFn = if isPlot then plotFrame else printStationdData beVerbose stations
  when isPlot $  (BS8.putStrLn . plotHeader) grid
  loop printFn tstart timeInterval printInterval (grid, tstart)
    where loop st t0 dt dt' (g, t) = Sim.integrate g t dt
                                     >>= printOutput st t0 dt'
                                     >>= loop st t0 dt dt'

main :: IO ()
main = execParser opts >>= validateOpt >>= run
  where
    opts = info (helper <*> Opt.parseOptions)
      ( fullDesc
     <> progDesc "Run weather simulation"
     <> header versionInfo )
