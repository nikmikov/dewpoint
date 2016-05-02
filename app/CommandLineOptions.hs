module CommandLineOptions(ProgramOpt, parseOptions) where

import Options.Applicative
import Data.List.Split(splitOn)

data Verbosity = Normal | Verbose deriving (Show)

data ProgramOpt = ProgramOpt {
  gridBoxSizeKilometers :: Int
  , stepIntervalMinutes :: Int
  , stationListIATACodes :: [String]
  , verbosity :: Verbosity
  } deriving (Show)


parseGridBoxSize :: Parser Int
parseGridBoxSize = option auto
               ( short 'c'
              <> long "gridbox-size"
              <> metavar "SIZE-KM"
              <> help "Simulation cell size in kilometers [1..100]"
              )

parseStepInterval :: Parser Int
parseStepInterval = option auto
               ( short 'i'
              <> long "interval"
              <> metavar "INTERVAL"
              <> help "Simulation step interval in minutes"
              )

parseStationList :: Parser [String]
parseStationList =  (splitOn ",")
                    <$> strOption
                    ( short 's'
                      <> long "stations"
                      <> metavar "STATION-LIST"
                      <> help "Comma separated list of stations as 3-letter IATA code " )


parseOptions :: Parser ProgramOpt
parseOptions = ProgramOpt <$> parseGridBoxSize
               <*> parseStepInterval
               <*> parseStationList
               <*> flag Normal Verbose
               ( long "verbose"
                 <> short 'v'
                 <> help "Enable verbose mode" )
