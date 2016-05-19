module CommandLineOptions(ProgramOpt(..), parseOptions, OutputType(..)) where

import Options.Applicative
import Data.List.Split(splitOn)
import Data.Time.Calendar(Day, fromGregorian)
import Data.Time.Format(defaultTimeLocale, parseTimeM, formatTime)

data Verbosity = Normal | Verbose deriving (Show, Eq)

data OutputType = Plot | Text deriving (Show, Eq)

outputTypeToString :: OutputType -> String
outputTypeToString t = case t of
  Plot -> "plot"
  Text -> "text"

data ProgramOpt = ProgramOpt {
   stationListIATACodes :: [String]
  , verbosity  :: Verbosity
  , startDate  :: Day
  , numDays    :: Integer
  , interval   :: Integer
  , outputInterval :: Integer
  , outputType :: OutputType
  , dataDir    :: String
  } deriving (Show)


parseStartDate :: Parser Day
parseStartDate = option  (str >>= parseTimeM True defaultTimeLocale "%F")
                 ( short 'd'
                   <> long "start-date"
                   <> metavar "START-DATE"
                   <> value (fromGregorian 2016 01 01)
                   <> showDefaultWith (formatTime defaultTimeLocale "%F")
                   <> help "Start date in the format yyyy-MM-dd"
                 )

parseNumDays :: Parser Integer
parseNumDays = option auto
               ( short 'n'
              <> long "num-days"
              <> value 30
              <> showDefault
              <> metavar "NUM-DAYS"
              <> help "Number of days to run"
              )

parseDataDir :: Parser String
parseDataDir = option auto
               ( long "data-dir"
              <> value "data"
              <> showDefault
              <> metavar "DATA-DIR"
              <> help "Directory with input data"
              )

parseStepInterval :: Parser Integer
parseStepInterval = option auto
               ( short 'i'
              <> long "interval"
              <> metavar "INTEGRATION-INTERVAL"
              <> value 10
              <> showDefault
              <> help "Simulation step interval in minutes"
              )

parseOutputInterval :: Parser Integer
parseOutputInterval = option auto
               ( short 'p'
              <> long "print-interval"
              <> metavar "PRINT-INTERVAL"
              <> value 180
              <> showDefault
              <> help "Output interval in minutes, should be multiple of integration interval"
              )


parseStationList :: Parser [String]
parseStationList =  filter (not . null) . (splitOn ",")
                    <$> strOption
                    ( short 's'
                      <> long "stations"
                      <> metavar "STATION-LIST"
                      <> value []
                      <> help "Comma separated list of stations as 3-letter IATA code " )

parseOutputType :: Parser OutputType
parseOutputType = option (str >>= parse')
               ( long "output"
                 <> short 'o'
                 <> value Text
                 <> showDefaultWith outputTypeToString
                 <> help "Output type: [plot | text]" )
  where parse' s = case s of
                       "plot" -> return Plot
                       "text" -> return Text
                       _      -> fail "output type on of : [plot | text] expected"

parseOptions :: Parser ProgramOpt
parseOptions = ProgramOpt
               <$> parseStationList
               <*> flag Normal Verbose
               ( long "verbose"
                 <> short 'v'
                 <> help "Enable verbose mode" )
               <*> parseStartDate
               <*> parseNumDays
               <*> parseStepInterval
               <*> parseOutputInterval
               <*> parseOutputType
               <*> parseDataDir
