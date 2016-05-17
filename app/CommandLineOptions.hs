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
  , outputType :: OutputType
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

parseStepInterval :: Parser Integer
parseStepInterval = option auto
               ( short 'i'
              <> long "interval"
              <> metavar "INTERVAL"
              <> value 10
              <> showDefault
              <> help "Simulation step interval in minutes"
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
               <*> parseOutputType
