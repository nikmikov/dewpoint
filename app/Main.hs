module Main where

import Paths_weathersim (version)
import Data.Version (showVersion)
import System.Environment (getArgs)
import CommandLineOptions
import Options.Applicative((<>), execParser, info, helper, fullDesc, progDesc, header)


-- | Return version string
versionInfo :: String
versionInfo = "Weather simulation v." ++ showVersion version

run :: ProgramOpt -> IO()
run opt = print opt

main :: IO ()
main = execParser opts >>= run
  where
    opts = info (helper <*> parseOptions)
      ( fullDesc
     <> progDesc "Print a greeting for TARGET"
     <> header versionInfo )
