#!/usr/bin/env python
''' Compare combinatorical and nominal b-tagging 
'''

# Standard imports and batch mode
import ROOT, os
from math import pi

# Analysis
import Analysis.Tools.syncer        as syncer

# RootTools
from RootTools.core.standard        import *

# TMB
from TMB.Tools.user                 import plot_directory
from TMB.Tools.delphesCutInterpreter import cutInterpreter

# Arguments
import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--logLevel',           action='store',      default='INFO',          nargs='?', choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'], help="Log level for logging")
argParser.add_argument('--plot_directory',     action='store',      default='bTag-comp')
#argParser.add_argument('--selection',          action='store',      default='singlelep-WHJet-onH')
#argParser.add_argument('--signal',             action='store',      default='WH', choices = ['WH', 'ZH'])
argParser.add_argument('--small',                                   action='store_true',     help='Run only on a small subset of the data?')
args = argParser.parse_args()

# Logger'singlelep-WHJet' if sample.name=='WH' else 'dilep-ZHJet-onZ'
import TMB.Tools.logger as logger
import RootTools.core.logger as logger_rt
logger    = logger.get_logger(   args.logLevel, logFile = None)
logger_rt = logger_rt.get_logger(args.logLevel, logFile = None)

plot_directory = os.path.join(plot_directory, args.plot_directory )
if args.small: plot_directory += "_small"

# Import samples
import TMB.Samples.pp_gen_v10 as samples

#sample_nom  = samples.DYBBJets
#sample_comb = samples.DYBBJets_comb
sample_nom  = samples.DYJets
sample_comb = samples.DYJets_comb
#sample_nom  = samples.DYJets_HT
#sample_comb = samples.DYJets_HT_comb

for sample in [sample_nom, sample_comb]:
    #if selectionString != "":
    #    sample.addSelectionString( selectionString )
    if args.small:
        sample.reduceFiles( to = 20 )

yield_nom   = sample_nom .getYieldFromDraw( "nBTag>=2" )
yield_comb  = sample_comb.getYieldFromDraw( "(1)", "combinatoricalBTagWeight2b" )

selection = "dilep-ZHJet-onZ-onH"
selectionString = cutInterpreter.cutString(selection)
for var, binning in [
        #("recoZ_pt",  [20,0,500]),
        #("nBTag",     [5,0,5]),
        ("ZH_Theta", [10,-2*pi,2*pi]),
    ]:

    histo_nom   = sample_nom .get1DHistoFromDraw( var, binning, selectionString+"&&nBTag>=2", "137*lumiweight1fb")
    histo_comb  = sample_comb.get1DHistoFromDraw( var, binning, selectionString, "combinatoricalBTagWeight2b*lumiweight1fb*137" )

    histo_nom .style = styles.errorStyle(ROOT.kBlue)
    histo_comb.style = styles.errorStyle(ROOT.kRed)

    plot = Plot.fromHisto(var, [[histo_comb], [histo_nom] ], texX = var)
    plotting.draw( plot, plot_directory = os.path.join( plot_directory, sample_nom.name+'-'+selection), copyIndexPHP = True)

syncer.sync()
