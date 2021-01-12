#!/usr/bin/env python
''' Analysis script for standard plots
'''

# Standard imports and batch mode
import ROOT, os
#ROOT.gROOT.SetBatch(True)
from math                           import sqrt, cos, sin, pi, isnan, sinh, cosh, log

# Analysis
import Analysis.Tools.syncer        as syncer

# RootTools
from RootTools.core.standard        import *

# TMB
from TMB.Tools.user                 import plot_directory
from TMB.Tools.helpers              import deltaPhi, getCollection, deltaR, mZ
from TMB.Tools.WeightInfo           import WeightInfo
from TMB.Tools.cutInterpreter       import cutInterpreter
from TMB.Samples.color              import color
# Arguments
import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--logLevel',           action='store',      default='INFO',          nargs='?', choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'], help="Log level for logging")
argParser.add_argument('--plot_directory',     action='store',      default='gen_v5')
argParser.add_argument('--selection',          action='store',      default=None)
argParser.add_argument('--samples',            action='store',      nargs = "*", default='WZ')
argParser.add_argument('--WC',                 action='store',      default='cWWW')
argParser.add_argument('--WCval',              action='store',      type=float, default=0)
argParser.add_argument('--small',                                   action='store_true',     help='Run only on a small subset of the data?')
args = argParser.parse_args()

# Logger
import TMB.Tools.logger as logger
import RootTools.core.logger as logger_rt
logger    = logger.get_logger(   args.logLevel, logFile = None)
logger_rt = logger_rt.get_logger(args.logLevel, logFile = None)

from TMB.Samples.gen_pp import *
samples = map( eval, args.samples )

if args.small: 
    for sample in samples:
        sample.reduceFiles( to = 1 ) 
for sample in samples:
    sample.color = getattr( color, sample.name )

ws = map( lambda s: WeightInfo(s.reweight_pkl), samples )

for i_w, w in enumerate(ws):
    w.set_order(2)

fisher_plots  = [ samples[i_w].get1DHistoFromDraw(
                  "TMath::Log10("+w.get_fisher_weight_string(args.WC, args.WC, **{args.WC:args.WCval})+")", 
                  [270, -25,2], 
                  selectionString = cutInterpreter.cutString(args.selection) if args.selection is not None else "(1)", 
                  weightString = "137*lumiweight1fb" 
                ) for i_w, w in enumerate(ws) ] 

def accumulate( histo, log_x_base = None):
    res = histo.Clone()
    res.Reset()
    for i_bin in reversed(range(1,histo.GetNbinsX()+1)):
        res.SetBinContent( i_bin, res.GetBinContent(i_bin+1) + histo.GetBinContent(i_bin)*10**histo.GetBinLowEdge(i_bin) )
    return res

fisher_plots_cumulative = map( accumulate, fisher_plots )

# determine quantiles
quantiles  = [ 0.1, 0.33, 0.5, 0.66, 0.9]
thresholds = {}
for i_h, h in enumerate(fisher_plots_cumulative):
    max_val         = h.GetBinContent(1)
    quantile_vals   = [ max_val*val for val in quantiles ]
    thresholds[i_h] = {}
    for i_bin in range(1,h.GetNbinsX()+1):
        for i_quantile_val, quantile_val in enumerate(quantile_vals):
            if h.GetBinContent( i_bin ) < quantile_val and not thresholds[i_h].has_key(quantiles[i_quantile_val]):
                thresholds[i_h][quantiles[i_quantile_val]] = 10**h.GetBinLowEdge( i_bin )

print "WC %s val %f" % (args.WC, args.WCval )
for i_sample, sample in enumerate( samples ):
    print "Thresholds for", sample.name, "at quantiles", quantiles, "are", map( lambda q: thresholds[i_sample][q] if thresholds[i_sample].has_key(q) else None, reversed(quantiles) ) 

for i_sample, sample in enumerate(samples):
    fisher_plots[i_sample].legendText            = sample.name
    fisher_plots[i_sample].style                 = styles.lineStyle( sample.color ) 
    fisher_plots_cumulative[i_sample].legendText = sample.name+" (acc.)"
    fisher_plots_cumulative[i_sample].style      = styles.lineStyle( sample.color, dashed = True) 

with open('missing.sh','a') as file_object:
    for i_sample, sample in enumerate(samples):
        if all( thresholds[i_sample].has_key(q) for q in quantiles):
            s_str = "python skim_plots.py --plot_directory {plot_directory} --WC {WC} --WCval {WCval} --FI_thresholds {FI_thresholds} --sample {sample}".format(
                plot_directory  = args.plot_directory,
                WC              = args.WC,
                WCval           = args.WCval,
                sample          = sample.name,
                selection       = args.selection,
                FI_thresholds   = " ".join( map( lambda q: str(thresholds[i_sample][q]), reversed(quantiles) )),
            )
            if args.selection is not None:
                s_str += "--selection "+args.selection
            if args.small:
                s_str+=" --small"
            file_object.write( s_str+'\n' )

#ipython -i skim_plots.py -- --plot_directory gen_v3_v2 --WCval 1 --FI_thresholds  0.07943282347242854 0.5011872336272756 0.6309573444801942 1.0  --WC cWWW --objects W g   --sample WA   --WG --selection nW1p-nG1p --small

p = Plot.fromHisto( 
    name = "fisher_%s_%3.2f"%(args.WC, args.WCval), 
    histos = [[h] for h in fisher_plots + fisher_plots_cumulative], 
    texX = "log_{10}(I_{F})" , texY = "" )

plotting.draw(p, 
    plot_directory = os.path.join(plot_directory, args.plot_directory+('_small' if args.small else '')), 
    ratio = None, logY = True, logX = False, yRange = (10**-3, "auto"),
    legend = ([0.2, 0.75, 0.8, 0.9],2)
    ) 
