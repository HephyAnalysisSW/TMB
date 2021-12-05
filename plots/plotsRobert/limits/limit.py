import ROOT
import os
import argparse
import array
from RootTools.core.Sample import Sample

from TMB.Tools.delphesCutInterpreter    import cutInterpreter
from Analysis.Tools                     import u_float
import Analysis.Tools.syncer            as syncer 
from Analysis.Tools.cardFileWriter      import cardFileWriter 
from RootTools.core.standard            import *
import TMB.Tools.helpers                as helpers

argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--logLevel',       action='store',      default='INFO',         nargs='?', choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'],             help="Log level for logging")
argParser.add_argument("--lumi",               action='store',      type=float,             default=137, help='Which lumi?')
argParser.add_argument('--config',             action='store', type=str, default = "ZH_delphes_bkgs", help="config")
argParser.add_argument('--config_module',      action='store', type=str, default = "TMB.BIT.configs", help = "config directory")
argParser.add_argument('--output_directory',   action='store', type=str,   default=os.path.expandvars('/mnt/hephy/cms/$USER/BIT/'))
argParser.add_argument('--small',              action='store_true', help="small?")
argParser.add_argument('--name',               action='store', type=str,   default='default', help="Name of the training")
argParser.add_argument('--input_directory',    action='store', type=str,   default=os.path.expandvars("/groups/hephy/cms/$USER/BIT/training-ntuple-ZH/MVA-training"))
argParser.add_argument('--selection',          action='store',      default='singlelep-WHJet-onH')

args = argParser.parse_args()

# Logging
import Analysis.Tools.logger as logger
logger = logger.get_logger(args.logLevel, logFile = None )
import RootTools.core.logger as logger_rt
logger_rt = logger_rt.get_logger(args.logLevel, logFile = None )

# card file writer
c                 = cardFileWriter.cardFileWriter()
c.setPrecision(3)

# MVA configuration
import importlib
configs = importlib.import_module(args.config_module)
config  = getattr( configs, args.config)

if args.small:
    args.name+='_small'

import TMB.Tools.user as user
# directories
plot_directory      = os.path.join( user. plot_directory, 'MVA', args.config, args.name)
cardfile_directory  = os.path.join( args.output_directory, 'cardfiles', args.config, args.name)
# saving
#if not os.path.exists(output_directory):
#    try:
#        os.makedirs(output_directory)
#    except OSError:
#        pass

import uproot
import awkward
import numpy as np
import pandas as pd

if args.small:
    args.name+='_small'

# get the training variable names
mva_variables = [ mva_variable[0] for mva_variable in config.mva_variables]

n_var_flat   = len(mva_variables)

# selection
selectionString = cutInterpreter.cutString(args.selection)

features           = {}
weight_derivatives = {}

bit_flavors = ['bkgs', 'nom']
bit          = {key:{} for key in bit_flavors}
bit_branches = {key:[ 'BIT_%s_%s'%( key, "_".join(der) ) for der in config.bit_derivatives] for key in bit_flavors }
n_der        = len(config.weight_derivative_combinations)
n_bits       = len(config.bit_derivatives)

n_small_samples = 10000 if args.small else None

for i_training_sample, training_sample in enumerate(config.training_samples):
    upfile_name = os.path.join(os.path.expandvars(args.input_directory), args.config, training_sample.name, training_sample.name+'.root')
    logger.info( "Loading upfile %i: %s from %s", i_training_sample, training_sample.name, upfile_name)
    upfile = uproot.open(upfile_name)

    features[training_sample.name]            = upfile["Events"].pandas.df(branches = mva_variables ).values[:n_small_samples]
    weight_derivatives[training_sample.name]  = upfile["Events"].pandas.df(branches = ["weight_derivatives"] ).values.reshape((-1,n_der))[:n_small_samples]
    for flavor in bit_flavors:
        bit[flavor][training_sample.name]     = upfile["Events"].pandas.df(branches = bit_branches[flavor]).values.reshape((-1,n_bits))[:n_small_samples]

flavor     = 'bkgs'
thresholds = range(0,20,3)+[float('inf')]
lumi       = 59.7

WC         = 'cHj3'
#WC         = 'cHW'

# find position of linear and quadratic estimated coefficient 
i_lin   = config.bit_derivatives.index((WC,))
i_quad  = config.bit_derivatives.index((WC, WC))
# find linear and quadratic weights
i_der_lin   = config.weight_derivative_combinations.index((WC,))
i_der_quad  = config.weight_derivative_combinations.index((WC, WC))

colors = [ROOT.kBlack, ROOT.kBlue, ROOT.kRed ]
def make_TH1F( h ):
    # remove infs from thresholds
    vals, thrs = h
    #thrs[0]  = thrs[1] - (thrs[2]-thrs[1])
    #thrs[-1] = thrs[-2] + (thrs[-2]-thrs[-3])
    histo = ROOT.TH1F("h","h",len(vals),0,len(vals))
    for i_v, v in enumerate(vals):
        histo.SetBinContent(i_v+1, v)
    return histo

binning_quantiles = [.0001, .001, .01, .025, .05, .1, .2, .3, .5, .6, .7, .8, .9, .95, .975, .99, .999, .9999]
binning_quantiles = [.01, .025, .05, .1, .2, .3, .5, .6, .7, .8, .9, .95, .975, .99,]

# nll values 
nll = {}
# loop over EFT param point

WC_vals = list(np.arange(-.1,.1,.01))
for WC_val in WC_vals:#np.arange(-1,1,.1):

    # set up card file
    cardFileName      = os.path.join(cardfile_directory, '%s_%3.2f.txt'%(WC,WC_val) )
    c.reset()

    # Make also histos
    histos = {}

    # fill card file & limits
    for i_training_sample, training_sample in enumerate(config.training_samples):

        histos[training_sample.name] = {}

        q_lin   = bit[flavor][training_sample.name][:,i_lin]   
        q_quad  = bit[flavor][training_sample.name][:,i_quad]  

        # compute test statistic. Scale to (r-1)/theta
        #q       = cHW*bit[flavor][training_sample.name][:,i_lin] + 0.5*cHW**2*bit[flavor][training_sample.name][:,i_quad]  
        q       = bit[flavor][training_sample.name][:,i_quad]  

        w_sm    = weight_derivatives[training_sample.name][:,0]
        w_bsm   = weight_derivatives[training_sample.name][:,0] + WC_val*weight_derivatives[training_sample.name][:,i_der_lin]+0.5*WC_val**2*weight_derivatives[training_sample.name][:,i_der_quad] 

        # determine binning from the test statistic distribution in the signal sample
        if i_training_sample==0:
            binning = [-float('inf')]+list(np.quantile( q, binning_quantiles))+[float('inf')]
            binning = helpers.remove_duplicates( binning )
            # zeros where we add up the observation (=SM expectation)
            observation = [0]*( len(binning)-1 )
            # define the bins
            for i_bin in range(len(binning)-1):
                c.addBin("Bin%i"%i_bin, [s.name for s in config.training_samples[1:]], "Bin%i"%i_bin)

            # add uncertainties
            c.addUncertainty('JEC',         'lnN') # correlated
            c.addUncertainty('JER',         'lnN') # correlated
            c.addUncertainty('btag_heavy',  'lnN') # uncorrelated, wait for offical recommendation
            c.addUncertainty('btag_light',  'lnN') # uncorrelated, wait for offical recommendation
            c.addUncertainty('trigger',     'lnN') # uncorrelated, statistics dominated
            c.addUncertainty('scale',       'lnN') # correlated.
            c.addUncertainty('PDF',         'lnN') # correlated.
            c.addUncertainty('Lumi',        'lnN')

        h_SM    = np.histogram(q, binning, weights = w_sm*float(lumi)/config.scale_weight)
        h_BSM   = np.histogram(q, binning, weights = w_bsm*float(lumi)/config.scale_weight)

        histos[training_sample.name]['SM']  = make_TH1F(h_SM)
        histos[training_sample.name]['BSM'] = make_TH1F(h_BSM)

        # loop over bins
        for i_b in range(len(binning)-1):

            expected          = max([0.01, round(h_BSM[0][i_b],3)])
            observation[i_b] += h_SM[0][i_b]
            name = training_sample.name if i_training_sample>0 else "signal"
            c.specifyExpectation("Bin%i"%i_b, name, expected )

            c.specifyUncertainty('JEC',         "Bin%i"%i_b, name, 1.09)
            c.specifyUncertainty('JER',         "Bin%i"%i_b, name, 1.01)
            c.specifyUncertainty('btag_heavy',  "Bin%i"%i_b, name, 1.04)
            c.specifyUncertainty('btag_light',  "Bin%i"%i_b, name, 1.04)
            c.specifyUncertainty('trigger',     "Bin%i"%i_b, name, 1.01)
            c.specifyUncertainty('scale',       "Bin%i"%i_b, name, 1.01) 
            c.specifyUncertainty('PDF',         "Bin%i"%i_b, name, 1.01)
            c.specifyUncertainty('Lumi',        "Bin%i"%i_b, name, 1.023)

        for i_b in range(len(binning)-1):
            c.specifyObservation("Bin%i"%i_b, int(round(observation[i_b],0)) )

    c.writeToFile(cardFileName)
    nll[WC_val] = c.calcNLL()
    print "nll['nll']+nll['nll0']", nll[WC_val]['nll']+nll[WC_val]['nll0']

    for i_training_samples, training_sample in enumerate(config.training_samples):
        histos[training_sample.name]['SM'].style       = styles.lineStyle(colors[i_training_samples], dashed = True )
        histos[training_sample.name]['SM'].legendText  = training_sample.name#+" (SM)" 
        histos[training_sample.name]['BSM'].style      = styles.lineStyle(colors[i_training_samples], dashed = False )
        histos[training_sample.name]['BSM'].legendText = training_sample.name#+" (cHW=%3.2f)"%cHW 

    plot = Plot.fromHisto(name = "q_%s_%3.2f"%(WC,WC_val), 
            histos = [[histos[s.name]['SM'] for s in  config.training_samples ], 
                     [ histos[s.name]['BSM'] for s in config.training_samples],],
            texX = "q_{%s=%3.2f}"%(WC,WC_val) , texY = "Number of events" )

    for log in [True]:
        plot_directory_ = os.path.join(plot_directory, ("log" if log else "lin"))
        plotting.draw(plot, plot_directory = plot_directory_, ratio = {'histos':[(1,0)], 'texY': 'Ratio'}, logY = log, logX = False, yRange = (10**-2,"auto"), 
            #yRange = (0, 0.5) if 'Mult' not in var else (0, 15 ),  
            legend = ([0.20,0.75,0.9,0.88],3), copyIndexPHP=True,
        )

min_NLL = min([nll[WC_val]['nll']+nll[WC_val]['nll0'] for WC_val in WC_vals] )
t = ROOT.TGraph( len(WC_vals), 
             array.array('d',WC_vals), 
             array.array('d',[ nll[WC_val]['nll']+nll[WC_val]['nll0'] - min_NLL for WC_val in WC_vals] )
    )

t.SetLineColor(ROOT.kBlue)
t.SetLineWidth(2)
#t.SetLineStyle(7)

l = ROOT.TLegend(0.55, 0.8, 0.8, 0.9)
l.SetFillStyle(0)
l.SetShadowColor(ROOT.kWhite)
l.SetBorderSize(0)

l.AddEntry( t, "q" )
c1 = ROOT.TCanvas()
t.Draw("AL")
t.SetTitle("")
t.GetXaxis().SetTitle(WC)
t.GetYaxis().SetTitle("NLL")

t.Draw("L")
l.Draw()

ROOT.gStyle.SetOptStat(0)
c1.SetTitle("")
c1.RedrawAxis()
c1.Print(os.path.join(plot_directory, "roc.png"))
c1.Print(os.path.join(plot_directory, "roc.pdf"))
c1.Print(os.path.join(plot_directory, "roc.root"))

syncer.sync()

