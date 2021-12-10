import ROOT
import os
import argparse
import array
from   RootTools.core.Sample import Sample
from   math import log
from   TMB.Tools.delphesCutInterpreter  import cutInterpreter
from   Analysis.Tools                   import u_float
import Analysis.Tools.syncer            as syncer 
from   Analysis.Tools.cardFileWriter    import cardFileWriter 
from   RootTools.core.standard          import *
import TMB.Tools.helpers                as helpers

argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--logLevel',       action='store',      default='INFO',         nargs='?', choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'],             help="Log level for logging")
argParser.add_argument("--lumi",               action='store',      type=float,             default=137, help='Which lumi?')
argParser.add_argument('--config',             action='store', type=str, default = "ZH_delphes_bkgs", help="config")
argParser.add_argument('--config_module',      action='store', type=str, default = "TMB.BIT.configs", help = "config directory")
argParser.add_argument('--output_directory',   action='store', type=str,   default=os.path.expandvars('/mnt/hephy/cms/$USER/BIT/'))
argParser.add_argument('--small',              action='store_true', help="small?")
argParser.add_argument('--name',               action='store', type=str,   default='v2', help="Name of the training")
argParser.add_argument('--input_directory',    action='store', type=str,   default=os.path.expandvars("/groups/hephy/cms/$USER/BIT/training-ntuple-ZH/MVA-training"))

args = argParser.parse_args()

# Logging
import Analysis.Tools.logger as logger
logger = logger.get_logger(args.logLevel, logFile = None )
import RootTools.core.logger as logger_rt
logger_rt = logger_rt.get_logger(args.logLevel, logFile = None )

# https://stackoverflow.com/questions/21844024/weighted-percentile-using-numpy
def weighted_quantile(values, quantiles, sample_weight=None, 
                      values_sorted=False, old_style=False):
    """ Very close to numpy.percentile, but supports weights.
    NOTE: quantiles should be in [0, 1]!
    :param values: numpy.array with data
    :param quantiles: array-like with many quantiles needed
    :param sample_weight: array-like of the same length as `array`
    :param values_sorted: bool, if True, then will avoid sorting of
        initial array
    :param old_style: if True, will correct output to be consistent
        with numpy.percentile.
    :return: numpy.array with computed quantiles.
    """
    values = np.array(values)
    quantiles = np.array(quantiles)
    if sample_weight is None:
        sample_weight = np.ones(len(values))
    sample_weight = np.array(sample_weight)
    assert np.all(quantiles >= 0) and np.all(quantiles <= 1), \
        'quantiles should be in [0, 1]'

    if not values_sorted:
        sorter = np.argsort(values)
        values = values[sorter]
        sample_weight = sample_weight[sorter]

    weighted_quantiles = np.cumsum(sample_weight) - 0.5 * sample_weight
    if old_style:
        # To be convenient with numpy.percentile
        weighted_quantiles -= weighted_quantiles[0]
        weighted_quantiles /= weighted_quantiles[-1]
    else:
        weighted_quantiles /= np.sum(sample_weight)
    return np.interp(quantiles, weighted_quantiles, values)

def PoissonLL( lam, n ):
    return lam - n*log(lam) + sum( [log(n_) for n_ in range(1,n+1)] )

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

#binning_quantiles = [ .001, .01, .025, .05, .1, .2, .3, .5, .6, .7, .8, .9, .95, .975, .99, .999]
#binning_quantiles = [.01, .025, .05, .1, .2, .3, .5, .6, .7, .8, .9, .95, .975, .99,]
#binning_quantiles = [.05, .1, .2, .3, .5, .6, .7, .8, .9, .95, ]
#binning_quantiles = [ ]
#binning_quantiles = [ .5 ]
binning_quantiles = [ i/200. for i in range(1,200) ] 

# loop over EFT param point

WC_vals = [i/50. for i in range(-50,21)] 
#WC_vals = [i/100. for i in range(-10,11)] 
#WC_vals = [i/20. for i in range(-10,11)] 
#WC_vals = [ -0.1, -0.05, 0, 0.05, 0.1 ]
#WC_vals = [ 0.1 ]

config.training_samples = config.training_samples

# test statistic
qs = {'linear'   :{'color':ROOT.kGreen}, 
      'quadratic':{'color':ROOT.kRed},
      'total'    :{'color':ROOT.kBlack},
    }

for q_name, q in qs.iteritems(): 
    # tGraphs 
    q['nll_tGraph'] = ROOT.TGraph( len(WC_vals) )
    q['nll_tGraph'].SetLineColor( q['color'] )
    q['nll_tGraph'].SetMarkerColor( q['color'] )
    q['nll_tGraph'].SetMarkerStyle( 0 )
    q['nll_tGraph'].SetLineWidth(2)
    q['unbinned_nll_tGraph'] = ROOT.TGraph( len(WC_vals) )
    q['unbinned_nll_tGraph'].SetLineColor( q['color'] )
    q['unbinned_nll_tGraph'].SetMarkerColor( q['color'] )
    q['unbinned_nll_tGraph'].SetMarkerStyle( 0 )
    q['unbinned_nll_tGraph'].SetLineWidth(2)
    # cardFileWriter
    #q['c'] = cardFileWriter.cardFileWriter()
    #q['c'].setPrecision(5)

    # test statistic dicts
    q['val'] = {}
    q['nll'] = {}
    q['poissonPrefitNLL'] = {}

for i_WC_val, WC_val in enumerate(WC_vals):#np.arange(-1,1,.1):

    # get weights and test statistics arrays 
    w_sm    = {}
    w_bsm   = {}
    
    q_lin   = {}
    q_tot   = {}
    q_quad  = {}

    for i_training_sample, training_sample in enumerate(config.training_samples):

        # compute weights
        w_sm  [training_sample.name] = weight_derivatives[training_sample.name][:,0]
        w_bsm [training_sample.name] = weight_derivatives[training_sample.name][:,0] + WC_val*weight_derivatives[training_sample.name][:,i_der_lin]+0.5*WC_val**2*weight_derivatives[training_sample.name][:,i_der_quad] 

        # compute test statistics. 
        qs['total']['val']    [training_sample.name] = bit[flavor][training_sample.name][:,i_lin] + 0.5*WC_val*bit[flavor][training_sample.name][:,i_quad]  
        qs['linear']['val']   [training_sample.name] = bit[flavor][training_sample.name][:,i_lin]  
        qs['quadratic']['val'][training_sample.name] = bit[flavor][training_sample.name][:,i_quad] 

        # make binning from signal sample
        if i_training_sample == 0: 
            for q_name, q in qs.iteritems():
                #unweighted_binning = [-float('inf')]+list(np.quantile( q['val'], binning_quantiles))+[float('inf')]
                #unweighted_binning = helpers.remove_duplicates( unweighted_binning )
                weighted_binning = [-float('inf')]+list(weighted_quantile( q['val'][training_sample.name], binning_quantiles, w_sm[training_sample.name]))+[float('inf')]
                q['binning'] = helpers.remove_duplicates( weighted_binning )
                q['observation'] = [0]*( len(q['binning'])-1 )

                # set up cardFileWriter
                #q['c'].reset()

                #for i_bin in range(len(q["binning"])-1):
                #    q['c'].addBin("Bin%i"%i_bin, [s.name for s in config.training_samples[1:]], "Bin%i"%i_bin)

                # add uncertainties
                #q['c'].addUncertainty('JEC',         'lnN') # correlated
                #q['c'].addUncertainty('JER',         'lnN') # correlated
                #q['c'].addUncertainty('btag_heavy',  'lnN') # uncorrelated, wait for offical recommendation
                #q['c'].addUncertainty('btag_light',  'lnN') # uncorrelated, wait for offical recommendation
                #q['c'].addUncertainty('trigger',     'lnN') # uncorrelated, statistics dominated
                #q['c'].addUncertainty('scale',       'lnN') # correlated.
                #q['c'].addUncertainty('PDF',         'lnN') # correlated.
                #q['c'].addUncertainty('Lumi',        'lnN')

                # place holders for histos
                q['h_SM']      = {}
                q['h_BSM']     = {}
                q['histo_SM']  = {}
                q['histo_BSM'] = {}


        # use linear binning for total test statistic
        qs['total']['binning'] = qs['linear']['binning']

        for q_name, q in qs.iteritems():
            #h_SM    = np.histogram(q, binning[q_name], weights = w_sm*float(lumi)/config.scale_weight)
            #h_BSM   = np.histogram(q, binning[q_name], weights = w_bsm*float(lumi)/config.scale_weight)
            #binning = q['binning']
            q['h_SM'][training_sample.name]      = np.histogram(q['val'][training_sample.name], q['binning'], weights = w_sm[training_sample.name]*float(lumi)/config.scale_weight)
            q['h_BSM'][training_sample.name]     = np.histogram(q['val'][training_sample.name], q['binning'], weights = w_bsm[training_sample.name]*float(lumi)/config.scale_weight)

            q['histo_SM'] [training_sample.name] = make_TH1F(q['h_SM'][training_sample.name])
            q['histo_BSM'][training_sample.name] = make_TH1F(q['h_BSM'][training_sample.name])

            #the MC estimate of the extended negativelog  likelihood (NOT ratio, i.e., SM LL is not subtracted)
    
    for q_name, q in qs.iteritems():
   
        q['nll'][WC_val] = 0.
        for i_b in range(len(q['binning'])-1):
            expectations = [ q['h_BSM'][training_sample.name][0][i_b] for training_sample in config.training_samples ]
            observations = [ q['h_SM'] [training_sample.name][0][i_b] for training_sample in config.training_samples ]
            #print i_b, PoissonLL( sum(expectations), int(round(sum(observations))) )
            q['nll'][WC_val] += PoissonLL( sum(expectations), int(round(sum(observations))) )
        q['nll_tGraph'].SetPoint( i_WC_val, WC_val, q['nll'][WC_val] )
        print "nll",q_name, "WC_val",WC_val,":",q['nll'][WC_val]
        for i_training_sample, training_sample in enumerate(config.training_samples):
            q['histo_SM'] [training_sample.name].style      = styles.lineStyle(colors[i_training_sample], dashed = True )
            q['histo_SM'] [training_sample.name].legendText = training_sample.name#+" (SM)" 
            q['histo_BSM'][training_sample.name].style      = styles.lineStyle(colors[i_training_sample], dashed = False )
            q['histo_BSM'][training_sample.name].legendText = training_sample.name#+" (cHW=%4.3f)"%cHW 

            #if q_name=="total":
            #    print 
            #    print q_name, training_sample.name, q['h_SM'][training_sample.name] 
            #    print 

            # loop over bins, filling the card
            #for i_b in range(len(binning[q_name])-1):
            #for i_b in range(len(q['binning'])-1):
                #expected               = max([0.01, round(q['h_BSM'][training_sample.name][0][i_b],4)])
                #q['observation'][i_b] += q['h_SM'][training_sample.name][0][i_b] 
                #name                   = training_sample.name if i_training_sample>0 else "signal"
                #q['c'].specifyExpectation("Bin%i"%i_b, name, expected )
                #q['c'].specifyUncertainty('JEC',         "Bin%i"%i_b, name, 1.09)
                #q['c'].specifyUncertainty('JER',         "Bin%i"%i_b, name, 1.01)
                #q['c'].specifyUncertainty('btag_heavy',  "Bin%i"%i_b, name, 1.04)
                #q['c'].specifyUncertainty('btag_light',  "Bin%i"%i_b, name, 1.04)
                #q['c'].specifyUncertainty('trigger',     "Bin%i"%i_b, name, 1.01)
                #q['c'].specifyUncertainty('scale',       "Bin%i"%i_b, name, 1.01) 
                #q['c'].specifyUncertainty('PDF',         "Bin%i"%i_b, name, 1.01)
                #q['c'].specifyUncertainty('Lumi',        "Bin%i"%i_b, name, 1.00001)

        plot = Plot.fromHisto(name = "q_%s_%s_%4.3f"%(q_name,WC,WC_val), 
                histos = [[q['histo_SM'] [s.name] for s in config.training_samples], 
                         [ q['histo_BSM'][s.name] for s in config.training_samples],],
                texX = "q_{%s=%4.3f}"%(WC,WC_val) , texY = "Number of events" )

        for log_ in [True]:
            plot_directory_ = os.path.join(plot_directory, ("log" if log_ else "lin"))
            plotting.draw(plot, plot_directory = plot_directory_, ratio = {'histos':[(1,0)], 'texY': 'Ratio'}, logY = log_, logX = False, yRange = (10**-2,"auto"), 
                #yRange = (0, 0.5) if 'Mult' not in var else (0, 15 ),  
                legend = ([0.20,0.75,0.9,0.88],3), copyIndexPHP=True, )

        #for i_b in range(len(q['binning'])-1):
        #    q['c'].specifyObservation("Bin%i"%i_b, int(round(q['observation'][i_b],0)) )

        #q['c'].writeToFile( os.path.join(cardfile_directory, '%s_%s_%4.3f.txt'%(q_name,WC,WC_val) ))
        #q['nll'][WC_val] = q['c'].poissonPrefitNLL() 
        #nll = q['c'].calcNLL()
        #print q['nll'][WC_val], nll 
        #q['nll_tGraph'].SetPoint( i_WC_val, WC_val, q['nll'][WC_val] )

    # compute unbinned LL
    unbinned_nll_tot  = 0
    unbinned_nll_quad = 0
    unbinned_nll_lin  = 0
    total_xsec_bsm    = 0
    total_xsec_sm     = 0
    # scale to lumi and compensate for the scaling in the config that we applied to all the training weights
    rescale           = float(lumi)/config.scale_weight
    for i_training_sample, training_sample in enumerate(config.training_samples):

         unbinned_nll_tot  += rescale*np.sum(np.nan_to_num(w_sm[training_sample.name]*( -np.log(1+WC_val*bit[flavor][training_sample.name][:,i_lin] + 0.5*WC_val**2*bit[flavor][training_sample.name][:,i_quad] ))))
         unbinned_nll_quad += rescale*np.sum(np.nan_to_num(w_sm[training_sample.name]*( -np.log(1+0.5*WC_val**2*bit[flavor][training_sample.name][:,i_quad] ))))
         unbinned_nll_lin  += rescale*np.sum(np.nan_to_num(w_sm[training_sample.name]*( -np.log(1+WC_val*bit[flavor][training_sample.name][:,i_lin] ))))

         total_xsec_bsm    += rescale*sum(w_bsm[training_sample.name])
         total_xsec_sm     += rescale*sum(w_bsm[training_sample.name])

    Poissonian_LL_term      = total_xsec_sm*log( total_xsec_bsm/total_xsec_sm )

    qs['total']     ['unbinned_nll_tGraph'].SetPoint( i_WC_val, WC_val, unbinned_nll_tot - Poissonian_LL_term )
    qs['quadratic'] ['unbinned_nll_tGraph'].SetPoint( i_WC_val, WC_val, unbinned_nll_quad- Poissonian_LL_term )
    qs['linear']    ['unbinned_nll_tGraph'].SetPoint( i_WC_val, WC_val, unbinned_nll_lin - Poissonian_LL_term )

for q_name, q in qs.iteritems():
    for key in [ 'nll_tGraph', 'unbinned_nll_tGraph']:
        min_y = min( [q[key].GetY()[i] for i in range(len(WC_vals))] )
        print "Subtracting likelihood offset for",q_name, key, "of", min_y
        for i in range(len(WC_vals)):
            q[key].SetPoint( i, q[key].GetX()[i], -min_y + q[key].GetY()[i])

for key in  [ 'nll_tGraph', 'unbinned_nll_tGraph']:
    c1 = ROOT.TCanvas()
    ROOT.gStyle.SetOptStat(0)
    c1.SetTitle("")

    l = ROOT.TLegend(0.55, 0.8, 0.8, 0.9)
    l.SetFillStyle(0)
    l.SetShadowColor(ROOT.kWhite)
    l.SetBorderSize(0)

    first = True
    for q_name, q in qs.iteritems(): 
        l.AddEntry( q[key], q_name )
        q[key].Draw("AL" if first else "L")
        q[key].SetTitle("")
        q[key].GetXaxis().SetTitle(WC)
        q[key].GetYaxis().SetTitle("NLL")

        first = False

    l.Draw()

    c1.RedrawAxis()
    c1.Print(os.path.join(plot_directory, key+".png"))
    c1.Print(os.path.join(plot_directory, key+".pdf"))
    c1.Print(os.path.join(plot_directory, key+".root"))

syncer.sync()

