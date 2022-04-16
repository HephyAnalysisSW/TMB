#!/usr/bin/env python
''' Analysis script for standard plots
'''

# Standard imports and batch mode
import ROOT, os, itertools
#ROOT.gROOT.SetBatch(True)
ROOT.gROOT.SetBatch(True)
c1 = ROOT.TCanvas() # do this to avoid version conflict in png.h with keras import ...
c1.Draw()
c1.Print('/tmp/delete.png')

import numpy as np
import copy
import operator
import random
from math                           import sqrt, cos, sin, pi, isnan, sinh, cosh, log, copysign
import array

# Analysis
import Analysis.Tools.syncer        as syncer
from   Analysis.Tools.WeightInfo    import WeightInfo
from   Analysis.Tools.helpers       import deltaPhi, deltaR, getObjDict

# RootTools
from RootTools.core.standard        import *

# TMB
import TMB.Tools.user               as user
from TMB.Tools.helpers              import deltaPhi, getCollection, deltaR, mZ
from TMB.Tools.delphesCutInterpreter import cutInterpreter
import TMB.Tools.VV_angles          as VV_angles
from TMB.Tools.genObjectSelection   import isGoodGenJet
import pickle

# Arguments
import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--logLevel',           action='store',      default='INFO',          nargs='?', choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'], help="Log level for logging")
argParser.add_argument('--plot_directory',     action='store',      default='delphes')
argParser.add_argument('--selection',          action='store',      default='singlelep-WHJet-onH')
argParser.add_argument('--DYsample',           action='store',      default='DYBBJets')
argParser.add_argument('--signal',             action='store',      default='WH', choices = ['WH', 'ZH', 'WH_nlo', 'ZH_nlo'])
argParser.add_argument('--small',                                   action='store_true',     help='Run only on a small subset of the data?')
argParser.add_argument('--overwrite',                               action='store_true',     help='Overwrite?')
argParser.add_argument('--show_derivatives',                        action='store_true',     help='Show also the derivatives?')
argParser.add_argument('--sign_reweight',                           action='store_true',     help='Apply sign(sin(2theta)sin(2Theta)) as weight ')
#argParser.add_argument('--combinatoricalBTags', action='store_true',   help="BTags combinatorical?")

argParser.add_argument("--nBins",              action="store",      type=int, default=30,                 help="Number of bins in each dimension")
argParser.add_argument("--nBinsTestStat",      action="store",      type=int, default=20,                 help="Number of bins for the test statistic")
argParser.add_argument("--WCs",                action="store",      nargs='*', default=["cHj3", -.08, .08, "cHW", 0,.5],                 help="Wilson coefficients")


args = argParser.parse_args()

import TMB.Tools.logger as logger
import RootTools.core.logger as logger_rt
logger    = logger.get_logger(   args.logLevel, logFile = None)
logger_rt = logger_rt.get_logger(args.logLevel, logFile = None)

plot_directory      = os.path.join(user.plot_directory, args.plot_directory,  args.signal )
results_directory   = os.path.join( user.results_directory, args.plot_directory,  args.signal ) 

if args.small: 
    plot_directory += "_small"
    results_directory += "_small"
#if args.combinatoricalBTags: 
#    plot_directory += "_combWeights"
#    results_directory += "_combWeights"

# Import samples
import TMB.Samples.pp_gen_v10 as samples
    
signal = getattr( samples, args.signal)
 
# WeightInfo
signal.weightInfo = WeightInfo(signal.reweight_pkl)
signal.weightInfo.set_order(2)
signal.read_variables = [VectorTreeVariable.fromString( "p[C/F]", nMax=200 )]

tex = {'cpW':'C_{HW}', 'cHW':'C_{HW}', 'cpq3i':'C_{HQ}^{(3)}', 'cHWtil':'C_{H#tilde{W}}', 'cHj3':'C_{HQ}^{(3)}'}

if signal.name.endswith("_nlo"):
    eft_configs = [
        {'color':ROOT.kBlack,       'param':{}, 'tex':"SM"},
        {'color':ROOT.kGreen-4,     'param':{'cpW':-1},  'tex':"c_{HW}=-1"},
        {'color':ROOT.kGreen+2,     'param':{'cpW':1},  'tex':"c_{HW}=1"},
        {'color':ROOT.kBlue-4,      'param':{'cpq3i':-.1}, 'tex':"c_{HQ}^{(3)}=-0.1"},
        {'color':ROOT.kBlue+2,      'param':{'cpq3i':.1},  'tex':"c_{HQ}^{(3)}=0.1"},
        ]
else:
    eft_configs = [
        {'color':ROOT.kBlack,       'param':{}, 'tex':"SM"},
        {'color':ROOT.kMagenta-4,   'param':{'cHWtil':-1},  'tex':"c_{H#tilde{W}}=-1"},
        {'color':ROOT.kMagenta+2,   'param':{'cHWtil':1},   'tex':"c_{H#tilde{W}}=1"},
        {'color':ROOT.kGreen-4,     'param':{'cHW':-1},  'tex':"c_{HW}=-1"},
        {'color':ROOT.kGreen+2,     'param':{'cHW':1},  'tex':"c_{HW}=1"},
        {'color':ROOT.kBlue-4,      'param':{'cHj3':-.1}, 'tex':"c_{HQ}^{(3)}=-0.1"},
        {'color':ROOT.kBlue+2,      'param':{'cHj3':.1},  'tex':"c_{HQ}^{(3)}=0.1"},
        ]

for eft in eft_configs:
    eft['func'] = signal.weightInfo.get_weight_func(**eft['param']) 
    eft['name'] = "_".join( ["signal"] + ( ["SM"] if len(eft['param'])==0 else [ "_".join([key, str(val)]) for key, val in sorted(eft['param'].iteritems())] ) )

if args.show_derivatives:
    eft_derivatives = [
        {'der':('cHW',),             'color':ROOT.kGreen+1,  'tex':'C_{HW}'},
        {'der':('cHW','cHW'),        'color':ROOT.kGreen+3,  'tex':'C^{2}_{HW}'},
        {'der':('cHWtil',),          'color':ROOT.kCyan+1,   'tex':'C_{H#tilde{W}}'},
        {'der':('cHWtil','cHWtil'),  'color':ROOT.kCyan+2,   'tex':'C^{2}_{#tilde{W}}'},
        {'der':('cHj3',),            'color':ROOT.kOrange-1, 'tex':'C_{Hq3}'},
        {'der':('cHj3','cHj3'),      'color':ROOT.kOrange-2, 'tex':'C^{2}_{Hq3}'},
        ]
else:
    eft_derivatives = []
    
for der in eft_derivatives:
    der['func'] = signal.weightInfo.get_diff_weight_func(der['der'])
    der['name'] = "_".join( ["derivative"] + list(der['der']) )

sequence = []
def make_eft_weights( event, sample):
    if sample.name!=signal.name:
        return
    SM_ref_weight         = eft_configs[0]['func'](event, sample)
    event.eft_weights     = [1] + [eft['func'](event, sample)/SM_ref_weight for eft in eft_configs[1:]]
    event.eft_derivatives = [der['func'](event, sample)/SM_ref_weight for der in eft_derivatives]


combinatoricalBTags = args.signal.startswith("ZH")
if args.signal.startswith("WH"):
    stack       = Stack( [samples.TTJets, samples.WJetsToLNu_HT] )
elif args.signal.startswith("ZH"):
    if args.DYsample == 'all':
        samples.DYJets_HT.addSelectionString("Sum$(genJet_matchBParton)==0")
        stack       = Stack( [samples.DYBBJets, samples.DYJets_HT] )
    else:
        stack       = Stack( [getattr( samples, args.DYsample )])
    plot_directory    += "_"+args.DYsample
    results_directory += "_"+args.DYsample
    if combinatoricalBTags: 
        if args.DYsample == 'all':
            samples.DYJets_HT_comb.addSelectionString("Sum$(genJet_matchBParton)==0")
            stack       = Stack( [samples.DYBBJets_comb, samples.DYJets_HT_comb] )
            stack[0][0].read_variables = [TreeVariable.fromString("combinatoricalBTagWeight2b/F")]
            stack[0][1].read_variables = [TreeVariable.fromString("combinatoricalBTagWeight2b/F")]
        else:
            stack       = Stack( [getattr( samples, stack[0][0].name+'_comb')] )
            stack[0][0].read_variables = [TreeVariable.fromString("combinatoricalBTagWeight2b/F")]

    #stack       = Stack( [samples.DYBBJets]) 
    #samples.DYBBJets.read_variables = [TreeVariable.fromString("combinatoricalBTagWeight2b/F")]

sequence.append( make_eft_weights )

eft_weights = [] 
eft_weights =  [[]]
for sample in stack.samples:
    sample.style = styles.fillStyle(sample.color)
    eft_weights[0].append( None )

for i_eft, eft in enumerate(eft_configs):
    stack.append( [signal] )
    #eft_weights.append( [get_eft_reweight(eft, signal.weightInfo)] )
    eft_weights.append( [lambda event, sample, i_eft=i_eft: event.eft_weights[i_eft]] )

for i_eft, eft in enumerate(eft_derivatives):
    stack.append( [signal] )
    #eft_weights.append( [get_eft_reweight(eft, signal.weightInfo)] )
    eft_weights.append( [lambda event, sample, i_eft=i_eft: event.eft_derivatives[i_eft]] )

if args.sign_reweight:
    #lumi_weight = lambda event, sample: lumi*event.lumiweight1fb*event.sign_reweight#*sin(2*event.VV_angles['Theta'])*sin(2*event.VV_angles['theta_V1'])
    weight_branches = ["lumiweight1fb", "sign_reweight"]
else:
    #lumi_weight = lambda event, sample: lumi*event.lumiweight1fb#*sin(2*event.VV_angles['Theta'])*sin(2*event.VV_angles['theta_V1'])
    weight_branches = ["lumiweight1fb"]

lumi  = 350
def weight_getter( branches ):
    getters = [ operator.attrgetter(branch) for branch in branches ]
    def getter( event, sample ):
#        for b, g in zip( branches, getters ):
#            print b, g(event)
#        print
        return reduce( operator.mul , [ g(event) for g in getters ], lumi ) 
    return getter

for sample in stack.samples:
    sample.weight = weight_getter( weight_branches + ( ["combinatoricalBTagWeight2b"] if "DY" in sample.name and combinatoricalBTags else [] ) )

# Read variables and sequences
jetVars          = ['pt/F', 'eta/F', 'phi/F', 'bTag/F', 'bTagPhys/I']

jetVarNames      = [x.split('/')[0] for x in jetVars]

lepVars          = ['pt/F','eta/F','phi/F','pdgId/I','isolationVar/F', 'isolationVarRhoCorr/F']
lepVarNames      = [x.split('/')[0] for x in lepVars]

# Training variables
read_variables = [\
    "nBTag/I", 
    "nBTag_loose/I",
    "recoMet_pt/F", "recoMet_phi/F",
    "genMet_pt/F", "genMet_phi/F",
    "recoZ_pt/F", "recoZ_eta/F", "recoZ_phi/F", "recoZ_mass/F", "recoZ_cosThetaStar/F", "recoZ_lldPhi/F", "recoZ_lldR/F", "recoZ_l1_index/I", "recoZ_l2_index/I",
    "nrecoJet/I",
    "recoJet[%s]"%(",".join(jetVars)),
    "nrecoLep/I",
    "recoLep[%s]"%(",".join(lepVars)),
    "ngenLep/I", "genLep[pt/F,eta/F,phi/F,pdgId/I,mother_pdgId/I]", 
    "lumiweight1fb/F",
    "genW[pt/F,eta/F,phi/F,l1_index/I,l2_index/I]", "ngenW/I",
    "evt/l", "run/I", "lumi/I",
    "H_dijet_mass/F", "H_pt/F", "H_j1_index/I", "H_j2_index/I", 
    "H_j1_index/I", "H_j2_index/I", 
    "WH_W_pt/F", "WH_dPhiMetLep/F", "WH_MT/F", "WH_nu_pt/F", "WH_nu_eta/F", "WH_nu_phi/F", "WH_nu_E/F", "WH_Theta/F", "WH_theta/F", "WH_phi/F", 
    "ZH_Theta/F", "ZH_theta/F", "ZH_phi/F",
]

preselection = [ 
    #("debug", "(evt==25857178)") 
    #("debug", "(genW_pt[0]>0)") 
]

selectionString  = "&&".join( [ c[1] for c in preselection] + ([cutInterpreter.cutString(args.selection)] if args.selection is not None else []))
subDirectory     =  '-'.join( [ c[0] for c in preselection] + ([args.selection] if args.selection is not None else []))
if subDirectory  == '': 
    subDirectory = 'inc'

for sample in stack.samples:
    if selectionString != "":
        sample.addSelectionString( selectionString )
    if args.small:
        #sample.reduceFiles( factor = 30 )
        sample.reduceFiles( to = 15 )

#BITs
import sys, os, time
sys.path.insert(0,os.path.expandvars("$CMSSW_BASE/src/BIT"))
from BoostedInformationTree import BoostedInformationTree
if signal.name.startswith('WH'):
    import TMB.BIT.configs.WH_delphes_bkgs as config
    #bits        = config.load("/groups/hephy/cms/robert.schoefbeck/BIT/models/WH_delphes/v2/")
    bits  = config.load("/groups/hephy/cms/robert.schoefbeck/BIT/models/WH_delphes_bkgs/first_try/")
elif signal.name.startswith('ZH'):
    import TMB.BIT.configs.ZH_delphes_bkgs_comb as config
    #bits        = config.load("/groups/hephy/cms/robert.schoefbeck/BIT/models/ZH_delphes/v2/")
    #bits_bkgs   = config.load("/groups/hephy/cms/robert.schoefbeck/BIT/models/ZH_delphes_bkgs/first_try/")
    bits   = config.load("/groups/hephy/cms/robert.schoefbeck/BIT/models/ZH_delphes_bkgs_comb/v2/")

for key, bit in bits.iteritems():
    bit.name = "BIT_"+"_".join( key )

bit_plots = [
    (('cHW',),             [20,-5,5]), 
    (('cHW','cHW'),        [30,-5,25]), 
    (('cHj3',),          [20,-5,5]), 
    (('cHj3','cHj3'),  [30,-5,25]), 
    (('cHW',),             [20,-1,1]), 
    (('cHW','cHW'),        [30,-5,25]), 
    (('cHj3',),            [20,-1,1]), 
    (('cHj3','cHj3'),      [30,-5,25]), 
]
sequence.extend( config.sequence )

features    = {sample.name:[] for sample in stack.samples}
predictions = {sample.name:{tuple(sorted(key)):[] for key in bits.keys()} for sample in stack.samples }

def bit_predict( event, sample ):

    for var, func in config.mva_variables:
        setattr( event, var, func(event, sample) )
    
    # get model inputs assuming lstm
    event.features = config.predict_inputs( event, sample )
    features[sample.name].append( event.features )
    for key, bit in bits.iteritems():
        #print has_lstm, flat_variables, lstm_jets
        prediction = bit.predict( event.features )
        
        setattr( event, bit.name, prediction )
        predictions[sample.name][tuple(sorted(key))].append( prediction )
#        if not prediction>-float('inf'):
#            print name, prediction, [[getattr( event, mva_variable) for mva_variable, _ in config.mva_variables]]
#            print "mva_m3", event.mva_m3, "m3", event.m3, "event.nJetGood", event.nJetGood
#            raise RuntimeError("Found NAN prediction?")

sequence.append( bit_predict )

for sample in stack.samples:
    config.add_weight_derivatives( sample )

weight_derivatives = {sample.name:{tuple(sorted(weight['comb'])):[] for weight in sample.weight_derivatives} for sample in stack.samples } 
def compute_weight_derivatives( event, sample ):
    #vector = [{weight['comb']:weight['func'](event, sample)} for weight in sample.weight_derivatives]
    #print vector
    for weight in sample.weight_derivatives:
        weight_derivatives[sample.name][tuple(sorted(weight['comb']))].append( weight['func'](event, sample) )

sequence.append( compute_weight_derivatives )

# load keras models
from keras.models import load_model

if signal.name.startswith('ZH'):
    keras_models = [
        ("ZH_TT_WJets", load_model("/groups/hephy/cms/robert.schoefbeck/TMB/models/ZH_TT_WJets/ZH_delphes_bkgs/multiclass_model.h5")),
    ]
else:
    keras_models = [
        ("WH_TT_WJets", load_model("/groups/hephy/cms/robert.schoefbeck/TMB/models/WH_TT_WJets/WH_delphes_bkgs/multiclass_model.h5")),
    ]

def keras_predict( event, sample ):

    # get model inputs assuming lstm
    for name, model in keras_models:
        prediction = model.predict( event.features.reshape(1,-1) )

        #print prediction
        for i_val, val in enumerate( prediction[0] ):
            setattr( event, name+'_'+config.training_samples[i_val].name, val)

sequence.append( keras_predict )

def make_sign_reweight( event, sample):
    event.sign_reweight = copysign( 1, sin(2*event.WH_Theta)*sin(2*event.WH_theta))
if args.sign_reweight:
    sequence.append( make_sign_reweight )

### Helpers
def addTransverseVector( p_dict ):
    ''' add a transverse vector for further calculations
    '''
    p_dict['vec2D'] = ROOT.TVector2( p_dict['pt']*cos(p_dict['phi']), p_dict['pt']*sin(p_dict['phi']) )

def addTLorentzVector( p_dict ):
    ''' add a TLorentz 4D Vector for further calculations
    '''
    p_dict['vecP4'] = ROOT.TLorentzVector( p_dict['pt']*cos(p_dict['phi']), p_dict['pt']*sin(p_dict['phi']),  p_dict['pt']*sinh(p_dict['eta']), p_dict['pt']*cosh(p_dict['eta']) )

# Use some defaults
Plot.setDefaults(stack = stack, weight = eft_weights, addOverFlowBin="upper")
 
plots        = []
plots2D      = []
postfix     = "_rw" if args.sign_reweight else "" 
postfix_tex = " weight sgn(sin(2#theta)sin(2#Theta))" if args.sign_reweight else "" 

for key_, binning in bit_plots:
    key = tuple(sorted(key))
    # 1D discriminator
    plots.append(Plot(
        name = bits[key].name+postfix,
        texX = bits[key].name+postfix_tex, texY = 'Number of Events / 10 GeV',
        attribute = lambda event, sample, var=bits[key].name: getattr(event, var),
        #binning=Binning.fromThresholds([0, 0.5, 1, 2,3,4,10]),
        binning   = binning,
        addOverFlowBin = 'upper',
    ))

for model_name, model in keras_models:
    for i_tr_s, tr_s in enumerate( config.training_samples ):
        disc_name = model_name+'_'+config.training_samples[i_tr_s].name
        plots.append(Plot(
            texX = disc_name, texY = 'Number of Events',
            name = "keras_"+disc_name, 
            attribute = lambda event, sample, disc_name=disc_name: getattr( event, disc_name ),
            binning=[50, 0, 1],
        ))

#features
for i_key, (key, _) in enumerate( config.mva_variables ):
    plots.append(Plot( name = key.replace("mva_", "")+postfix,
      texX = config.plot_options[key]['tex']+postfix_tex, texY = 'Number of Events',
      attribute = lambda event, sample, i_key=i_key: event.features[i_key],
      binning   =  config.plot_options[key]['binning'],
    ))

# Text on the plots
def drawObjects( hasData = False ):
    tex = ROOT.TLatex()
    tex.SetNDC()
    tex.SetTextSize(0.04)
    tex.SetTextAlign(11) # align right
    lines = [
      (0.15, 0.95, 'CMS Preliminary' if hasData else "CMS Simulation"), 
      #(0.45, 0.95, 'L=%3.1f fb{}^{-1} (13 TeV) Scale %3.2f'% ( lumi_scale, dataMCScale ) ) if plotData else (0.45, 0.95, 'L=%3.1f fb{}^{-1} (13 TeV)' % lumi_scale)
    ]
    return [tex.DrawLatex(*l) for l in lines] 

# draw function for plots
def drawPlots(plots, subDirectory=''):
  for log in [False, True]:
    plot_directory_ = os.path.join(plot_directory, subDirectory)
    plot_directory_ = os.path.join(plot_directory_, "log") if log else os.path.join(plot_directory_, "lin")
    for plot in plots:
        if  type(plot)==Plot2D:
            plotting.draw2D( plot,
                       plot_directory = plot_directory_,
                       logX = False, logY = False, logZ = log,
                       drawObjects = drawObjects(),
                       copyIndexPHP = True,
#                       oldColors = True,
                       ) 
        else:
            if not max(l[0].GetMaximum() for l in plot.histos): continue # Empty plot
            subtr = 0 #if args.show_derivatives else len(eft_configs)
            plotting.draw(plot,
              plot_directory = plot_directory_,
              ratio =  None,
              logX = False, logY = log, sorting = False,
              yRange = (0.03, "auto") if log else "auto",
              scaling = {},
              legend =  ( (0.17,0.9-0.05*(sum(map(len, plot.histos))-subtr)/2,1.,0.9), 2),
              drawObjects = drawObjects( ),
              copyIndexPHP = True,
            )

# Make plots & np arrays or load these
results_file = os.path.join( results_directory, "data.pkl")
if os.path.exists( results_file ) and not args.overwrite:
    features, weight_derivatives, predictions =  pickle.load( file(results_file) )
    logger.info( "Loaded data from %s" % results_file )
else:
    plotting.fill(plots+plots2D, read_variables = read_variables, sequence = sequence, max_events = -1 if args.small else -1)

    offset = 1 
    for plot in plots:
        for i_eft, eft in enumerate(eft_configs):
            plot.histos[i_eft+offset][0].legendText = eft['tex']
            plot.histos[i_eft+offset][0].style      = styles.lineStyle(eft['color'],width=2)
            plot.histos[i_eft+offset][0].SetName(eft['name'])
        for i_eft, eft in enumerate(eft_derivatives):
            if args.show_derivatives:
                plot.histos[i_eft+offset+len(eft_configs)][0].legendText = eft['tex']
                plot.histos[i_eft+offset+len(eft_configs)][0].style = styles.lineStyle(eft['color'],width=2,dashed=True)
            else:
                plot.histos[i_eft+offset+len(eft_configs)][0].legendText = None
                plot.histos[i_eft+offset+len(eft_configs)][0].style = styles.invisibleStyle()
            plot.histos[i_eft+offset+len(eft_configs)][0].SetName(eft['name'])

    #plot_phi_subtr.histos = plot_phi_subtr.histos[1:]

    drawPlots(plots+plots2D, subDirectory = subDirectory)

    features           = {key:np.array(features[key]) for key in features.keys()}
    weight_derivatives = {key1:{key2:np.array(weight_derivatives[key1][key2]) for key2 in weight_derivatives[key1].keys()} for key1 in weight_derivatives.keys()}
    predictions        = {key1:{key2:np.array(predictions[key1][key2]) for key2 in predictions[key1].keys()} for key1 in predictions.keys()}

    if not os.path.exists( results_directory ):
        os.makedirs( results_directory )
    pickle.dump( [features, weight_derivatives, predictions], file( results_file, 'w'))
    logger.info( "Dumped data to %s" % results_file )
    syncer.sync()

logger.info( "Done with prefix %s and selectionString %s", args.selection, cutInterpreter.cutString(args.selection) )

def make_q_event( order, predictions, **kwargs ):
    eft      = kwargs
    if order not in ["lin", "quad", "total"]:
        raise RuntimeError("Order %s not known" % order )
    #index = {}
    #for i_comb, comb in enumerate(config.bit_derivatives):
    #    index[comb] = i_comb
    #    index[tuple(reversed(comb))] = i_comb
    result = np.zeros(len(predictions.values()[0]))
    if order in ["lin", "total"]:
        for coeff in eft.keys():
            result += eft[coeff]*predictions[(coeff,)]
    if order in ["quad", "total"]:
        for coeff1 in eft.keys():
            for coeff2 in eft.keys():
                result += .5*eft[coeff1]*eft[coeff2]*predictions[tuple(sorted([coeff1, coeff2]))]
    return result

def make_weights( weight_derivatives, **kwargs ):
    eft      = kwargs

    result = copy.deepcopy(weight_derivatives[()])
    for coeff1 in eft.keys():
        result += eft[coeff1]*weight_derivatives[(coeff1,)]
        for coeff2 in eft.keys():
            result += .5*eft[coeff1]*eft[coeff2]*weight_derivatives[tuple(sorted([coeff1, coeff2]))]
    return result

def make_cdf_map( x, y ):
    import scipy.interpolate
    map__ = scipy.interpolate.interp1d(x, y, kind='linear')
    max_x, min_x = max(x), min(x)
    max_y, min_y = max(y), min(y)
    def map_( x_ ):
        x__ = np.array(x_)
        result = np.zeros_like(x__).astype('float')
        result[x__>max_x] = max_y
        result[x__<min_x] = min_y
        vals = (x__>=min_x) & (x__<=max_x)
        result[vals] = map__(x__[vals])
        return result

    return map_

def getContours( h, level):
    _h     = h.Clone()
    ctmp = ROOT.TCanvas()
    _h.SetContour(1,array.array('d', [level]))
    _h.Draw("contzlist")
    _h.GetZaxis().SetRangeUser(0.0001,1)
    ctmp.Update()
    contours = ROOT.gROOT.GetListOfSpecials().FindObject("contours")
    return contours.At(0).Clone()


n_toys = 50000

# do not make the following inconsistent
levels          = [ 0.95, 0.68]
quantile_levels = [ 0.025, 0.16, .5, 1-0.16, 1-0.025 ]

exp_nll_ratio = {}
# interpret args
WC1, WC1_min, WC1_max, WC2, WC2_min, WC2_max = map( lambda x:float(x[1]) if x[0] in [1,2,4,5] else x[1], enumerate(args.WCs) )

step1 = (WC1_max-WC1_min)/args.nBins
step2 = (WC2_max-WC2_min)/args.nBins
WC1_vals = np.arange(WC1_min, WC1_max+step1, (WC1_max-WC1_min)/args.nBins)
WC2_vals = np.arange(WC2_min, WC2_max+step2, (WC2_max-WC2_min)/args.nBins)

test_statistics = ["total", "lin", "quad"]
exp_nll_ratio = {}
power         = {}
for test_statistic in test_statistics:
    print "Test statistic", test_statistic
    power[test_statistic] = {level:ROOT.TH2D("power_"+test_statistic, "power_"+test_statistic, len(WC1_vals)-1, array.array('d', WC1_vals), len(WC2_vals)-1, array.array('d', WC2_vals)) for level in levels}
    exp_nll_ratio[test_statistic] = ROOT.TH2D("exp_nll_ratio_"+test_statistic, "exp_nll_ratio_"+test_statistic, len(WC1_vals)-1, array.array('d', WC1_vals), len(WC2_vals)-1, array.array('d', WC2_vals))

    for i_WC1_val, WC1_val in enumerate( WC1_vals ):
        for i_WC2_val, WC2_val in enumerate( WC2_vals ):
            if WC1_val==WC2_val==0: continue

            eft = {WC1:WC1_val, WC2:WC2_val}

            for sample in [signal]: #stack.samples:#"bkg"]:

                q_event = make_q_event( test_statistic, predictions[sample.name], **eft )

                # compute BSM weights
                w_sm   = lumi/float(config.scale_weight)*make_weights( weight_derivatives[sample.name] )
                w_bsm  = lumi/float(config.scale_weight)*make_weights( weight_derivatives[sample.name], **eft)

                q_event_argsort     = np.argsort(q_event)
                q_event_argsort_inv = np.argsort(q_event_argsort)
                cdf_sm = np.cumsum(w_sm[q_event_argsort])
                cdf_sm/=cdf_sm[-1]

                # map to the SM CDF of q
                if sample.name==signal.name: 
                    cdf_map = make_cdf_map( q_event[q_event_argsort], cdf_sm )

                #q_event_cdf = cdf_sm[q_event_argsort_inv] #uniformly distributed under the SM hypothesis
                q_event_cdf = cdf_map( q_event )

                #min_, max_ = 0, 1 
                binning = np.linspace(0, 1, args.nBinsTestStat+1)

                np_histo_sm  = np.histogram(q_event_cdf, bins=binning, weights = w_sm )[0]
                np_histo_bsm = np.histogram(q_event_cdf, bins=binning, weights = w_bsm )[0]

                # Expectation_BSM( -Log( Prod_i( Pois_i( n_i, lambda_i(theta))/Pois_i( n_i, lambda_i(0)) ) ))
                exp_nll_ratio_ = 2*np.sum(np_histo_sm - np_histo_bsm - np_histo_bsm*np.log(np_histo_sm/np_histo_bsm))
                exp_nll_ratio[test_statistic].SetBinContent( exp_nll_ratio[test_statistic].FindBin( WC1_val, WC2_val ), exp_nll_ratio_)

                binned_toys_sm    = np.random.poisson(lam=np_histo_sm, size=(n_toys, len(np_histo_sm)))
                binned_toys_theta = np.random.poisson(lam=np_histo_bsm, size=(n_toys, len(np_histo_bsm)))

                q_theta_given_sm    = [ np.sum( toy_ll ) for toy_ll in 2*(np_histo_sm - np_histo_bsm  - binned_toys_sm*np.log(np_histo_sm/np_histo_bsm))]
                q_theta_given_theta = [ np.sum( toy_ll ) for toy_ll in 2*(np_histo_sm - np_histo_bsm  - binned_toys_theta*np.log(np_histo_sm/np_histo_bsm))]

                if True:
                    n = float(len(q_theta_given_sm))
                    mean_q_theta_given_sm     = np.sum(q_theta_given_sm)/n
                    sigma_q_theta_given_sm    = sqrt( np.sum((q_theta_given_sm - mean_q_theta_given_sm)**2)/(n-1) )
                    q_theta_given_sm    = (q_theta_given_sm - mean_q_theta_given_sm)/sigma_q_theta_given_sm
                    q_theta_given_theta = (q_theta_given_theta - mean_q_theta_given_sm)/sigma_q_theta_given_sm

                print i_WC1_val, WC1_val,  i_WC2_val, WC2_val,  "sqrt(2NLL)", sqrt(abs(exp_nll_ratio_))

               # Exclusion: The null hypothesis is the BSM point, the alternate is the SM.
                quantiles_theta = np.quantile( q_theta_given_theta, quantile_levels )
                quantiles_sm    = np.quantile( q_theta_given_sm, quantile_levels )
                size_           = np.sum(np.histogram( q_theta_given_sm, quantiles_sm)[0])/float(n_toys)
                #power_histo     = np.histogram( q_theta_given_theta, quantiles_sm)
                for i_level, level in enumerate(levels):
                    #if level != 0.68: continue
                    power_toy_count = np.count_nonzero((q_theta_given_sm>=quantiles_theta[i_level]) & (q_theta_given_sm<quantiles_theta[-1-i_level]))
                    power_ = 1. - power_toy_count/float(n_toys)
                    power[test_statistic][level].SetBinContent( power[test_statistic][level].FindBin( WC1_val, WC2_val ), power_ )
                #power_ = 1-np.sum(np.histogram( q_theta_given_theta, quantiles_sm)[0])/float(n_toys)
                    print "theta", round(WC1_val,3), round(WC2_val,3), "level", level, "size", quantile_levels[-1-i_level] - quantile_levels[i_level], "power", round(power_,3), test_statistic, WC1, WC2

colors   = { 'quad':ROOT.kRed, 'lin':ROOT.kBlue, 'total':ROOT.kBlack}
nll_levels = [2.27, 5.99]
contours = { key:{level:getContours( exp_nll_ratio[key], level) for level in nll_levels} for key in exp_nll_ratio.keys() }
contour_objects = []
for test_statistic in contours.keys():
    for level in contours[test_statistic].keys():
        for i_tgraph in range(contours[test_statistic][level].GetSize()):
            tgr = contours[test_statistic][level].At(i_tgraph)
            print i_tgraph, tgr, test_statistic, level

            tgr.SetLineColor(colors[test_statistic])
            tgr.SetLineWidth(2)
            tgr.SetLineStyle(ROOT.kDashed if level!=nll_levels[-1] else 1)
            tgr.SetMarkerStyle(0)
            contour_objects.append( tgr )

for test_statistic in test_statistics:
    plot2D = Plot2D.fromHisto(name = "exp_nll_ratio_%s_%s_vs_%s_lumi_%3.2f_nBinsTestStat_%i"%(test_statistic, WC1, WC2, lumi, args.nBinsTestStat), histos = [[exp_nll_ratio[test_statistic]]], texX = tex[WC1], texY = tex[WC2] )
    plotting.draw2D(plot2D, 
        plot_directory = os.path.join( plot_directory, subDirectory, "log"), 
        logY = False, logX = False, logZ = True, 
        copyIndexPHP=True, 
        drawObjects = contour_objects, 
        zRange = (0.05,25),
        histModifications = [lambda h:ROOT.gStyle.SetPalette(58)],
        )

colors   = { 'quad':ROOT.kRed, 'lin':ROOT.kBlue, 'total':ROOT.kBlack}

contours = { key:{level:getContours( power[key][level], 0.5 ) for level in levels } for key in power.keys()}
contour_objects = []
for test_statistic in contours.keys():
    for level in contours[test_statistic].keys():
        for i_tgraph in range(contours[test_statistic][level].GetSize()):
            tgr = contours[test_statistic][level].At(i_tgraph)
            print i_tgraph, tgr, test_statistic, level

            tgr.SetLineColor(colors[test_statistic])
            tgr.SetLineWidth(2)
            tgr.SetLineStyle(ROOT.kDashed if level!=0.95 else 1)
            tgr.SetMarkerStyle(0)
            contour_objects.append( tgr )

for test_statistic in test_statistics:
    for level in levels:
        plot2D = Plot2D.fromHisto(name = "power_%s_%s_vs_%s_lumi_%3.2f_level_%3.2f_nBinsTestStat_%i"%(test_statistic, WC1, WC2, lumi, level, args.nBinsTestStat), histos = [[power[test_statistic][level]]], texX = tex[WC1], texY = tex[WC2] )
        plotting.draw2D(plot2D, 
            plot_directory = os.path.join( plot_directory, subDirectory, "log"), 
            logY = False, logX = False, logZ = True, 
            copyIndexPHP=True, 
            drawObjects = contour_objects, 
            zRange = (0.05,1),
            histModifications = [lambda h:ROOT.gStyle.SetPalette(58)],
            )

syncer.sync()
