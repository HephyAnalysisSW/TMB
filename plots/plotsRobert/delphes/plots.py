#!/usr/bin/env python
''' Analysis script for standard plots
'''

# Standard imports and batch mode
import ROOT, os, itertools
#ROOT.gROOT.SetBatch(True)
import copy
import operator
import random
from math                           import sqrt, cos, sin, pi, isnan, sinh, cosh, log, copysign

# Analysis
import Analysis.Tools.syncer        as syncer
from   Analysis.Tools.WeightInfo    import WeightInfo
from   Analysis.Tools.helpers       import deltaPhi, deltaR, getObjDict

# RootTools
from RootTools.core.standard        import *

# TMB
from TMB.Tools.user                 import plot_directory
from TMB.Tools.helpers              import deltaPhi, getCollection, deltaR, mZ
from TMB.Tools.delphesCutInterpreter import cutInterpreter
import TMB.Tools.VV_angles          as VV_angles
from TMB.Tools.genObjectSelection   import isGoodGenJet

# Arguments
import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--logLevel',           action='store',      default='INFO',          nargs='?', choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'], help="Log level for logging")
argParser.add_argument('--plot_directory',     action='store',      default='delphes')
argParser.add_argument('--selection',          action='store',      default='singlelep-WHJet-onH')
argParser.add_argument('--signal',             action='store',      default='WH', choices = ['WH', 'ZH'])
argParser.add_argument('--small',                                   action='store_true',     help='Run only on a small subset of the data?')
argParser.add_argument('--no_bkgs',                                 action='store_true',     help='Do not show bkg plots?')
argParser.add_argument('--show_derivatives',                        action='store_true',     help='Show also the derivatives?')
argParser.add_argument('--sign_reweight',                           action='store_true',     help='Apply sign(sin(2theta)sin(2Theta)) as weight ')
args = argParser.parse_args()

# Logger'singlelep-WHJet' if sample.name=='WH' else 'dilep-ZHJet-onZ'
import TMB.Tools.logger as logger
import RootTools.core.logger as logger_rt
logger    = logger.get_logger(   args.logLevel, logFile = None)
logger_rt = logger_rt.get_logger(args.logLevel, logFile = None)

plot_directory = os.path.join(plot_directory, args.plot_directory,  args.signal )
if args.small: plot_directory += "_small"

# Import samples
import TMB.Samples.pp_gen_v10 as samples
    
signal = getattr( samples, args.signal)
 
# WeightInfo
signal.weightInfo = WeightInfo(signal.reweight_pkl)
signal.weightInfo.set_order(2)
signal.read_variables = [VectorTreeVariable.fromString( "p[C/F]", nMax=200 )]

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
        {'der':('cHW',),             'color':ROOT.kGreen+1,  'tex':'c_{HW}'},
        {'der':('cHW','cHW'),        'color':ROOT.kGreen+3,  'tex':'c^{2}_{HW}'},
        {'der':('cHWtil',),          'color':ROOT.kCyan+1,   'tex':'c_{H#tilde{W}}'},
        {'der':('cHWtil','cHWtil'),  'color':ROOT.kCyan+2,   'tex':'c^{2}_{#tilde{W}}'},
        {'der':('cHj3',),            'color':ROOT.kOrange-1, 'tex':'c_{Hq3}'},
        {'der':('cHj3','cHj3'),      'color':ROOT.kOrange-2, 'tex':'c^{2}_{Hq3}'},
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

if args.no_bkgs:
    stack = Stack( )
else: 
    if args.signal == "WH":
        stack       = Stack( [samples.TTJets, samples.WJetsToLNu_HT] )
    elif args.signal == "ZH":
        stack       = Stack( [samples.DYBBJets]) 

sequence.append( make_eft_weights )

eft_weights = [] 
if not args.no_bkgs:
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

lumi  = 59.7
if args.sign_reweight:
    lumi_weight = lambda event, sample: lumi*event.lumiweight1fb*event.sign_reweight#*sin(2*event.VV_angles['Theta'])*sin(2*event.VV_angles['theta_V1'])
else:
    lumi_weight = lambda event, sample: lumi*event.lumiweight1fb#*sin(2*event.VV_angles['Theta'])*sin(2*event.VV_angles['theta_V1'])

for sample in stack.samples:
    sample.weight = lumi_weight

# Read variables and sequences
jetVars          = ['pt/F', 'eta/F', 'phi/F', 'bTag/F', 'bTagPhys/I']

jetVarNames      = [x.split('/')[0] for x in jetVars]

lepVars          = ['pt/F','eta/F','phi/F','pdgId/I','isolationVar/F', 'isolationVarRhoCorr/F']
lepVarNames      = [x.split('/')[0] for x in lepVars]

# Training variables
read_variables = [\
#    "nBTag/I", 
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
if signal.name == 'WH':
    import TMB.BIT.configs.WH_delphes as config
    bits        = config.load("/groups/hephy/cms/robert.schoefbeck/BIT/models/WH_delphes/first_try/")
    bits_bkgs   = config.load("/groups/hephy/cms/robert.schoefbeck/BIT/models/WH_delphes_bkgs/first_try/")
elif signal.name == 'ZH':
    import TMB.BIT.configs.ZH_delphes as config
    bits        = config.load("/groups/hephy/cms/robert.schoefbeck/BIT/models/ZH_delphes/first_try/")
    bits_bkgs   = config.load("/groups/hephy/cms/robert.schoefbeck/BIT/models/ZH_delphes_bkgs/first_try/")

models = [
    ("BIT_cHW",             bits[('cHW',)],             [20,-5,5]), 
    ("BIT_cHW_cHW",         bits[('cHW','cHW')],        [30,-5,25]), 
#    ("BIT_cHWtil",          bits[('cHWtil',)],          [20,-5,5]), 
#    ("BIT_cHWtil_cHWtil",   bits[('cHWtil','cHWtil')],  [30,-5,25]), 
    ("BIT_cHj3",            bits[('cHj3',)],          [20,-5,5]), 
    ("BIT_cHj3_cHj3",       bits[('cHj3','cHj3')],  [30,-5,25]), 
    ("BIT_bkgs_cHW",             bits_bkgs[('cHW',)],             [20,-1,1]), 
    ("BIT_bkgs_cHW_cHW",         bits_bkgs[('cHW','cHW')],        [30,-5,25]), 
#    ("BIT_bkgs_cHWtil",          bits_bkgs[('cHWtil',)],          [20,-1,1]), 
#    ("BIT_bkgs_cHWtil_cHWtil",   bits_bkgs[('cHWtil','cHWtil')],  [30,-5,25]), 
    ("BIT_bkgs_cHj3",            bits_bkgs[('cHj3',)],            [20,-1,1]), 
    ("BIT_bkgs_cHj3_cHj3",       bits_bkgs[('cHj3','cHj3')],      [30,-5,25]), 
]
sequence.extend( config.sequence )

def bit_predict( event, sample ):

    for var, func in config.mva_variables:
        setattr( event, var, func(event, sample) )
    
    # get model inputs assuming lstm
    event.features = config.predict_inputs( event, sample)
    for name, model, _ in models:
        #print has_lstm, flat_variables, lstm_jets
        prediction = model.predict( event.features )
        setattr( event, name, prediction )
#        if not prediction>-float('inf'):
#            print name, prediction, [[getattr( event, mva_variable) for mva_variable, _ in config.mva_variables]]
#            print "mva_m3", event.mva_m3, "m3", event.m3, "event.nJetGood", event.nJetGood
#            raise RuntimeError("Found NAN prediction?")

sequence.append( bit_predict )

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

for model_name, _, binning in models:

    # 1D discriminator
    plots.append(Plot(
        name = model_name+postfix,
        texX = model_name+postfix_tex, texY = 'Number of Events / 10 GeV',
        attribute = lambda event, sample, model_name=model_name: getattr(event, model_name),
        #binning=Binning.fromThresholds([0, 0.5, 1, 2,3,4,10]),
        binning   = binning,
        addOverFlowBin = 'upper',
    ))

#    # 2D with background 
#    var = "mva_H_pt"
#    i_key = [v[0] for v in config.mva_variables].index(var)
#    plots2D.append(Plot2D(
#        name = "bkg2D_H_pt_"+model_name+postfix,
#        texX = model_name, texY = config.plot_options[var]['tex'],
#        stack = Stack(*stack[0:1]),
#        attribute = (
#            lambda event, sample, model_name=model_name: getattr(event, model_name),
#            lambda event, sample, i_key=i_key: event.features[i_key],
#            ),
#        #binning=Binning.fromThresholds([0, 0.5, 1, 2,3,4,10]),
#        binning   = binning+config.plot_options[var]['binning'],
#    ))
#
#    # 2D with signal
#    plots2D.append(Plot2D(
#        name = "sig2D_H_pt_"+model_name+postfix,
#        texX = model_name, texY = config.plot_options[var]['tex'],
#        stack = Stack(*stack[1:2]),
#        attribute = (
#            lambda event, sample, model_name=model_name: getattr(event, model_name),
#            lambda event, sample, i_key=i_key: event.features[i_key],
#            ),
#        #binning=Binning.fromThresholds([0, 0.5, 1, 2,3,4,10]),
#        binning   = binning+config.plot_options[var]['binning'],
#    ))

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

plotting.fill(plots+plots2D, read_variables = read_variables, sequence = sequence, max_events = -1 if args.small else -1)

#plot_phi_subtr = copy.deepcopy(plots[subtr_phi_ind])
#plot_phi_subtr.name = plot_phi_subtr.name.replace("pos_","subtr_")
#plot_phi_subtr.texX = plot_phi_subtr.texX.replace("(pos.)","(subtr)")
#for i_h_, h_ in enumerate(plots[subtr_phi_ind].histos):
#    for i_h, h in enumerate(h_):
#        h_sub = plots[subtr_phi_ind+1].histos[i_h_][i_h].Clone()
#        h_sub.Scale(-1)
#        plot_phi_subtr.histos[i_h_][i_h].Add(h_sub)
#
#plots.append( plot_phi_subtr )

#color EFT
offset = 0 if args.no_bkgs else 1
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

logger.info( "Done with prefix %s and selectionString %s", args.selection, cutInterpreter.cutString(args.selection) )

syncer.sync()

