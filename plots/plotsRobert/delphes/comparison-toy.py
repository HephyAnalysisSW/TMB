#!/usr/bin/env python
''' Comparison of toy and delphes 
'''

# Standard imports and batch mode
import ROOT, sys, os, itertools
#ROOT.gROOT.SetBatch(True)
import copy
import operator
import random
from math                           import sqrt, cos, sin, pi, isnan, sinh, cosh, log, copysign
import numpy as np
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
import TMB.Tools.stat_helpers       as stat_helpers

# Arguments
import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--logLevel',           action='store',      default='INFO',          nargs='?', choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'], help="Log level for logging")
argParser.add_argument('--plot_directory',     action='store',      default='BIT_VH_9')
argParser.add_argument('--signal',             action='store',      default='ZH', choices = ['WH', 'ZH'])
argParser.add_argument('--small',                                   action='store_true',     help='Run only on a small subset of the data?')
args = argParser.parse_args()

# Logger
import TMB.Tools.logger as logger
import RootTools.core.logger as logger_rt
logger    = logger.get_logger(   args.logLevel, logFile = None)
logger_rt = logger_rt.get_logger(args.logLevel, logFile = None)

plot_directory = os.path.join(plot_directory, args.plot_directory,  "comparison-toy", args.signal )
if args.small: plot_directory += "_small"

# Import samples
import TMB.Samples.pp_gen_v10 as samples
    
signal = getattr( samples, args.signal)
if args.signal == 'WH':
    selection = "singlelep-WHJet-onH-realW-ptW300"
elif args.signal == 'ZH':
    selection = "dilep-ZHJet-onZ-onH-ptZ300"

# WeightInfo
signal.weightInfo = WeightInfo(signal.reweight_pkl)
signal.weightInfo.set_order(2)
signal.read_variables = [VectorTreeVariable.fromString( "p[C/F]", nMax=200 )]

eft_configs = [
    {'color':ROOT.kBlack,       'delphes':{}, 'tex':"SM"},
#    {'color':ROOT.kMagenta-4,   'delphes':{'cHWtil':-1},  'tex':"c_{H#tilde{W}}=-1"},
#    {'color':ROOT.kMagenta+2,   'delphes':{'cHWtil':1},   'tex':"c_{H#tilde{W}}=1"},
#    {'color':ROOT.kGreen-4,     'delphes':{'cHW':-1},     'tex':"c_{HW}=-1"},
#    {'color':ROOT.kGreen+2,     'delphes':{'cHW':1},      'tex':"c_{HW}=1"},
#    {'color':ROOT.kBlue-4,      'delphes':{'cHj3':-.1},   'tex':"c_{HQ}^{(3)}=-0.1"},
#    {'color':ROOT.kBlue+2,      'delphes':{'cHj3':.1},    'tex':"c_{HQ}^{(3)}=0.1"},
#    {'color':ROOT.kMagenta-4,  'delphes':{'cHW':-2.5},  'tex':"c_{HW}=-2.5"},
#    {'color':ROOT.kMagenta+2,  'delphes':{'cHW':-2},   'tex':"c_{HW}=-2"},
#    {'color':ROOT.kGreen-4,    'delphes':{'cHW':-1.5},     'tex':"c_{HW}=-1.5"},
#    {'color':ROOT.kGreen+2,    'delphes':{'cHW':-1},      'tex':"c_{HW}=-1"},
#    {'color':ROOT.kBlue-4,     'delphes':{'cHW':-.5},   'tex':"c_{HW}=-.5"},
#    {'color':ROOT.kBlue+2,     'delphes':{'cHW':.5},    'tex':"c_{HW}=.5"},
#    {'color':ROOT.kOrange-4,   'delphes':{'cHW':1.},   'tex':"c_{HW}=1"},
#    {'color':ROOT.kOrange+2,   'delphes':{'cHW':1.5},    'tex':"c_{HW}=1.5"},
#    {'color':ROOT.kRed-4,      'delphes':{'cHW':2},   'tex':"c_{HW}=2"},
#    {'color':ROOT.kRed+2,      'delphes':{'cHW':2.5},    'tex':"c_{HW}=2.5"},
    {'color':ROOT.kMagenta-4,  'delphes':{'cHWtil':-2.5},  'tex':"c_{H#tilde{W}}=-2.5"},
    {'color':ROOT.kMagenta+2,  'delphes':{'cHWtil':-2},   'tex':"c_{H#tilde{W}}=-2"},
    {'color':ROOT.kGreen-4,    'delphes':{'cHWtil':-1.5},     'tex':"c_{H#tilde{W}}=-1.5"},
    {'color':ROOT.kGreen+2,    'delphes':{'cHWtil':-1},      'tex':"c_{H#tilde{W}}=-1"},
    {'color':ROOT.kBlue-4,     'delphes':{'cHWtil':-.5},   'tex':"c_{H#tilde{W}}=-.5"},
    {'color':ROOT.kBlue+2,     'delphes':{'cHWtil':.5},    'tex':"c_{H#tilde{W}}=.5"},
    {'color':ROOT.kOrange-4,   'delphes':{'cHWtil':1.},   'tex':"c_{H#tilde{W}}=1"},
    {'color':ROOT.kOrange+2,   'delphes':{'cHWtil':1.5},    'tex':"c_{H#tilde{W}}=1.5"},
    {'color':ROOT.kRed-4,      'delphes':{'cHWtil':2},   'tex':"c_{H#tilde{W}}=2"},
    {'color':ROOT.kRed+2,      'delphes':{'cHWtil':2.5},    'tex':"c_{H#tilde{W}}=2.5"},
#    {'color':ROOT.kMagenta-4,  'delphes':{'cHj3':-.25},  'tex':"c_{HQ}^{(3)}=-.25"},
#    {'color':ROOT.kMagenta+2,  'delphes':{'cHj3':-.2},   'tex':"c_{HQ}^{(3)}=-.2"},
#    {'color':ROOT.kGreen-4,    'delphes':{'cHj3':-.15},  'tex':"c_{HQ}^{(3)}=-.15"},
#    {'color':ROOT.kGreen+2,    'delphes':{'cHj3':-.1},   'tex':"c_{HQ}^{(3)}=-.1"},
#    {'color':ROOT.kBlue-4,     'delphes':{'cHj3':-.05},  'tex':"c_{HQ}^{(3)}=-.5"},
#    {'color':ROOT.kBlue+2,     'delphes':{'cHj3':.05},   'tex':"c_{HQ}^{(3)}=.5"},
#    {'color':ROOT.kOrange-4,   'delphes':{'cHj3':.1},    'tex':"c_{HQ}^{(3)}=.1"},
#    {'color':ROOT.kOrange+2,   'delphes':{'cHj3':.15},   'tex':"c_{HQ}^{(3)}=.15"},
#    {'color':ROOT.kRed-4,      'delphes':{'cHj3':.2},    'tex':"c_{HQ}^{(3)}=.2"},
#    {'color':ROOT.kRed+2,      'delphes':{'cHj3':.25},   'tex':"c_{HQ}^{(3)}=.25"},
    ]

# import the toy VH model
sys.path.insert(0,os.path.expandvars("$CMSSW_BASE/src/BIT"))
import VH_models
toy_model = getattr(VH_models, "WH_Spannowsky" if args.signal=="WH" else "ZH_Spannowsky")

for eft in eft_configs:
    eft['func'] = signal.weightInfo.get_weight_func(**eft['delphes']) 
    eft['name'] = "_".join( ["signal"] + ( ["SM"] if len(eft['delphes'])==0 else [ "_".join([key, str(val)]) for key, val in sorted(eft['delphes'].iteritems())] ) )
    eft['toy']  = copy.deepcopy(eft['delphes'])
    if eft['toy'].has_key( 'cHj3' ):
        eft['toy']['cHQ3'] = eft['toy']['cHj3'] 
        del eft['toy']['cHj3']
    eft['toy'] = toy_model.make_eft( **eft['toy'] )
    

toy_feature_names = toy_model.feature_names

nEvents  = 100000
toy_features = toy_model.getEvents(nEvents)
nEvents  = len(toy_features)
for eft in eft_configs:
    eft['toy_weights'] = toy_model.getWeights(toy_features, eft=eft['toy'])
    logger.info( "Produced toy weights for %r"%eft['toy'] )
logger.info("Created toy data set of size %i" % nEvents )

sequence = []
def make_eft_weights( event, sample):
    if sample.name!=signal.name:
        return
    SM_ref_weight         = eft_configs[0]['func'](event, sample)
    event.eft_weights     = [1] + [eft['func'](event, sample)/SM_ref_weight for eft in eft_configs[1:]]
sequence.append( make_eft_weights )

stack       = Stack( )
eft_weights = [] 

for i_eft, eft in enumerate(eft_configs):
    stack.append( [signal] )
    eft_weights.append( [lambda event, sample, i_eft=i_eft: event.eft_weights[i_eft]] )
    #eft_weights.append( [get_eft_reweight(eft, signal.weightInfo)] )

lumi  = 137
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
]

selectionString  = "&&".join( [ c[1] for c in preselection] + [cutInterpreter.cutString(selection)] )
subDirectory     =  '-'.join( [ c[0] for c in preselection] + [selection])

for sample in stack.samples:
    if selectionString != "":
        sample.addSelectionString( selectionString )
    if args.small:
        #sample.reduceFiles( factor = 30 )
        sample.reduceFiles( to = 15 )

#BITs
import sys, os, time
from BoostedInformationTree import BoostedInformationTree
if signal.name == 'WH':
    import TMB.BIT.configs.WH_delphes as config
    bits        = config.load("/groups/hephy/cms/robert.schoefbeck/BIT/models/WH_delphes/v2/")
    bits_bkgs   = config.load("/groups/hephy/cms/robert.schoefbeck/BIT/models/WH_delphes_bkgs/first_try/")
elif signal.name == 'ZH':
    import TMB.BIT.configs.ZH_delphes as config
    bits        = config.load("/groups/hephy/cms/robert.schoefbeck/BIT/models/ZH_delphes/v2/")
    bits_bkgs   = config.load("/groups/hephy/cms/robert.schoefbeck/BIT/models/ZH_delphes_bkgs/first_try/")

models = [
#    ("BIT_cHW",             bits[('cHW',)],             [20,-5,5]), 
#    ("BIT_cHW_cHW",         bits[('cHW','cHW')],        [30,-5,25]), 
#    ("BIT_cHj3",            bits[('cHj3',)],          [20,-5,5]), 
#    ("BIT_cHj3_cHj3",       bits[('cHj3','cHj3')],  [30,-5,25]), 
#    ("BIT_bkgs_cHW",             bits_bkgs[('cHW',)],             [20,-1,1]), 
#    ("BIT_bkgs_cHW_cHW",         bits_bkgs[('cHW','cHW')],        [30,-5,25]), 
#    ("BIT_bkgs_cHj3",            bits_bkgs[('cHj3',)],            [20,-1,1]), 
#    ("BIT_bkgs_cHj3_cHj3",       bits_bkgs[('cHj3','cHj3')],      [30,-5,25]), 
]
sequence.extend( config.sequence )

def bit_predict( event, sample ):

    for var, func in config.mva_variables:
        setattr( event, var, func(event, sample) )
    
    # get model inputs assuming lstm
    event.features = config.predict_inputs( event, sample)
#    for name, model, _ in models:
#        #print has_lstm, flat_variables, lstm_jets
#        prediction = model.predict( event.features )
#        setattr( event, name, prediction )

sequence.append( bit_predict )

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

#for model_name, _, binning in models:
#
#    # 1D discriminator
#    plots.append(Plot(
#        name = model_name,
#        texX = model_name, texY = 'Number of Events / 10 GeV',
#        attribute = lambda event, sample, model_name=model_name: getattr(event, model_name),
#        #binning=Binning.fromThresholds([0, 0.5, 1, 2,3,4,10]),
#        binning   = binning,
#        addOverFlowBin = 'upper',
#    ))

correspondence = {
    ('mva_Z_pt' if args.signal=="ZH" else 'mva_W_pt'):{'toy_name':'pT', 'features':toy_features[:, toy_model.feature_names.index('pT')]},
}

#features
for i_key, (key, _) in enumerate( config.mva_variables ):
    if not key in correspondence.keys(): continue

    plots.append(Plot( name = key.replace("mva_", ""),
      texX = config.plot_options[key]['tex'], texY = 'Number of Events',
      attribute = lambda event, sample, i_key=i_key: event.features[i_key],
      binning   =  config.plot_options[key]['binning'],
    ))

    binning = config.plot_options[key]['binning']
    plots[-1].toy_histos = [ stat_helpers.make_TH1F( np.histogram( 
            correspondence[key]['features'], 
            bins    = np.linspace(binning[1], binning[2], binning[0]), 
            weights = eft_config['toy_weights'][()] )) for eft_config in eft_configs ]

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

## draw function for plots
#def drawPlots(plots, subDirectory=''):
#  for log in [False, True]:
#    plot_directory_ = os.path.join(plot_directory, subDirectory)
#    plot_directory_ = os.path.join(plot_directory_, "log") if log else os.path.join(plot_directory_, "lin")
#    for plot in plots:
#            if not max(l[0].GetMaximum() for l in plot.histos): continue # Empty plot
#            subtr = 0 
#            plotting.draw(plot,
#              plot_directory = plot_directory_,
#              ratio =  None,
#              logX = False, logY = log, sorting = False,
#              yRange = (0.03, "auto") if log else "auto",
#              scaling = {},
#              legend =  ( (0.17,0.9-0.05*(sum(map(len, plot.histos))-subtr)/2,1.,0.9), 2),
#              drawObjects = drawObjects( ),
#              copyIndexPHP = True,
#            )

plotting.fill(plots, read_variables = read_variables, sequence = sequence, max_events = -1 if args.small else -1)

#color EFT
#draw_plots = []
for plot in plots:
    scale = plot.histos[0][0].Integral() / plot.toy_histos[0].Integral()
    logger.info( "Scaling %s. Delphes: %3.2f for a lumi of %3.2f. Toy normalistion: %3.2f. Scale factor: %3.2f", plot.name, plot.histos[0][0].Integral(), lumi, plot.toy_histos[0].Integral(), scale )

    tot_xsec_toy     = ROOT.TGraph(len(eft_configs))
    tot_xsec_delphes = ROOT.TGraph(len(eft_configs))

    for i_eft, eft in enumerate(eft_configs):
        plot.histos[i_eft][0].legendText = eft['tex']
        plot.histos[i_eft][0].style      = styles.lineStyle(eft['color'],width=2)
        plot.histos[i_eft][0].SetName(eft['name'])

        # scale according to SM
        plot.toy_histos[i_eft].Scale( scale )
        plot.toy_histos[i_eft].legendText = eft['tex'] + ' (toy)'
        plot.toy_histos[i_eft].style      = styles.lineStyle(eft['color'],width=2, dashed=True)
        plot.toy_histos[i_eft].SetName(eft['name'])

        param_val    = 0 if i_eft==0 else eft['delphes'].values()[0]
        xsec_toy     = plot.toy_histos[i_eft].Integral()
        xsec_delphes = plot.histos[i_eft][0].Integral()

        if len(eft['delphes'])==1:
            WC = eft['delphes'].keys()[0]
            tot_xsec_toy.GetXaxis().SetTitle( WC )

        print eft['delphes'], xsec_toy, xsec_delphes
        tot_xsec_toy    .SetPoint(i_eft, param_val, xsec_toy ) 
        tot_xsec_delphes.SetPoint(i_eft, param_val, xsec_delphes) 

    draw_plot = Plot.fromHisto( plot.name, plot.histos+[[h] for h in plot.toy_histos], texX = plot.texX, texY = plot.texY )
    draw_plot.stack = None
    for log in [False, True]:
        plot_directory_ = os.path.join(plot_directory, "log") if log else os.path.join(plot_directory, "lin")
        plotting.draw( 
            plot = draw_plot,
            plot_directory = plot_directory_,
            #ratio          = {'yRange':(0.6,1.4)} if len(plot.stack)>=2 else None,
            logY = log, logX = False, sorting = False,
            legend         = ( (0.15,0.7,0.9,0.92),2),
            drawObjects    = [], #quantile_lines + drawObjects,
            copyIndexPHP   = True,
            extensions     = ["png"],
          )

    c1 = ROOT.TCanvas()
    tot_xsec_toy.SetMarkerColor(ROOT.kGreen)
    tot_xsec_toy.Draw("AP")
    tot_xsec_delphes.SetMarkerColor(ROOT.kRed)
    tot_xsec_delphes.Draw("Psame")
    l = ROOT.TLegend(0.7, 0.15, 1.0, 0.25) 
    l.AddEntry( tot_xsec_toy, "toy" )
    l.AddEntry( tot_xsec_delphes, "delphes" )
    l.Draw()
    c1.Print(os.path.join(plot_directory, "log", "xsec_%s_%s.png" %(WC, plot.name)))


#plot_phi_subtr.histos = plot_phi_subtr.histos[1:]

#drawPlots(draw_plots, subDirectory = subDirectory)

logger.info( "Done with %s", selection )

syncer.sync()

