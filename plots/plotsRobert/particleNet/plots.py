#!/usr/bin/env python

# Standard imports and batch mode
import ROOT, os, itertools
#ROOT.gROOT.SetBatch(True)
ROOT.gROOT.SetBatch(True)

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
argParser.add_argument('--selection',          action='store',      default=None)
argParser.add_argument('--sample',             action='store',      default='WH', choices = ['WH', 'ZH', 'WH_nlo', 'ZH_nlo'])
argParser.add_argument('--small',                                   action='store_true',     help='Run only on a small subset of the data?')
argParser.add_argument('--show_derivatives',                        action='store_true',     help='Show also the derivatives?')

args = argParser.parse_args()

# Logger'singlelep-WHJet' if sample.name=='WH' else 'dilep-ZHJet-onZ'
import TMB.Tools.logger as logger
import RootTools.core.logger as logger_rt
logger    = logger.get_logger(   args.logLevel, logFile = None)
logger_rt = logger_rt.get_logger(args.logLevel, logFile = None)

plot_directory = os.path.join(plot_directory, args.plot_directory,  args.sample )
if args.small: plot_directory += "_small"

# Import samples
import TMB.Samples.genTopJets_v1 as samples
    
sample = getattr( samples, args.sample)
 
# WeightInfo
sample.weightInfo = WeightInfo(sample.reweight_pkl)
sample.weightInfo.set_order(2)
sample.read_variables = [VectorTreeVariable.fromString( "p[C/F]", nMax=200 )]

eft_configs = [
    {'color':ROOT.kBlack,       'param':{}, 'tex':"SM"},
    {'color':ROOT.kMagenta-4,   'param':{'cHWtil':-1},  'tex':"c_{H#tilde{W}}=-1",  'binning':[20,0,1.5]},
    {'color':ROOT.kMagenta+2,   'param':{'cHWtil':1},   'tex':"c_{H#tilde{W}}=1",   'binning':[20,0,1.5]},
    {'color':ROOT.kGreen-4,     'param':{'cHW':-1},     'tex':"c_{HW}=-1",          'binning':[20,0,1.5]},
    {'color':ROOT.kGreen+2,     'param':{'cHW':1},      'tex':"c_{HW}=1",           'binning':[20,0,1.5]},
    {'color':ROOT.kBlue-4,      'param':{'cHj3':-.1},   'tex':"c_{HQ}^{(3)}=-0.1",  'binning':[20,0,50]},
    {'color':ROOT.kBlue+2,      'param':{'cHj3':.1},    'tex':"c_{HQ}^{(3)}=0.1",   'binning':[20,0,50]},
    ]

for eft in eft_configs:
    eft['func'] = sample.weightInfo.get_weight_func(**eft['param']) 
    eft['name'] = "_".join( ["sample"] + ( ["SM"] if len(eft['param'])==0 else [ "_".join([key, str(val)]) for key, val in sorted(eft['param'].iteritems())] ) )

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
    der['func'] = sample.weightInfo.get_diff_weight_func(der['der'])
    der['name'] = "_".join( ["derivative"] + list(der['der']) )

sequence = []
def make_eft_weights( event, sample):
    SM_ref_weight         = eft_configs[0]['func'](event, sample)
    event.eft_weights     = [1] + [eft['func'](event, sample)/SM_ref_weight for eft in eft_configs[1:]]
    event.eft_derivatives = [der['func'](event, sample)/SM_ref_weight for der in eft_derivatives]

stack = Stack( [sample] )

sequence.append( make_eft_weights )

eft_weights = [] 
if not args.no_bkgs:
    eft_weights =  [[]]
    for sample in stack.samples:
        if not hasattr( sample, "style" ):
            sample.style = styles.fillStyle(sample.color)
        eft_weights[0].append( None )

for i_eft, eft in enumerate(eft_configs):
    stack.append( [sample] )
    eft_weights.append( [lambda event, sample, i_eft=i_eft: event.eft_weights[i_eft]] )

for i_eft, eft in enumerate(eft_derivatives):
    stack.append( [sample] )
    eft_weights.append( [lambda event, sample, i_eft=i_eft: event.eft_derivatives[i_eft]] )

#if args.sign_reweight:
#    weight_branches = ["lumiweight1fb", "sign_reweight"]
#else:
#    #lumi_weight = lambda event, sample: lumi*event.lumiweight1fb#*sin(2*event.VV_angles['Theta'])*sin(2*event.VV_angles['theta_V1'])
#    weight_branches = ["lumiweight1fb"]
##lumi  = 300
#def weight_getter( branches ):
#    getters = [ operator.attrgetter(branch) for branch in branches ]
#    def getter( event, sample ):
##        for b, g in zip( branches, getters ):
##            print b, g(event)
##        print
#        return reduce( operator.mul , [ g(event) for g in getters ], lumi ) 
#    return getter
#
#for sample in stack.samples:
#    sample.weight = weight_getter( weight_branches + ( ["combinatoricalBTagWeight2b"] if "DY" in sample.name and args.combinatoricalBTags else [] ) )

# Read variables and sequences
#jetVars          = ['pt/F', 'eta/F', 'phi/F', 'bTag/F', 'bTagPhys/I']
#jetVarNames      = [x.split('/')[0] for x in jetVars]

read_variables = [\
    "partonTop_pt/F", "partonTop_eta/F", "partonTop_phi/F", "partonTop_mass/F", "partonTop_pdgId/I", 
    "genQ1_pt/F", "genQ1_eta/F", "genQ1_phi/F", "genQ1_mass/F", "genQ1_pdgId/I", 
    "genQ2_pt/F", "genQ2_eta/F", "genQ2_phi/F", "genQ2_mass/F", "genQ2_pdgId/I", 
    "genb_pt/F", "genQ2_eta/F", "genQ2_phi/F", "genQ2_mass/F", "genQ2_pdgId/I", 
    "genW_pt/F", "genW_eta/F", "genW_phi/F", "genW_mass/F", "genW_pdgId/I", 
    "genTop_pt/F", "genTop_eta/F", "genTop_phi/F", "genTop_mass/F", "genTop_pdgId/I", 
    "genJet_pt/F", "genJet_eta/F", "genJet_phi/F", "genJet_mass/F", "genJet_nConstituents/I", 
    "ne/I", "nmu/I", "nch/I", "nph/I", "nchh/I", "nneh/I",
]

preselection = [ 
]

selectionString  = "&&".join( [ c[1] for c in preselection] + ([cutInterpreter.cutString(args.selection)] if args.selection is not None else []))
subDirectory     =  '-'.join( [ c[0] for c in preselection] + ([args.selection] if args.selection is not None else []))
if subDirectory  == '': 
    subDirectory = 'inc'

for sample in stack.samples:
    if selectionString != "":
        sample.addSelectionString( selectionString )
    if args.small:
        sample.reduceFiles( to = 1 )

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

#for eft_config in eft_configs:
#    param = eft_config['param']
#    if len(param)!=1: continue
#    wc_, val_ = list(param.iteritems())[0]
#    name =  "opt_%s_%f"%(wc_, val_)
#    plots.append(Plot(
#        texX = "q(%s=%3.2f)"%(wc_, val_), texY = 'Number of Events',
#        name =  name, 
#        attribute = lambda event, sample, disc_name=name: getattr( event, disc_name ),
#        binning=eft_config['binning'],
#    ))

# Text on the plots
def drawObjects( hasData = False ):
    tex = ROOT.TLatex()
    tex.SetNDC()
    tex.SetTextSize(0.04)
    tex.SetTextAlign(11) # align right
    lines = [
      #(0.15, 0.95, 'CMS Preliminary' if hasData else "Delphes Simulation"), 
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

