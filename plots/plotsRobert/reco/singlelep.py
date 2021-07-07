#!/usr/bin/env python
''' Analysis script for standard plots
'''
#
# Standard imports and batch mode
#
import ROOT, os
ROOT.gROOT.SetBatch(True)
c1 = ROOT.TCanvas() # do this to avoid version conflict in png.h with keras import ...
c1.Draw()
c1.Print('/tmp/delete.png')

import itertools
import copy
import array
import operator
from math                                import sqrt, cos, sin, pi, atan2, cosh, log

# RootTools
from RootTools.core.standard             import *

# tWZ
from TMB.Tools.user                      import plot_directory
from TMB.Tools.cutInterpreter            import cutInterpreter
from tWZ.Tools.objectSelection           import lepString # probably will merge TMB and tWZ repos 
# Analysis
from Analysis.Tools.helpers              import deltaPhi, deltaR
from Analysis.Tools.puProfileCache       import *
from Analysis.Tools.puReweighting        import getReweightingFunction
from Analysis.Tools.WeightInfo                import WeightInfo
import Analysis.Tools.syncer             as     syncer
import numpy as np

# Arguments
import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--logLevel',       action='store',      default='INFO', nargs='?', choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'], help="Log level for logging")
argParser.add_argument('--small',                             action='store_true', help='Run only on a small subset of the data?', )
#argParser.add_argument('--sorting',                           action='store', default=None, choices=[None, "forDYMB"],  help='Sort histos?', )
argParser.add_argument('--plot_directory', action='store', default='FI-test-v6')
argParser.add_argument('--WC',                 action='store',      default='cpQ3')
argParser.add_argument('--WCval',              action='store',      nargs = '*',             type=float,    default=[1.0],  help='Values of the Wilson coefficient for the distribution.')
argParser.add_argument('--WCval_FI',           action='store',      nargs = '*',             type=float,    default=[0.0, 1.0],  help='Values of the Wilson coefficient to show FI for.')
argParser.add_argument('--era',            action='store', type=str, default="Autumn18")
argParser.add_argument('--sample',        action='store', type=str, default="ttG_noFullyHad_fast")
argParser.add_argument('--selection',      action='store', default='singlelep-photon')
args = argParser.parse_args()

# Logger
import TMB.Tools.logger as logger
import RootTools.core.logger as logger_rt
logger    = logger.get_logger(   args.logLevel, logFile = None)
logger_rt = logger_rt.get_logger(args.logLevel, logFile = None)

if args.small:                        args.plot_directory += "_small"

import TMB.Samples.pp_tWZ_v6 as samples

sample = getattr( samples, args.sample )

lumi_scale = 137

sample.weight = lambda event, sample: lumi_scale*event.weight 


if args.small:
    sample.reduceFiles( to = 1 )

# WeightInfo
w = WeightInfo(sample.reweight_pkl)
w.set_order(2)

# define which Wilson coefficients to plot
FIs = []

params =  [ {'legendText':'SM',  'color':ROOT.kBlue, 'WC':{}} ]
params += [ {'legendText':'%s %3.2f'%(args.WC, wc), 'color':ROOT.kOrange+i_wc,  'WC':{args.WC:wc} } for i_wc, wc in enumerate(args.WCval)]

for i_param, param in enumerate(params):
    param['sample'] = sample
    param['style']  = styles.lineStyle( param['color'] )

stack = Stack(*[ [ param['sample'] ] for param in params ] )
weight= [ [ w.get_weight_func(**param['WC']) ] for param in params ]
FIs.append( ( ROOT.kGray+1, "WC@SM",   w.get_fisher_weight_string(args.WC,args.WC, **{args.WC:0})) )
for i_WCval, WCval in enumerate(args.WCval_FI):
    FIs.append( ( ROOT.kGray+i_WCval,   "WC@WC=%3.2f"%WCval, w.get_fisher_weight_string(args.WC,args.WC, **{args.WC:WCval})) )

import re
# make stack and construct weights
def coeff_getter( coeff):
    def getter_( event, sample):
        #print "Filling with", coeff, event.p_C[coeff]
        return event.p_C[coeff]
    return getter_

def add_fisher_plot( plot, fisher_string, color, legendText):
    # the coefficients we need
    required_coefficients = map( lambda s: int(s.replace('[','').replace(']','')), list(set(re.findall(r'\[[0-9][0-9]*\]', fisher_string ) )) )
    # copy the plot
    fisher_plot = copy.deepcopy( plot )
    # modify weights etc
    fisher_plot.name+="_fisher_coeffs"
    fisher_plot.weight = [ [ coeff_getter(coeff) ] for coeff in required_coefficients ]
    fisher_plot.stack  = Stack( *[[sample] for coeff in required_coefficients] )
    # for computing the FI
    fisher_plot.coefficients = required_coefficients
    fisher_plot.formula      = re.sub(r'p_C\[([0-9][0-9]*)\]', r'{lambda_\1}', fisher_string)
    # pass through information for plotting
    fisher_plot.color        = color
    fisher_plot.legendText   = legendText
    # add the fisher_plot to the original plot so we know which is which
    if hasattr( plot, "fisher_plots"):
        plot.fisher_plots.append( fisher_plot )
    else:
        plot.fisher_plots = [fisher_plot]
    return fisher_plot

            
# Read variables and sequences
sequence       = []

read_variables = [
    "weight/F", "year/I", "met_pt/F", "met_phi/F", "nBTag/I", "nJetGood/I", "PV_npvsGood/I",
    "l1_pt/F", "l1_eta/F" , "l1_phi/F", "l1_mvaTOP/F", "l1_mvaTOPWP/I", "l1_index/I", 
    #"l2_pt/F", "l2_eta/F" , "l2_phi/F", "l2_mvaTOP/F", "l2_mvaTOPWP/I", "l2_index/I",
    #"l3_pt/F", "l3_eta/F" , "l3_phi/F", "l3_mvaTOP/F", "l3_mvaTOPWP/I", "l3_index/I",
    "JetGood[pt/F,eta/F,phi/F]",
    "lep[pt/F,eta/F,phi/F,pdgId/I,muIndex/I,eleIndex/I]",
#    "Z1_l1_index/I", "Z1_l2_index/I", "nonZ1_l1_index/I", "nonZ1_l2_index/I", 
#    "Z1_phi/F", "Z1_pt/F", "Z1_mass/F", "Z1_cosThetaStar/F", "Z1_eta/F", "Z1_lldPhi/F", "Z1_lldR/F",
    "Muon[pt/F,eta/F,phi/F,dxy/F,dz/F,ip3d/F,sip3d/F,jetRelIso/F,miniPFRelIso_all/F,pfRelIso03_all/F,mvaTOP/F,mvaTTH/F,pdgId/I,segmentComp/F,nStations/I,nTrackerLayers/I]",
    "Electron[pt/F,eta/F,phi/F,dxy/F,dz/F,ip3d/F,sip3d/F,jetRelIso/F,miniPFRelIso_all/F,pfRelIso03_all/F,mvaTOP/F,mvaTTH/F,pdgId/I,vidNestedWPBitmap/I]",
    "np/I", VectorTreeVariable.fromString("p[C/F]",nMax=500),
    "photon_pt/F",
    "photon_eta/F",
    "photon_phi/F",
    "photonJetdR/F", "photonLepdR/F",
]

read_variables_MC = ['reweightBTag_SF/F', 'reweightPU/F', 'reweightL1Prefire/F', 'reweightLeptonSF/F'] #'reweightTrigger/F']
# define 3l selections

#MVA
##from Analysis.TMVA.Reader    import Reader
#import TMB.MVA.configs.ttG as mva_config 
#
#sequence.extend( mva_config.sequence )
#read_variables.extend( mva_config.read_variables )

def discriminator_getter(name):
    def _disc_getter( event, sample ):
        return getattr( event, name )
    return _disc_getter

mu_string  = lepString('mu','VL')
ele_string = lepString('ele','VL')
def getLeptonSelection():
    return "Sum$({mu_string})+Sum$({ele_string})==1".format(mu_string=mu_string,ele_string=ele_string)

# Getter functor for lepton quantities
def lep_getter( branch, index, abs_pdg = None, functor = None, debug=False):
    if functor is not None:
        if abs_pdg == 13:
            def func_( event, sample ):
                if debug:
                    print "Returning", "Muon_%s"%branch, index, abs_pdg, "functor", functor, "result",
                    print functor(getattr( event, "Muon_%s"%branch )[event.lep_muIndex[index]]) if abs(event.lep_pdgId[index])==abs_pdg else float('nan')
                return functor(getattr( event, "Muon_%s"%branch )[event.lep_muIndex[index]]) if abs(event.lep_pdgId[index])==abs_pdg else float('nan')
        else:
            def func_( event, sample ):
                return functor(getattr( event, "Electron_%s"%branch )[event.lep_eleIndex[index]]) if abs(event.lep_pdgId[index])==abs_pdg else float('nan')
    else:
        if abs_pdg == 13:
            def func_( event, sample ):
                if debug:
                    print "Returning", "Muon_%s"%branch, index, abs_pdg, "functor", functor, "result",
                    print getattr( event, "Muon_%s"%branch )[event.lep_muIndex[index]] if abs(event.lep_pdgId[index])==abs_pdg else float('nan')
                return getattr( event, "Muon_%s"%branch )[event.lep_muIndex[index]] if abs(event.lep_pdgId[index])==abs_pdg else float('nan')
        else:
            def func_( event, sample ):
                return getattr( event, "Electron_%s"%branch )[event.lep_eleIndex[index]] if abs(event.lep_pdgId[index])==abs_pdg else float('nan')
    return func_

yields     = {}
allPlots   = {}

#yt_TWZ_filter.scale = lumi_scale * 1.07314

# Use some defaults
Plot.setDefaults(stack = stack, weight = weight, selectionString = cutInterpreter.cutString(args.selection))

plots        = []
fisher_plots = []

plots.append(Plot(
  name = 'nVtxs', texX = 'vertex multiplicity', texY = 'Number of Events',
  attribute = TreeVariable.fromString( "PV_npvsGood/I" ),
  binning=[50,0,50],
  addOverFlowBin='upper',
))

plots.append(Plot(
    name = 'M3',
    texX = 'M_{3} (GeV)', texY = 'Number of Events / 10 GeV',
    attribute = TreeVariable.fromString( "m3/F" ),
    binning=[30,0,300],
    addOverFlowBin='upper',
))

plots.append(Plot(
    name = 'photonJetdR',
    texX = '#Delta R(#gamma, jets)', texY = 'Number of Events',
    attribute = TreeVariable.fromString( "photonJetdR/F" ),
    binning=[30,0,6],
    addOverFlowBin='upper',
))

plots.append(Plot(
    name = 'photonLepdR',
    texX = '#Delta R(#gamma, leptons)', texY = 'Number of Events',
    attribute = TreeVariable.fromString( "photonLepdR/F" ),
    binning=[30,0,6],
    addOverFlowBin='upper',
))

plots.append(Plot(
    name = 'mT',
    texX = 'M_{T} (GeV)', texY = 'Number of Events / 10 GeV',
    attribute = lambda event, sample: sqrt(2.*event.met_pt*event.l1_pt*(1-cos(event.met_phi-event.l1_phi))),
    binning=[30,0,300],
    addOverFlowBin='upper',
))

plots.append(Plot(
    name = 'photon_pt',
    texX = 'p_{T}(#gamma) (GeV)', texY = 'Number of Events / 20 GeV',
    attribute = lambda event, sample:event.photon_pt,
    binning=[15,0,300],
    addOverFlowBin='upper',
))
for color, legendText, fisher_string in FIs:
    fisher_plots.append( add_fisher_plot( plots[-1], fisher_string, color, legendText ) )

plots.append(Plot(
    name = 'l1_pt',
    texX = 'p_{T}(l_{1}) (GeV)', texY = 'Number of Events / 20 GeV',
    attribute = lambda event, sample:event.l1_pt,
    binning=[15,0,300],
    addOverFlowBin='upper',
))
for color, legendText, fisher_string in FIs:
    fisher_plots.append( add_fisher_plot( plots[-1], fisher_string, color, legendText ) )

plots.append(Plot(
    name = 'l1_eta',
    texX = '#eta(l_{1})', texY = 'Number of Events',
    attribute = lambda event, sample: event.l1_eta,
    binning=[20,-3,3],
))

plots.append(Plot(
    name = 'l1_mvaTOP',
    texX = 'MVA_{TOP}(l_{1})', texY = 'Number of Events',
    attribute = lambda event, sample: event.l1_mvaTOP,
    binning=[20,-1,1],
))

plots.append(Plot(
    name = 'l1_mvaTOPWP',
    texX = 'MVA_{TOP}(l_{1}) WP', texY = 'Number of Events',
    attribute = lambda event, sample: event.l1_mvaTOPWP,
    binning=[5,0,5],
))

plots.append(Plot(
    texX = 'E_{T}^{miss} (GeV)', texY = 'Number of Events / 20 GeV',
    attribute = TreeVariable.fromString( "met_pt/F" ),
    binning=[400/20,0,400],
    addOverFlowBin='upper',
))

plots.append(Plot(
  texX = 'N_{jets}', texY = 'Number of Events',
  attribute = TreeVariable.fromString( "nJetGood/I" ), #nJetSelected
  binning=[8,-0.5,7.5],
))

plots.append(Plot(
  texX = 'N_{b-tag}', texY = 'Number of Events',
  attribute = TreeVariable.fromString( "nBTag/I" ), #nJetSelected
  binning=[4,-0.5,3.5],
))

plots.append(Plot(
  texX = 'p_{T}(leading jet) (GeV)', texY = 'Number of Events / 30 GeV',
  name = 'jet0_pt', attribute = lambda event, sample: event.JetGood_pt[0],
  binning=[600/30,0,600],
))

plots.append(Plot(
  texX = 'p_{T}(subleading jet) (GeV)', texY = 'Number of Events / 30 GeV',
  name = 'jet1_pt', attribute = lambda event, sample: event.JetGood_pt[1],
  binning=[600/30,0,600],
))

# Text on the plots
def drawObjects( hasData = False ):
    tex = ROOT.TLatex()
    tex.SetNDC()
    tex.SetTextSize(0.04)
    tex.SetTextAlign(11) # align right
    lines = [
      (0.15, 0.95, 'CMS Preliminary' if hasData else "CMS Simulation"),
      (0.45, 0.95, 'L=%3.1f fb{}^{-1} (13 TeV)' % lumi_scale),
    ]
    return [tex.DrawLatex(*l) for l in lines]

# draw function for plots
def drawPlots(plots):
  for log in [False, True]:
    plot_directory_ = os.path.join(plot_directory, args.plot_directory, args.sample, args.WC, args.selection, ("log" if log else "lin") )
    for plot in plots:
      if not max(l[0].GetMaximum() for l in plot.histos): continue # Empty plot

      len_FI = len(plot.fisher_plots) if hasattr(plot, "fisher_plots") else 0
      plotting.draw(plot,
        plot_directory = plot_directory_,
        ratio = {'histos':[(i,0) for i in range(1,len(plot.histos)-len_FI)], 'yRange':(0.1,1.9)} ,
        logX = False, logY = log, sorting = False,
        yRange = (0.03, "auto") if log else "auto",
        scaling = {},
        legend =  ( (0.17,0.9-0.05*sum(map(len, plot.histos))/2,1.,0.9), 2),
        drawObjects = drawObjects( ),
        copyIndexPHP = True,
      )

plotting.fill(plots+fisher_plots, read_variables = read_variables, sequence = sequence, max_events = -1 if args.small else -1)

for plot in plots:
    for i_h, hl in enumerate(plot.histos):
        # dress up
        hl[0].legendText = params[i_h]['legendText']
        hl[0].style = params[i_h]['style']

    # calculate & add FI histo
    if hasattr(plot, "fisher_plots" ):
        for fisher_plot in plot.fisher_plots:

            # make empty histo
            FI = plot.histos[0][0].Clone()
            FI.Reset()

            # calculate FI in each bin
            for i_bin in range(1, FI.GetNbinsX()+1):
                yields = {'lambda_%i'%coeff : fisher_plot.histos[i_coeff][0].GetBinContent(i_bin) for i_coeff, coeff in enumerate(fisher_plot.coefficients)}
                try:
                    FI.SetBinContent( i_bin,  eval(fisher_plot.formula.format(**yields)) )
                except ZeroDivisionError:
                    pass

            # dress up
            FI.legendText = fisher_plot.legendText
            FI.style = styles.lineStyle( fisher_plot.color )

            # scale the FI & indicate log(I) in the plot
            if FI.Integral()>0:
                FI.legendText += "(%3.2f)" % log(FI.Integral(),10)
                FI.Scale(plot.histos[0][0].Integral()/FI.Integral())

            # add the FI histos back to the plot and fake the stack so that the draw command does not get confused
            plot.histos.append( [FI] )
            stack.append( [sample] )

drawPlots(plots)

syncer.sync()

logger.info( "Done with prefix %s and selectionString %s", args.selection, cutInterpreter.cutString(args.selection) )
