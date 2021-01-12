#!/usr/bin/env python
# Standard imports and batch mode
import ROOT, os, itertools, operator
#ROOT.gROOT.SetBatch(True)
import copy
from math                           import sqrt, cos, sin, pi, isnan, sinh, cosh, log

# Analysis
import Analysis.Tools.syncer        as syncer

# RootTools
from RootTools.core.standard        import *

# TMB
from TMB.Tools.user                 import plot_directory
from TMB.Tools.helpers              import getCollection
from TMB.Tools.WeightInfo           import WeightInfo
from TMB.Tools.cutInterpreter       import cutInterpreter
from TMB.Tools.VV_angles            import VV_angles

# Arguments
import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--logLevel',           action='store',      default='INFO',          nargs='?', choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'], help="Log level for logging")
argParser.add_argument('--plot_directory',     action='store',      default='fisher')
argParser.add_argument('--selection',          action='store',      default=None)
argParser.add_argument('--sample',             action='store',      default='WZ')
argParser.add_argument('--WC',                 action='store',      default=None)
argParser.add_argument('--WCval',              action='store',      type=float,    default=0.0,  help='Value of the Wilson coefficient')
argParser.add_argument('--small',                                   action='store_true',     help='Run only on a small subset of the data?')
args = argParser.parse_args()

# Logger
import TMB.Tools.logger as logger
import RootTools.core.logger as logger_rt
logger    = logger.get_logger(   args.logLevel, logFile = None)
logger_rt = logger_rt.get_logger(args.logLevel, logFile = None)

if args.small: args.plot_directory += "_small"

# Import samples
from TMB.Samples.gen_pp import *
sample = eval(args.sample) 

# WeightInfo
w = WeightInfo(sample.reweight_pkl)
w.set_order(2)

if args.WC is None:
    wc_kwargs = {}
else:
    wc_kwargs = {args.WC:args.WCval}
sample.weight =  w.get_weight_func(**wc_kwargs)

stack = Stack([sample])

# Read variables and sequences
read_variables = [
    "genMet_pt/F", "genMet_phi/F", "np/I",
    "ngenJet/I", "genJet[pt/F,eta/F,phi/F,matchBParton/I]", 
    "ngenLep/I", "genLep[pt/F,eta/F,phi/F,pdgId/I,mother_pdgId/I]", 
    "ngenTop/I", "genTop[pt/F,eta/F,phi/F]", 
    "ngenZ/I", "genZ[pt/F,phi/F,eta/F,daughter_pdgId/I,l1_index/I,l2_index/I]",
    "ngenW/I", "genW[pt/F,phi/F,eta/F,daughter_pdgId/I,l1_index/I,l2_index/I]",
    "ngenPhoton/I", "genPhoton[pt/F,phi/F,eta/F]"
]
read_variables.append( VectorTreeVariable.fromString('p[C/F]', nMax=2000) )

preselection = [ 
    #( "SMP-20-005",  "genPhoton_pt[0]>300&&genMet_pt>80&&genLep_pt[genW_l1_index[0]]>80&&sqrt(acos(cos(genLep_phi[genW_l1_index[0]]-genPhoton_phi[0]))**2+(genLep_eta[genW_l1_index[0]]-genPhoton_eta[0])**2)>3.0"),
    #( "SMP-20-005-light",  "genPhoton_pt[0]>150&genMet_pt>30&&genLep_pt[genW_l1_index[0]]>30&&sqrt(acos(cos(genLep_phi[genW_l1_index[0]]-genPhoton_phi[0]))**2+(genLep_eta[genW_l1_index[0]]-genPhoton_eta[0])**2)>3.0"),
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

reader = sample.treeReader( variables = map( lambda v: TreeVariable.fromString(v) if type(v)==type("") else v, read_variables) )
reader.start()
fi = None
while reader.run():
    fi_ = w.get_fisherInformation_matrix( 
        map( operator.itemgetter('C'), getCollection( reader.event, "p", ["C"], "np" ) ),
        w.variables, **wc_kwargs)[1]
    if fi is None:
        fi= fi_
    else:
        fi+=fi_

FI = ROOT.TH2D('FI', 'FI', len(w.variables), 0, len(w.variables), len(w.variables), 0, len(w.variables))
for i in range( FI.GetNbinsX() ):
    FI.GetXaxis().SetBinLabel(i+1, w.variables[i] )
    FI.GetYaxis().SetBinLabel(i+1, w.variables[i] )
    for j in range( FI.GetNbinsY() ):
        FI.SetBinContent( i+1, j+1, abs(fi[i][j]) )

if args.WC is not None:
    name = "FI_%s_%3.2f"%(args.WC, args.WCval)
else:
    name = "FI_SM"

plot = Plot2D.fromHisto(name, [[FI]] , texX = "", texY = "")
plot.drawOption = "COLZ,TEXT"
plotting.draw2D(plot,
    plot_directory = os.path.join( plot_directory, args.plot_directory, sample.name), 
    logZ=True, zRange = (10**-4.5,1),
    canvasModifications = [lambda c: ROOT.tdrStyle.SetPaintTextFormat("3.2f")])
