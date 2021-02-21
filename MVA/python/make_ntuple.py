#!/usr/bin/env python

# General
import os
import ROOT

# Analysis
#import Analysis.Tools.syncer
# RootTools
from RootTools.core.standard import *
# TopEFT
from TMB.Tools.cutInterpreter    import cutInterpreter

# MVA configuration
import TMB.MVA.configs  as configs 

import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--selection',          action='store', type=str,   default='nLep1-nPhoton1')
argParser.add_argument('--sample',             action='store', type=str,   default='ttG_noFullyHad_fast')
argParser.add_argument('--config',             action='store', type=str,   default='ttG')
argParser.add_argument('--output_directory',   action='store', type=str,   default='.')
argParser.add_argument('--small',              action='store_true')

args = argParser.parse_args()

#Logger
import TMB.Tools.logger as logger
logger = logger.get_logger("INFO", logFile = None )
import Analysis.Tools.logger as logger_an
logger_an = logger_an.get_logger("INFO", logFile = None )
#import RootTools.core.logger as logger_rt
#logger_rt = logger_rt.get_logger("DEBUG", logFile = None )

# get samples
import TMB.Samples.pp_tWZ as samples
sample = getattr(samples, args.sample)

if args.small:
    sample.reduceFiles(to=1)

# selection
if args.selection == None:
    selectionString = "(1)"
else:
    selectionString = cutInterpreter.cutString( args.selection )
sample.setSelectionString( selectionString )

count  = int(sample.getYieldFromDraw( weightString="(1)" )["val"])
logger.info( "Found %i events for sample %s", count, sample.name )

#config
config = getattr( configs, args.config)

# where the output goes
output_file  = os.path.join( args.output_directory, sample.name + ("_small" if args.small else "") + ".root" )

# reader
reader = sample.treeReader( \
    #variables = map( TreeVariable.fromString, config.read_variables),
    variables =  config.read_variables,
    sequence  = config.sequence,
    )
reader.start()

#filler
def filler( event ):

    r = reader.event

    # fill extra variables
    #event.isTraining = isTraining
    #event.isSignal   = isSignal
    # write mva variables
    for name, func in config.all_mva_variables.iteritems():
#                setattr( event, name, func(reader.event) )
        setattr( event, name, func(reader.event, sample=None) )

    if hasattr(config, "FIs"):
        for FI_name, FI in config.FIs.iteritems():
            #print( event, 'mva_'+FI_name, FI['func']( [r.p_C[i] for i in range(r.np) ] ) )
            setattr( event, 'FI_'+FI_name, FI['func']( [r.p_C[i] for i in range(r.np) ] )[1][0][0] )

# Create a maker. Maker class will be compiled. 

mva_variables = ["%s/F"%var for var in config.all_mva_variables.keys()]
if hasattr( config, "FIs"):
    FI_variables = ["FI_%s/F"%var for var in config.FIs.keys() ]
else:
    FI_variables = [] 

maker = TreeMaker(
    sequence  = [ filler ],
    variables = map(TreeVariable.fromString, 
#          ["isTraining/I", "isSignal/I"] + 
          mva_variables+FI_variables,  
        ),
    treeName = "Events"
    )
maker.start()

logger.info( "Starting event loop" )
counter=0
while reader.run():

    #for func in config.sequence:
    #    func(reader.event)

    ## determine whether training or test
    #isTraining = self.samples[i_sample].training_test_list.pop(0)
    #isSignal   = (i_sample == 0)

    maker.run()
    counter += 1
    if counter%10000 == 0:
        logger.info("Written %i events.", counter)

nEventsTotal = maker.tree.GetEntries()

tmp_directory = ROOT.gDirectory
dirname = os.path.dirname(output_file)
if not os.path.exists(dirname):
    os.makedirs(dirname)

outputfile = ROOT.TFile.Open(output_file, 'recreate')
maker.tree.Write()
outputfile.Close()
tmp_directory.cd()
logger.info( "Written %s", output_file)
#
#      # Destroy the TTree
maker.clear()
logger.info( "Written %i events to %s",  nEventsTotal, output_file )
