#!/usr/bin/env python

# General
import os
import ROOT

# Analysis
#import Analysis.Tools.syncer
# RootTools
from RootTools.core.standard import *

# TMB
from TMB.Tools.cutInterpreter    import cutInterpreter
from TMB.Tools.helpers import getVarValue, getObjDict

# MVA configuration
import TMB.MVA.configs  as configs 

import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--logLevel', action='store', nargs='?',  choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'],   default='INFO', help="Log level for logging" )
argParser.add_argument('--selection',          action='store', type=str,   default='singlelep-photon')
argParser.add_argument('--sample',             action='store', type=str,   default='ttG_noFullyHad_fast')
argParser.add_argument('--sample_file',        action='store', type=str,   default='$CMSSW_BASE/src/TMB/Samples/python/pp_tWZ_v6.py')
argParser.add_argument('--config',             action='store', type=str,   default='ttG')
argParser.add_argument('--output_directory',   action='store', type=str,   default='.')
argParser.add_argument('--small',              action='store_true')

args = argParser.parse_args()

#Logger
import TMB.Tools.logger as logger
logger = logger.get_logger(args.logLevel, logFile = None )
import Analysis.Tools.logger as logger_an
logger_an = logger_an.get_logger(args.logLevel, logFile = None )
import RootTools.core.logger as logger_rt
logger_rt = logger_rt.get_logger(args.logLevel, logFile = None )

#import TMB.Samples.pp_tWZ_v6 as samples
# get samples
import imp
samples = imp.load_source('samples', os.path.expandvars(args.sample_file))
sample = getattr(samples, args.sample)

subDir = args.config
if args.small:
    sample.reduceFiles(to=1)
    subDir += '_small'

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
output_file  = os.path.join( args.output_directory, "MVA-training", subDir, sample.name, args.selection, sample.name + ".root" )

# reader
reader = sample.treeReader( \
    #variables = map( TreeVariable.fromString, config.read_variables),
    variables =  config.read_variables,
    sequence  = config.sequence,
    )
reader.start()

def fill_vector_collection( event, collection_name, collection_varnames, objects, maxN = 100):
    setattr( event, "n"+collection_name, len(objects) )
    for i_obj, obj in enumerate(objects[:maxN]):
        for var in collection_varnames:
            if var in obj.keys():
                if type(obj[var]) == type("string"):
                    obj[var] = int(ord(obj[var]))
                if type(obj[var]) == type(True):
                    obj[var] = int(obj[var])
                getattr(event, collection_name+"_"+var)[i_obj] = obj[var]

#filler
def filler( event ):

    r = reader.event

    # copy scalar variables
    for name, func in config.all_mva_variables.iteritems():
        setattr( event, name, func(reader.event, sample=None) )

    # copy vector variables
    for name, vector_var in config.mva_vector_variables.iteritems():
        objs = [ getObjDict( reader.event, vector_var['name']+'_', vector_var['varnames'], i ) for i in range(int(getVarValue(reader.event, 'n'+vector_var['name']))) ] 
        objs = filter( vector_var['selector'], objs ) 

        fill_vector_collection( event, name, vector_var['varnames'], objs )

    # fill FIs
    if hasattr(config, "FIs"):
        for FI_name, FI in config.FIs.iteritems():
            #print( event, 'mva_'+FI_name, FI['func']( [r.p_C[i] for i in range(r.np) ] ) )
            setattr( event, 'FI_'+FI_name, FI['func']( [r.p_C[i] for i in range(r.np) ] )[1][0][0] )

# Create a maker. Maker class will be compiled. 

# scalar variables
mva_variables = ["%s/F"%var for var in config.all_mva_variables.keys()]

# vector variables, if any
for name, vector_var in config.mva_vector_variables.iteritems():
    #mva_variables.append( VectorTreeVariable.fromString(name+'['+','.join(vector_var['vars'])+']') )
    mva_variables.append( name+'['+','.join(vector_var['vars'])+']' )
# FIs
if hasattr( config, "FIs"):
    FI_variables = ["FI_%s/F"%var for var in config.FIs.keys() ]
else:
    FI_variables = [] 


maker = TreeMaker(
    sequence  = [ filler ],
    variables = map(TreeVariable.fromString, 
          mva_variables+FI_variables,  
        ),
    treeName = "Events"
    )
maker.start()

logger.info( "Starting event loop" )
counter=0
while reader.run():

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
