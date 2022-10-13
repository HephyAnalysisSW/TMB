#!/usr/bin/env python
''' Make flat ntuple from GEN data tier 
'''
#
# Standard imports and batch mode
#
import ROOT
import os, sys
ROOT.gROOT.SetBatch(True)
import itertools
from math                                import sqrt, cos, sin, pi, acos, cosh, sinh
import imp

#RootTools
from RootTools.core.standard             import *

#Analysis
from Analysis.Tools.WeightInfo        import WeightInfo
from Analysis.Tools.leptonJetArbitration    import cleanJetsAndLeptons

# TMB
from TMB.Tools.user                   import skim_output_directory
from TMB.Tools.GenSearch              import GenSearch
from TMB.Tools.helpers                import deltaPhi, deltaR, deltaR2, cosThetaStar, closestOSDLMassToMZ, checkRootFile
from TMB.Tools.HyperPoly              import HyperPoly
from TMB.Tools.DelphesProducer        import DelphesProducer
from TMB.Tools.genObjectSelection     import isGoodGenJet, isGoodGenLepton, isGoodGenPhoton, genJetId
from TMB.Tools.DelphesObjectSelection import isGoodRecoMuon, isGoodRecoElectron, isGoodRecoLepton, isGoodRecoJet, isGoodRecoPhoton

# CMSSW FastJet, Wrappers
#import fastjet
#ak8 = fastjet.JetDefinition(fastjet.antikt_algorithm, 0.8, fastjet.E_scheme)
softDrop       = ROOT.SoftDropWrapper(0.0, 0.1, 0.8, 200)
nSubjettiness  = ROOT.NsubjettinessWrapper( 1, 0.8, 0, 6 )

ecf  = ROOT.ECFWrapper()
ecfs = [
 ('ecf1',       ( 1, 1., 1., "ECF" )),
 ('ecf2',       ( 2, 1., 1., "ECF" )),
 ('ecf3',       ( 3, 1., 1., "ECF" )),
 ('ecfC1',      ( 1, 1., 1., "C" )),
 ('ecfC2',      ( 2, 1., 1., "C" )),
 ('ecfC3',      ( 3, 1., 1., "C" )),
 ('ecfD',       ( 2, 1., 1., "D" )),
 ('ecfDbeta2',  ( 2, 2., 2., "D" )),
 ('ecfM1',      ( 1, 1., 1., "M" )),
 ('ecfM2',      ( 2, 1., 1., "M" )),
 ('ecfM3',      ( 3, 1., 1., "M" )),
 ('ecfM1beta2', ( 1, 2., 2., "M" )),
 ('ecfM2beta2', ( 2, 2., 2., "M" )),
 ('ecfM3beta2', ( 3, 2., 2., "M" )),
 ('ecfN1',      ( 1, 1., 1., "N" )),
 ('ecfN2',      ( 2, 1., 1., "N" )),
 ('ecfN3',      ( 3, 1., 1., "N" )),
 ('ecfN1beta2', ( 1, 2., 2., "N" )),
 ('ecfN2beta2', ( 2, 2., 2., "N" )),
 ('ecfN3beta2', ( 3, 2., 2., "N" )),
 ('ecfU1',      ( 1, 1., 1., "U" )),
 ('ecfU2',      ( 2, 1., 1., "U" )),
 ('ecfU3',      ( 3, 1., 1., "U" )),
 ('ecfU1beta2', ( 1, 2., 2., "U" )),
 ('ecfU2beta2', ( 2, 2., 2., "U" )),
 ('ecfU3beta2', ( 3, 2., 2., "U" )),
]

for _, args in ecfs:
    ecf.addECF( *args )


#from TMB.Tools.UpgradeJECUncertainty  import UpgradeJECUncertainty 
#
# Arguments
# 
import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--logLevel',           action='store',      default='INFO',          nargs='?', choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'], help="Log level for logging")
argParser.add_argument('--small',              action='store_true', help='Run only on a small subset of the data?')#, default = True)
argParser.add_argument('--overwrite',          action='store',      nargs='?', choices = ['none', 'all', 'target'], default = 'none', help='Overwrite?')#, default = True)
argParser.add_argument('--targetDir',          action='store',      default='v1')
argParser.add_argument('--sample',             action='store',      default='tt1LepHad', help="Name of the sample loaded from fwlite_benchmarks. Only if no inputFiles are specified")
argParser.add_argument('--inputFiles',         action='store',      nargs = '*', default=[])
argParser.add_argument('--delphesEra',         action='store',      default = None, choices = ["RunII", "RunIICentral", "RunIInoDelphesIso", "RunIIPileUp", "PhaseII"], help="specify delphes era")
argParser.add_argument('--addReweights',       action='store_true',   help="Add reweights?")
argParser.add_argument('--nJobs',              action='store',      nargs='?', type=int, default=1,  help="Maximum number of simultaneous jobs.")
argParser.add_argument('--job',                action='store',      nargs='?', type=int, default=0,  help="Run only job i")
argParser.add_argument('--removeDelphesFiles', action='store_true',   help="remove Delphes file after postprocessing?")
argParser.add_argument('--interpolationOrder', action='store',      nargs='?', type=int, default=2,  help="Interpolation order for EFT weights.")
args = argParser.parse_args()

#
# Logger
#
import TMB.Tools.logger as _logger
import RootTools.core.logger as _logger_rt
logger    = _logger.get_logger(   args.logLevel, logFile = None)
logger_rt = _logger_rt.get_logger(args.logLevel, logFile = None)

# Load sample either from 
if len(args.inputFiles)>0:
    logger.info( "Input files found. Ignoring 'sample' argument. Files: %r", args.inputFiles)
    sample = FWLiteSample( args.targetSampleName, args.inputFiles)
else:
    sample_file = "$CMSSW_BASE/python/TMB/Samples/genTopJets_v1.py"
    samples = imp.load_source( "samples", os.path.expandvars( sample_file ) )
    sample = getattr( samples, args.sample )
    logger.debug( 'Loaded sample %s with %i files.', sample.name, len(sample.files) )

maxEvents = -1
if args.small: 
    args.targetDir += "_small"
    maxEvents       = 500 
    sample.files=sample.files[:1]

# output directory
output_directory = os.path.join(skim_output_directory, 'gen', args.targetDir, sample.name) 

if not os.path.exists( output_directory ): 
    try:
        os.makedirs( output_directory )
    except OSError:
        pass
    logger.info( "Created output directory %s", output_directory )

# upgrade JEC
# upgradeJECUncertainty = UpgradeJECUncertainty()

# variables
variables = []

# Load reweight pickle file if supposed to keep weights. 
extra_variables = []
if args.addReweights:
    weightInfo = WeightInfo( sample.reweight_pkl )
    weightInfo.set_order( args.interpolationOrder ) 
    # Determine coefficients for storing in vector
    # Sort Ids wrt to their position in the card file

    # weights from base base points 
    weight_base      = TreeVariable.fromString( "weight[base/F]")
    weight_base.nMax = weightInfo.nid
    extra_variables.append(weight_base)

    # coefficients for the weight parametrization
    param_vector      = TreeVariable.fromString( "p[C/F]" )
    param_vector.nMax = HyperPoly.get_ndof(weightInfo.nvar, args.interpolationOrder)
    hyperPoly         = HyperPoly( args.interpolationOrder )
    extra_variables.append(param_vector)
    extra_variables.append(TreeVariable.fromString( "chi2_ndof/F"))
    def interpret_weight(weight_id):
        str_s = weight_id.rstrip('_nlo').split('_')
        res={}
        for i in range(len(str_s)/2):
            res[str_s[2*i]] = float(str_s[2*i+1].replace('m','-').replace('p','.'))
        return res

    # Suddenly only lower case weight.id ... who on earth does such things?
    weightInfo_data_lower = {k.lower():val for k, val in weightInfo.data.iteritems()}
    weightInfo_data_lower.update(weightInfo.data)

    target_coeff = {}
    for i_comb, comb in enumerate(weightInfo.combinations[1:]):
        variables.append( "target_%s/F"%("_".join(comb)) )
        target_coeff[i_comb] = "_".join(comb)

# Run only job number "args.job" from total of "args.nJobs"
if args.nJobs>1:
    n_files_before = len(sample.files)
    sample = sample.split(args.nJobs)[args.job]
    n_files_after  = len(sample.files)
    logger.info( "Running job %i/%i over %i files from a total of %i.", args.job, args.nJobs, n_files_after, n_files_before)

max_jet_abseta = 5.1

products = {
    'lhe':{'type':'LHEEventProduct', 'label':("externalLHEProducer")},
    'gen':{'type':'GenEventInfoProduct', 'label':'generator'},
    'gp':{'type':'vector<reco::GenParticle>', 'label':("genParticles")},
    'ak8GenJets':{'type':'vector<reco::GenJet>', 'label':("ak8GenJetsNoNu")},
}

def varnames( vec_vars ):
    return [v.split('/')[0] for v in vec_vars.split(',')]

def vecSumPt(*args):
    return sqrt( sum([o['pt']*cos(o['phi']) for o in args],0.)**2 + sum([o['pt']*sin(o['phi']) for o in args],0.)**2 )

def addIndex( collection ):
    for i  in range(len(collection)):
        collection[i]['index'] = i

# EDM standard variables
variables  += ["run/I", "lumi/I", "evt/l"]

if args.delphesEra is not None:
    pass

if args.addReweights:
    variables.append('rw_nominal/F')
    # Lumi weight 1fb / w_0
    variables.append("ref_lumiweight1fb/F")

variables += ["partonTop_pt/F", "partonTop_eta/F", "partonTop_phi/F", "partonTop_mass/F", "partonTop_pdgId/I"]
variables += ["genQ1_pt/F", "genQ1_eta/F", "genQ1_phi/F", "genQ1_mass/F", "genQ1_pdgId/I"]
variables += ["genQ2_pt/F", "genQ2_eta/F", "genQ2_phi/F", "genQ2_mass/F", "genQ2_pdgId/I"]
variables += ["genb_pt/F", "genb_eta/F", "genb_phi/F", "genb_mass/F", "genb_pdgId/I"]
variables += ["genW_pt/F", "genW_eta/F", "genW_phi/F", "genW_mass/F"]
variables += ["genTop_pt/F", "genTop_eta/F", "genTop_phi/F", "genTop_mass/F"]
variables += ["gen_theta/F", "gen_phi/F"]
variables += ["genJet_pt/F", "genJet_eta/F", "genJet_phi/F", "genJet_mass/F", "genJet_nConstituents/I", "genJet_isMuon/I", "genJet_isElectron/I", "genJet_isPhoton/I"]

variables += ['genJet_SDmass/F', 
              'genJet_SDsubjet0_eta/F', 'genJet_SDsubjet0_deltaEta/F', 'genJet_SDsubjet0_phi/F', 'genJet_SDsubjet0_deltaPhi/F', 'genJet_SDsubjet0_deltaR/F', 'genJet_SDsubjet0_mass/F', 
              'genJet_SDsubjet1_eta/F', 'genJet_SDsubjet1_deltaEta/F', 'genJet_SDsubjet1_phi/F', 'genJet_SDsubjet1_deltaPhi/F', 'genJet_SDsubjet1_deltaR/F', 'genJet_SDsubjet1_mass/F', 
              'genJet_tau1/F', 'genJet_tau2/F', 'genJet_tau3/F', 'genJet_tau4/F', 'genJet_tau21/F', 'genJet_tau32/F']

for i_ecf, (name, _) in enumerate( ecfs ):
    variables.append( "genJet_%s/F"%name )

variables += ["dR_genJet_Q1/F", "dR_genJet_Q2/F", "dR_genJet_W/F", "dR_genJet_b/F", "dR_genJet_top/F", "dR_genJet_maxQ1Q2b/F"]

variables += ["label_lin_0/I", "label_lin_1/I", "label_lin_2/I", "label_lin_3/I"]
variables += ["label_quad_0/I", "label_quad_1/I", "label_quad_2/I", "label_quad_3/I"]

variables += ["gen_cand_sum_pt/F"]

categories = [
    {'name':'e',   'func':lambda p:abs(p.pdgId())==11}, #electrons 
    {'name':'mu',  'func':lambda p:abs(p.pdgId())==13}, #muons 
    {'name':'ph',  'func':lambda p:p.pdgId()==22}, #photons 
    {'name':'chh', 'func':lambda p:abs(p.pdgId())>100 and p.charge()!=0 }, #charged hadrons 
    {'name':'neh', 'func':lambda p:abs(p.pdgId())>100 and p.charge()==0 }, # neutral hadrons
]

cand_vars           =  "pt/F,etarel/F,phirel/F,eta/F,phi/F,pdgId/I"
cand_varnames       =  varnames( cand_vars )
cand_ch_vars        =  "pt/F,etarel/F,phirel/F,eta/F,phi/F,pdgId/I,charge/I"
cand_ch_varnames    =  varnames( cand_ch_vars )
for cat in categories:
    if cat['name'] in ['ph', 'neh']:
        variables.append( VectorTreeVariable.fromString("%s[%s]"%(cat['name'], cand_vars), nMax=1000 ) )
    else:
        variables.append( VectorTreeVariable.fromString("%s[%s]"%(cat['name'], cand_ch_vars), nMax=1000 ) )

def fill_vector_collection( event, collection_name, collection_varnames, objects):
    setattr( event, "n"+collection_name, len(objects) )
    for i_obj, obj in enumerate(objects):
        for var in collection_varnames:
            getattr(event, collection_name+"_"+var)[i_obj] = obj[var]
def fill_vector( event, collection_name, collection_varnames, obj):
    for var in collection_varnames:
        try:
            setattr(event, collection_name+"_"+var, obj[var] )
        except TypeError as e:
            logger.error( "collection_name %s var %s obj[var] %r", collection_name, var,  obj[var] )
            raise e
        except KeyError as e:
            logger.error( "collection_name %s var %s obj[var] %r", collection_name, var,  obj[var] )
            raise e

logger.info( "Running over files: %s", ", ".join(sample.files ) )

if args.delphesEra is not None:
    if args.delphesEra == 'RunII':
        from TMB.Tools.DelphesReader          import DelphesReader
        delphesCard = 'delphes_card_CMS'
    elif args.delphesEra == 'RunIICentral':
        from TMB.Tools.DelphesReader          import DelphesReader
        delphesCard = 'delphes_card_CMS_Central'
    elif args.delphesEra == 'RunIInoDelphesIso':
        from TMB.Tools.DelphesReader          import DelphesReader
        delphesCard = 'delphes_card_CMS_noLepIso'
    elif args.delphesEra == 'RunIIPileUp':
        from TMB.Tools.DelphesReader          import DelphesReader
        delphesCard = 'delphes_card_CMS_PileUp'
    elif args.delphesEra == 'PhaseII':
        from TMB.Tools.DelphesReaderCMSHLLHC  import DelphesReader
        delphesCard = 'CMS_PhaseII/CMS_PhaseII_200PU_v03'
readers = []

# FWLite reader if this is an EDM file
fwliteReader = sample.fwliteReader( products = products )
readers.append( fwliteReader )

# Delphes reader if we read Delphes
if args.delphesEra is not None:
    delphes_file = os.path.join( output_directory, 'delphes', sample.name+'.root' )
    if      ( not os.path.exists( delphes_file )) or \
            ( os.path.exists( delphes_file ) and not checkRootFile( delphes_file, checkForObjects=["Delphes"])) or \
            args.overwrite in ['all']:
        logger.debug( "Reproducing delphes file %s", delphes_file)
        delphesProducer = DelphesProducer( card = delphesCard )
        delphesProducer.produce( sample.files, delphes_file)
    delphesReader = DelphesReader( Sample.fromFiles( delphes_file, delphes_file, treeName = "Delphes" ) ) # RootTools version
    readers.append( delphesReader )

def addTLorentzVector( p_dict ):
    ''' add a TLorentz 4D Vector for further calculations
    '''
    p_dict['vecP4'] = ROOT.TLorentzVector( p_dict['pt']*cos(p_dict['phi']), p_dict['pt']*sin(p_dict['phi']),  p_dict['pt']*sinh(p_dict['eta']), p_dict['pt']*cosh(p_dict['eta']) )

tmp_dir     = ROOT.gDirectory
#post_fix = '_%i'%args.job if args.nJobs > 1 else ''
output_filename =  os.path.join(output_directory, sample.name + '.root')

_logger.   add_fileHandler( output_filename.replace('.root', '.log'), args.logLevel )
_logger_rt.add_fileHandler( output_filename.replace('.root', '_rt.log'), args.logLevel )

if os.path.exists( output_filename ) and checkRootFile( output_filename, checkForObjects=["Events"]) and args.overwrite =='none' :
    logger.info( "File %s found. Quit.", output_filename )
    sys.exit(0)

output_file = ROOT.TFile( output_filename, 'recreate')
output_file.cd()
maker = TreeMaker(
    #sequence  = [ filler ],
    variables = [ (TreeVariable.fromString(x) if type(x)==str else x) for x in variables ] + extra_variables,
    treeName = "Events"
    )

tmp_dir.cd()

gRandom = ROOT.TRandom3()
def filler( event ):

    event.run, event.lumi, event.evt = fwliteReader.evt
    if fwliteReader.position % 100==0: logger.info("At event %i/%i", fwliteReader.position, fwliteReader.nEvents)

    if args.addReweights:
        event.nweight = weightInfo.nid
        lhe_weights = fwliteReader.products['lhe'].weights()
        weights      = []
        param_points = []
        for weight in lhe_weights:
            # Store nominal weight (First position!)
            weight_id = weight.id.rstrip('_nlo')
            if weight_id in ['rwgt_1','dummy']: 
                event.rw_nominal = weight.wgt
            #print "Hello weight", weight_id, ( weight_id.lower() in weightInfo_data_lower.keys()) 
            if not weight_id.lower() in weightInfo_data_lower.keys(): 
                continue
            pos = weightInfo_data_lower[weight_id]
            #print "pos", weight.wgt, event.weight_base[pos]
            event.weight_base[pos] = weight.wgt
            weights.append( weight.wgt )
            interpreted_weight = interpret_weight(weight_id.lower()) 
            #for var in weightInfo.variables:
            #    getattr( event, "rw_"+var )[pos] = interpreted_weight[var]
            # weight data for interpolation
            if not hyperPoly.initialized: param_points.append( tuple(interpreted_weight[var.lower()] for var in weightInfo.variables) )

        # get list of values of ref point in specific order
        ref_point_coordinates = [weightInfo.ref_point_coordinates[var] for var in weightInfo.variables]

        # Initialize with Reference Point

        if not hyperPoly.initialized: 
            #print "evt,run,lumi", event.run, event.lumi, event.evt
            #print "ref point", ref_point_coordinates, "param_points", param_points
            #for i_p, p in enumerate(param_points):
                #print "weight", i_p, weights[i_p], " ".join([ "%s=%3.2f"%( weightInfo.variables[i], p[i]) for i in range(len(p)) if p[i]!=0])
            hyperPoly.initialize( param_points, ref_point_coordinates )

        coeff = hyperPoly.get_parametrization( weights )
        # = HyperPoly(weight_data, args.interpolationOrder)
        event.np = hyperPoly.ndof
        event.chi2_ndof = hyperPoly.chi2_ndof(coeff, weights)
        #logger.debug( "chi2_ndof %f coeff %r", event.chi2_ndof, coeff )
        if event.chi2_ndof>10**-6: logger.warning( "chi2_ndof is large: %f", event.chi2_ndof )
        for n in xrange(hyperPoly.ndof):
            event.p_C[n] = coeff[n]
            if n>0:
                setattr( event, "target_"+target_coeff[n-1], coeff[n]/coeff[0] )
        # lumi weight / w0

    index_lin  = weightInfo.combinations.index(('ctWRe',))
    lin_rel     = event.p_C[index_lin]/event.p_C[0]
    event.label_lin_0 = lin_rel<=-0.3 
    event.label_lin_1 = lin_rel>-0.3 and lin_rel<=-0.2 
    event.label_lin_2 = lin_rel>-0.2 and lin_rel<=-0.1 
    event.label_lin_3 = lin_rel>-0.1  
    index_quad = weightInfo.combinations.index(('ctWRe','ctWRe'))
    quad_rel    = event.p_C[index_quad]/event.p_C[0]
    event.label_quad_0 = quad_rel<=0.01 
    event.label_quad_1 = quad_rel>0.01 and quad_rel<=0.02
    event.label_quad_2 = quad_rel>0.02 and quad_rel<=0.03
    event.label_quad_3 = quad_rel>0.03  

    # All gen particles
    gp        = fwliteReader.products['gp']
    #gpPacked  = fwliteReader.products['gpPacked']

    # for searching
    search  = GenSearch( gp )

    # all gen-tops
    gen_tops = filter( lambda p:abs(p.pdgId())==6 and search.isLast(p), gp)

    # jets
    ak8GenJets = fwliteReader.products['ak8GenJets']
    genJets = filter( genJetId, ak8GenJets )

    # hadronic gen_tops
    W, hadronic_gen_top = None, None
    for i_gen_top, gen_top in enumerate(gen_tops):
        t_daughters = [ gen_top.daughter(i) for i in range(gen_top.numberOfDaughters())]
        logger.debug( 'top',i_gen_top,"pdgId", gen_top.pdgId(),"d_pdgid", [d.pdgId() for d in t_daughters] )

        # we must always have a b
        b = next( (search.descend(b_) for b_ in t_daughters if abs(b_.pdgId()) == 5), None)
        W = next( (search.descend(W_) for W_ in t_daughters if abs(W_.pdgId()) == 24), None)

        # W resonance
        if W is not None and abs(W.daughter(0).pdgId()) in [1,2,3,4]:
            quark1 = W.daughter(0)
            quark2 = W.daughter(1)
            hadronic_gen_top = gen_top
            break
        # 3-body
        elif len(t_daughters)==3:
            quark1 = next( (q for q in t_daughters if abs(q.pdgId()) in [1,2,3,4]), None)
            quark2 = next( (q for q in t_daughters if abs(q.pdgId()) in [1,2,3,4] and not q.pdgId()==quark1.pdgId()), None)
            if quark1 is not None and quark2 is not None:
                hadronic_gen_top = gen_top
                break
        # leptonic
        else:
            quark1, quark2 = None, None

    if hadronic_gen_top is not None and b is not None:

        # quark2 is the anti-quark
        if quark1.pdgId()<0:
            quark1, quark2 = quark2, quark1

        event.partonTop_pt   = hadronic_gen_top.pt()
        event.partonTop_eta  = hadronic_gen_top.eta()
        event.partonTop_phi  = hadronic_gen_top.phi()
        event.partonTop_mass = hadronic_gen_top.mass()
        event.partonTop_pdgId= hadronic_gen_top.pdgId()

        event.genQ1_pt   = quark1.pt()
        event.genQ1_eta  = quark1.eta()
        event.genQ1_phi  = quark1.phi()
        event.genQ1_mass = quark1.mass()
        event.genQ1_pdgId= quark1.pdgId()

        event.genQ2_pt   = quark2.pt()
        event.genQ2_eta  = quark2.eta()
        event.genQ2_phi  = quark2.phi()
        event.genQ2_mass = quark2.mass()
        event.genQ2_pdgId= quark2.pdgId()

        event.genb_pt    = b.pt()
        event.genb_eta   = b.eta()
        event.genb_phi   = b.phi()
        event.genb_mass  = b.mass()
        event.genb_pdgId = b.pdgId()

        # 4-vectors
        vec_quark1, vec_quark2, vec_b, vec_W, vec_top = ROOT.TLorentzVector(), ROOT.TLorentzVector(), ROOT.TLorentzVector(), ROOT.TLorentzVector(), ROOT.TLorentzVector()

        vec_quark1.SetPtEtaPhiM(quark1.pt(),quark1.eta(),quark1.phi(),quark1.mass())
        vec_quark2.SetPtEtaPhiM(quark2.pt(),quark2.eta(),quark2.phi(),quark2.mass())
        vec_b.SetPtEtaPhiM(b.pt(),b.eta(),b.phi(),b.mass())
        vec_t = vec_quark1 + vec_quark2 + vec_b
        vec_W = vec_quark1 + vec_quark2

        event.genW_pt   = vec_W.Pt()
        event.genW_eta  = vec_W.Eta()
        event.genW_phi  = vec_W.Phi()
        event.genW_mass = vec_W.M()
        
        event.genTop_pt   = vec_t.Pt()
        event.genTop_eta  = vec_t.Eta()
        event.genTop_phi  = vec_t.Phi()
        event.genTop_mass = vec_t.M()

        # compute theta and phi
        beam = ROOT.TLorentzVector()
        beam.SetPxPyPzE(0,0,6500,6500)

        boost_t = vec_t.BoostVector()

        vec_W.Boost(-boost_t)
        vec_quark1.Boost(-boost_t)
        vec_quark2.Boost(-boost_t)

        n_scatter = ((beam.Vect().Unit()).Cross(vec_W.Vect())).Unit()
        n_decay   = (vec_quark1.Vect().Cross(vec_quark2.Vect())).Unit()
        sign_flip =  1 if ( ((n_scatter.Cross(n_decay))*(vec_W.Vect())) > 0 ) else -1

        try:
            event.gen_phi = sign_flip*acos(n_scatter.Dot(n_decay))
        except ValueError:
            event.gen_phi = -100

        boost_W = vec_W.BoostVector()

        vec_quark1.Boost(-boost_W)

        gen_theta = float('nan')
        try:
            event.gen_theta = (vec_W).Angle(vec_quark1.Vect())
        except ValueError:
            event.gen_theta = -100

        #print "theta,phi",gen_theta,gen_phi

        # match gen jet
        matched_genJet = next( (j for j in genJets if deltaR( {'eta':j.eta(),'phi':j.phi()}, {'eta':hadronic_gen_top.eta(), 'phi':hadronic_gen_top.phi()}) < 0.6), None)

        if matched_genJet:
            gen_particles = list(matched_genJet.getJetConstituentsQuick())

            event.genJet_pt  = matched_genJet.pt()
            event.genJet_eta = matched_genJet.eta()
            event.genJet_phi = matched_genJet.phi()
            event.genJet_mass= matched_genJet.mass()
            event.genJet_nConstituents  = len( gen_particles )
            event.genJet_isMuon         = matched_genJet.isMuon()
            event.genJet_isElectron     = matched_genJet.isElectron()
            event.genJet_isPhoton       = matched_genJet.isPhoton()

            event.gen_cand_sum_pt = sum([c.pt() for c in gen_particles],0)

            event.dR_genJet_Q1  = deltaR( {'phi':event.genQ1_phi, 'eta':event.genQ1_eta},   {'phi':event.genJet_phi, 'eta':event.genJet_eta} )
            event.dR_genJet_Q2  = deltaR( {'phi':event.genQ2_phi, 'eta':event.genQ2_eta},   {'phi':event.genJet_phi, 'eta':event.genJet_eta} )
            event.dR_genJet_W   = deltaR( {'phi':event.genW_phi, 'eta':event.genW_eta},     {'phi':event.genJet_phi, 'eta':event.genJet_eta} )
            event.dR_genJet_b   = deltaR( {'phi':event.genb_phi, 'eta':event.genb_eta},     {'phi':event.genJet_phi, 'eta':event.genJet_eta} )
            event.dR_genJet_top = deltaR( {'phi':event.genTop_phi, 'eta':event.genTop_eta}, {'phi':event.genJet_phi, 'eta':event.genJet_eta} )

            event.dR_genJet_maxQ1Q2b = max( [ event.dR_genJet_Q1, event.dR_genJet_Q2, event.dR_genJet_b ] )

            count = 0 
            for cat in categories:
                cands = filter( cat['func'], gen_particles )
                cands_list = [ {'pt':c.pt(), 'eta':c.eta(), 'phi':c.phi(), 'charge':c.charge(), 'pdgId':c.pdgId()} for c in cands ]
                for p in cands_list:
                    p['phirel'] = deltaPhi(event.genJet_phi, p['phi'], returnAbs=False)
                    #p['phirel'] = acos(cos(p['phi'] - event.genJet_phi) )
                    p['etarel'] = p['eta'] - event.genJet_eta
                    #print (p['pdgId'], p['charge'])

                cands_list.sort( key = lambda p:-p['pt'] )
                addIndex( cands_list )
                if cat['name'] in ['ph', 'neh']:
                    fill_vector_collection( event, cat['name'], cand_varnames, cands_list)
                else:
                    fill_vector_collection( event, cat['name'], cand_ch_varnames, cands_list)
                count+=len(cands)
            assert count == len( gen_particles ), "Missing a gen particle in categorization!!"

            #clustSeq      = fastjet.ClusterSequence( map( make_pseudoJet, matched_genJet.getJetConstituentsQuick() ), ak8 )
            #sortedJets    = fastjet.sorted_by_pt(clustSeq.inclusive_jets())

            # softdrop
            genCandsVec = ROOT.vector("TLorentzVector")()
            for p in matched_genJet.getJetConstituentsQuick() :
                genCandsVec.push_back( ROOT.TLorentzVector( p.p4().Px(), p.p4().Py(), p.p4().Pz(), p.p4().E()) )
            genSDJets = softDrop.result( genCandsVec )
            if genSDJets.size()>=1:
                genSDJet = genSDJets[0]
                event.genJet_SDmass = genSDJet.m() # softdrop mass

                genSDSubJets = genSDJet.pieces()
                if len(genSDSubJets)>0:
                    event.genJet_SDsubjet0_eta      = genSDSubJets[0].eta()
                    event.genJet_SDsubjet0_deltaEta = genSDSubJets[0].eta() - matched_genJet.eta()
                    event.genJet_SDsubjet0_phi      = genSDSubJets[0].phi()
                    event.genJet_SDsubjet0_deltaPhi = deltaPhi(genSDSubJets[0].phi(), matched_genJet.phi(), returnAbs=False)
                    event.genJet_SDsubjet0_deltaR   = sqrt( event.genJet_SDsubjet0_deltaPhi**2 + event.genJet_SDsubjet0_deltaEta**2 ) 
                    event.genJet_SDsubjet0_mass     = genSDSubJets[0].m()
                if len(genSDSubJets)>1:
                    event.genJet_SDsubjet1_eta      = genSDSubJets[1].eta()
                    event.genJet_SDsubjet1_deltaEta = genSDSubJets[1].eta() - matched_genJet.eta()
                    event.genJet_SDsubjet1_phi      = genSDSubJets[1].phi()
                    event.genJet_SDsubjet1_deltaPhi = deltaPhi(genSDSubJets[1].phi(), matched_genJet.phi(), returnAbs=False)
                    event.genJet_SDsubjet1_deltaR   = sqrt( event.genJet_SDsubjet1_deltaPhi**2 + event.genJet_SDsubjet1_deltaEta**2 ) 
                    event.genJet_SDsubjet1_mass     = genSDSubJets[1].m()
            
            ns_tau = nSubjettiness.getTau( 4, genCandsVec )
            event.genJet_tau1 = ns_tau[0]
            event.genJet_tau2 = ns_tau[1]
            event.genJet_tau3 = ns_tau[2]
            event.genJet_tau4 = ns_tau[3]
            event.genJet_tau21 = ns_tau[1]/ns_tau[0] if ns_tau[0]>0 else 0
            event.genJet_tau32 = ns_tau[2]/ns_tau[1] if ns_tau[1]>0 else 0

            ecf.setParticles( genCandsVec )
            result = ecf.result()
            for i_ecf, (name, _) in enumerate( ecfs ):
                setattr( event, "genJet_%s"%name, result[i_ecf] )

            # only fill if we have everything. This is a per-top ntuple, not per-events
            maker.fill()

    maker.event.init()

counter = 0
for reader in readers:
    reader.start()
maker.start()

while readers[0].run( ):
    for reader in readers[1:]:
        reader.run()

    filler( maker.event )
         
    counter += 1
    if counter == maxEvents:  break

logger.info( "Done with running over %i events.", readers[0].nEvents )

output_file.cd()
maker.tree.Write()
output_file.Close()

logger.info( "Written output file %s", output_filename )

##cleanup delphes file:
if os.path.exists( output_filename ) and args.delphesEra is not None and args.removeDelphesFiles:
    os.remove( delphes_file )
    logger.info( "Removing Delphes file %s", delphes_file )
