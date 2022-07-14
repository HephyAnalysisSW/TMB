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
import  TMB.Tools.VV_angles           as VV_angles

#from TMB.Tools.UpgradeJECUncertainty  import UpgradeJECUncertainty 
#
# Arguments
# 
import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--logLevel',           action='store',      default='INFO',          nargs='?', choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'], help="Log level for logging")
argParser.add_argument('--small',              action='store_true', help='Run only on a small subset of the data?')#, default = True)
argParser.add_argument('--miniAOD',            action='store_true', help='Run on miniAOD?')#, default = True)
argParser.add_argument('--overwrite',          action='store',      nargs='?', choices = ['none', 'all', 'target'], default = 'none', help='Overwrite?')#, default = True)
argParser.add_argument('--targetDir',          action='store',      default='v7')
argParser.add_argument('--sample',             action='store',      default='fwlite_ttZ_ll_LO_scan', help="Name of the sample loaded from fwlite_benchmarks. Only if no inputFiles are specified")
argParser.add_argument('--inputFiles',         action='store',      nargs = '*', default=[])
argParser.add_argument('--delphesEra',         action='store',      default = None, choices = ["RunII", "RunIICentral", "RunIInoDelphesIso", "RunIIPileUp", "PhaseII"], help="specify delphes era")
argParser.add_argument('--targetSampleName',   action='store',      default=None, help="Name of the sample in case inputFile are specified. Otherwise ignored")
argParser.add_argument('--nJobs',              action='store',      nargs='?', type=int, default=1,  help="Maximum number of simultaneous jobs.")
argParser.add_argument('--job',                action='store',      nargs='?', type=int, default=0,  help="Run only job i")
argParser.add_argument('--addReweights',       action='store_true',   help="Add reweights?")
argParser.add_argument('--combinatoricalBTags', action='store_true',   help="BTags combinatorical?")
argParser.add_argument('--removeDelphesFiles', action='store_true',   help="remove Delphes file after postprocessing?")
argParser.add_argument('--interpolationOrder', action='store',      nargs='?', type=int, default=3,  help="Interpolation order for EFT weights.")
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
    if args.miniAOD:
        sample_file = "$CMSSW_BASE/python/TMB/Samples/miniAOD.py"
    elif args.targetDir.startswith('v4'):
        sample_file = "$CMSSW_BASE/python/TMB/Samples/gen_v4.py"
    elif args.targetDir.startswith('v5'):
        sample_file = "$CMSSW_BASE/python/TMB/Samples/gen_v5.py"
    elif args.targetDir.startswith('v6'):
        sample_file = "$CMSSW_BASE/python/TMB/Samples/gen_v6.py"
    elif args.targetDir.startswith('v7'):
        sample_file = "$CMSSW_BASE/python/TMB/Samples/gen_v7.py"
    elif args.targetDir.startswith('v8'):
        sample_file = "$CMSSW_BASE/python/TMB/Samples/gen_v8.py"
    elif args.targetDir.startswith('v9'):
        sample_file = "$CMSSW_BASE/python/TMB/Samples/gen_v9.py"
    elif args.targetDir.startswith('v10'):
        sample_file = "$CMSSW_BASE/python/TMB/Samples/gen_v10.py"
    samples = imp.load_source( "samples", os.path.expandvars( sample_file ) )
    sample = getattr( samples, args.sample )
    logger.debug( 'Loaded sample %s with %i files.', sample.name, len(sample.files) )

maxEvents = -1
if args.small: 
    args.targetDir += "_small"
    maxEvents       = 500 
    sample.files=sample.files[:1]

xsec = sample.xsec if hasattr( sample, "xsec" ) else sample.xSection 
nEvents = sample.nEvents
lumiweight1fb = xsec * 1000. / nEvents

# output directory
output_directory = os.path.join(skim_output_directory, 'gen', args.targetDir, sample.name) 

if not os.path.exists( output_directory ): 
    try:
        os.makedirs( output_directory )
    except OSError:
        pass
    logger.info( "Created output directory %s", output_directory )

# Load reweight pickle file if supposed to keep weights. 
extra_variables = []
if args.addReweights:
    weightInfo = WeightInfo( sample.reweight_pkl )
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


# Run only job number "args.job" from total of "args.nJobs"
if args.nJobs>1:
    n_files_before = len(sample.files)
    sample = sample.split(args.nJobs)[args.job]
    n_files_after  = len(sample.files)
    logger.info( "Running job %i/%i over %i files from a total of %i.", args.job, args.nJobs, n_files_after, n_files_before)

max_jet_abseta = 5.1

if args.miniAOD:
    products = {
        'lhe':{'type':'LHEEventProduct', 'label':("externalLHEProducer")},
        'gen':{'type':'GenEventInfoProduct', 'label':'generator'},
        'gp':{'type':'vector<reco::GenParticle>', 'label':("prunedGenParticles")},
        #'gpPacked':{'type':'vector<pat::PackedGenParticle>', 'label':("packedGenParticles")},
        #'gp':{'type':'vector<pat::PackedGenParticle>', 'label':("packedGenParticles")},
        'genJets':{'type':'vector<reco::GenJet>', 'label':("slimmedGenJets")},
        'genMET':{'type':'vector<reco::GenMET>',  'label':("genMetTrue")},
    }
else:
    products = {
        'lhe':{'type':'LHEEventProduct', 'label':("externalLHEProducer")},
        'gen':{'type':'GenEventInfoProduct', 'label':'generator'},
        'gp':{'type':'vector<reco::GenParticle>', 'label':("genParticles")},
        'genJets':{'type':'vector<reco::GenJet>', 'label':("ak4GenJets")},
        'genMET':{'type':'vector<reco::GenMET>',  'label':("genMetTrue")},
    }

def varnames( vec_vars ):
    return [v.split('/')[0] for v in vec_vars.split(',')]

def vecSumPt(*args):
    return sqrt( sum([o['pt']*cos(o['phi']) for o in args],0.)**2 + sum([o['pt']*sin(o['phi']) for o in args],0.)**2 )

def addIndex( collection ):
    for i  in range(len(collection)):
        collection[i]['index'] = i

# upgrade JEC
# upgradeJECUncertainty = UpgradeJECUncertainty()

# variables
variables = []

# Lumi weight 1fb
variables += ["lumiweight1fb/F"]

# EDM standard variables
variables  += ["run/I", "lumi/I", "evt/l"]

# MET
variables += ["genMet_pt/F", "genMet_phi/F"]
# jet vector
jet_read_vars       =  "pt/F,eta/F,phi/F,isMuon/I,isElectron/I,isPhoton/I"
jet_read_varnames   =  varnames( jet_read_vars )
jet_write_vars      = jet_read_vars+',matchBParton/I,matchCParton/I' 
jet_write_varnames  =  varnames( jet_write_vars )
variables += ["genJet[%s]"%jet_write_vars]
variables += ["genBj0_%s"%var for var in jet_write_vars.split(',')]
variables += ["genBj1_%s"%var for var in jet_write_vars.split(',')]
#VBS (gen)
variables += ["genVBSj0_%s"%var for var in jet_write_vars.split(',')]
variables += ["genVBSj1_%s"%var for var in jet_write_vars.split(',')]
variables += ["genVBSj0_index/I", "genVBSj1_index/I","genVBS_dEta/F", "genVBS_dPhi/F", "genVBS_mjj/F"]
# lepton vector 
lep_vars       =  "pt/F,eta/F,phi/F,pdgId/I,status/I"
lep_extra_vars =  "mother_pdgId/I,grandmother_pdgId/I"
lep_varnames   =  varnames( lep_vars ) 
lep_all_varnames = lep_varnames + varnames(lep_extra_vars)
variables     += ["genLep[%s]"%(','.join([lep_vars, lep_extra_vars]))]
## associated jet indices
#variables += [ "genBjLeadlep_index/I", "genBjLeadhad_index/I" ]
#variables += [ "genBjNonZlep_index/I", "genBjNonZhad_index/I" ]
# top vector
top_vars       =  "pt/F,eta/F,phi/F,pdgId/I,mass/F"
top_varnames   =  varnames( top_vars ) 
variables     += ["genTop[%s]"%top_vars]

if args.combinatoricalBTags:
    variables += ["combinatoricalBTagWeight2b/F"]

# to be stored for each boson
boson_varnames = [ 'pt', 'phi', 'eta', 'mass', 'status']
# Z vector from gen collection
boson_all_varnames = boson_varnames + ['cosThetaStar', 'daughter_pdgId','l1_index', 'l2_index', 'mother_pdgId', 'grandmother_pdgId']
variables     += ["genZ[pt/F,phi/F,eta/F,mass/F,status/I,cosThetaStar/F,daughter_pdgId/I,mother_pdgId/I,grandmother_pdgId/I,l1_index/I,l2_index/I]"]
variables     += ["genW[pt/F,phi/F,eta/F,mass/F,status/I,cosThetaStar/F,daughter_pdgId/I,mother_pdgId/I,grandmother_pdgId/I,l1_index/I,l2_index/I]"]
variables     += ["LHE_genZ_pt/F", "LHE_genZ_eta/F", "LHE_genZ_phi/F", "LHE_genZ_mass/F"]
variables     += ["LHE_genH_pt/F", "LHE_genH_eta/F", "LHE_genH_phi/F", "LHE_genH_mass/F"]
variables     += ["x1/F", "x2/F"]
# Z vector from genleps
# gamma vector
gen_photon_vars = "pt/F,phi/F,eta/F,mass/F,mother_pdgId/I,grandmother_pdgId/I,isISR/I,relIso04/F"#,minLeptonDR/F,minJetDR/F"
variables     += ["genPhoton[%s]"%gen_photon_vars]
gen_photon_varnames = varnames( gen_photon_vars )

if args.delphesEra is not None:
    # reconstructed bosons
    variables     += ["recoZ_l1_index/I", "recoZ_l2_index/I", "recoNonZ_l1_index/I", "recoNonZ_l2_index/I",  "recoZ_pt/F", "recoZ_eta/F", "recoZ_phi/F", "recoZ_mass/F", "recoZ_lldPhi/F", "recoZ_lldR/F", "recoZ_cosThetaStar/F"]

    # reconstructed leptons
    recoLep_vars       = "pt/F,eta/F,phi/F,pdgId/I,isolationVar/F,isolationVarRhoCorr/F,sumPtCharged/F,sumPtNeutral/F,sumPtChargedPU/F,sumPt/F,ehadOverEem/F,genMatched/I"
    variables         += ["recoLep[%s]"%recoLep_vars]
    recoLep_varnames  = varnames( recoLep_vars )
    # generated jets 
    variables += ["ndelphesGenJet/I", "delphesGenJet[pt/F,eta/F,phi/F]"]
    # reconstructed jets
    recoJet_vars    = 'pt/F,eta/F,phi/F,bTag/F,bTagPhys/I,nCharged/I,nNeutrals/I,matchGenBJet/I'#,pt_JEC_up/F,pt_JEC_up/F'

    btagWPs = ["loose"]#, "medium", "tight"] #, "looswMTD", "mediumMTD", "tightMTD"]
    default_btagWP = "loose"
    variables.append( "nBTag/I" )
    #variables.append( "nrecoJets_JEC_down/I" )
    #variables.append( "nrecoJets_JEC_up/I" )
    #variables.append( "nBTag_JEC_down/I" )
    #variables.append( "nBTag_JEC_up/I" )
    for btagWP in btagWPs:
        variables.append( "nBTag_"+btagWP+"/I" )
        recoJet_vars += ',bTag_'+btagWP+"/I"

    variables += ["recoJet[%s]"%recoJet_vars]
    recoJet_varnames = varnames( recoJet_vars )
    variables += ["recoBj0_%s"%var for var in recoJet_vars.split(',')]
    variables += ["recoBj1_%s"%var for var in recoJet_vars.split(',')]

    # associated jet indices
    variables += [ "recoBjNonZlep_index/I", "recoBjNonZhad_index/I" ]
    variables += [ "recoBjLeadlep_index/I", "recoBjLeadhad_index/I" ]

    #VBS (gen)
    variables += ["recoVBSj0_%s"%var for var in recoJet_vars.split(',')]
    variables += ["recoVBSj1_%s"%var for var in recoJet_vars.split(',')]
    variables += ["recoVBSj0_index/I", "recoVBSj1_index/I", "recoVBS_dEta/F", "recoVBS_dPhi/F", "recoVBS_mjj/F"]

    # reconstructed photons
    recoPhoton_vars = 'pt/F,eta/F,phi/F,isolationVar/F,isolationVarRhoCorr/F,sumPtCharged/F,sumPtNeutral/F,sumPtChargedPU/F,sumPt/F,ehadOverEem/F,minLeptonDR/F,minLeptonPt/F,minJetDR/F'#genIndex/I
    variables      += ["recoPhoton[%s]"%recoPhoton_vars]
    recoPhoton_varnames = varnames( recoPhoton_vars )

    variables      += ["recoMet_pt/F", "recoMet_phi/F"]
    variables      += ["delphesGenMet_pt/F", "delphesGenMet_phi/F"]

    # Systematics
    from TMB.Tools.bTagEff.delphesBTaggingEff import getBTagSF_1a
    #variables      += ["reweight_BTag_B/F", "reweight_BTag_L/F"]
    variables      += ["reweight_id_mu/F", "reweight_id_ele/F"]

    # VH variables
    variables += [
        "H_dijet_mass/F", "H_pt/F", "H_j1_index/I", "H_j2_index/I",
        "WH_W_pt/F", "WH_dPhiMetLep/F", "WH_MT/F", "WH_nu_pt/F", "WH_nu_eta/F", "WH_nu_phi/F", "WH_nu_E/F", "WH_Theta/F", "WH_theta/F", "WH_phi/F",
        "ZH_l1_index/I", "ZH_l2_index/I", "ZH_Theta/F", "ZH_theta/F", "ZH_phi/F",
    ]

if args.addReweights:
    variables.append('rw_nominal/F')
    # Lumi weight 1fb / w_0
    variables.append("ref_lumiweight1fb/F")


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
    variables = [ TreeVariable.fromString(x) for x in variables ] + extra_variables,
    treeName = "Events"
    )

tmp_dir.cd()

gRandom = ROOT.TRandom3()
def filler( event ):

    event.lumiweight1fb = lumiweight1fb

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
            print "evt,run,lumi", event.run, event.lumi, event.evt
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
            #print n, "p_C", n, coeff[n]
        # lumi weight / w0
        event.ref_lumiweight1fb = event.lumiweight1fb / coeff[0]

    # parton x1 x2
    event.x1 = fwliteReader.products['gen'].pdf().x.first
    event.x2 = fwliteReader.products['gen'].pdf().x.second

    # MET
    genMet = {'pt':fwliteReader.products['genMET'][0].pt(), 'phi':fwliteReader.products['genMET'][0].phi()}
    event.genMet_pt  = genMet['pt']
    event.genMet_phi = genMet['phi'] 

    # All gen particles
    gp        = fwliteReader.products['gp']
    #gpPacked  = fwliteReader.products['gpPacked']

    # for searching
    search  = GenSearch( gp )

    # find heavy objects before they decay
    genTops = map( lambda t:{var: getattr(t, var)() for var in top_varnames}, filter( lambda p:abs(p.pdgId())==6 and search.isLast(p),  gp) )

    genTops.sort( key = lambda p:-p['pt'] )
    fill_vector_collection( event, "genTop", top_varnames, genTops ) 

    # generated leptons from SM bosons
    genLeps    = [ (search.ascend(l), l) for l in filter( lambda p:abs(p.pdgId()) in [11, 12, 13, 14, 15, 16]  and abs(p.mother(0).pdgId()) in [22, 23, 24, 25], gp)]
    genLeps.sort( key = lambda p: -p[1].pt() )
    genLeps_from_bosons =  [first for first, last in genLeps]
    genLeps_dict = [ {var: getattr(last, var)() for var in lep_varnames} for first, last in genLeps ]
    addIndex( genLeps_dict )
    for i_genLep, (first, last) in enumerate(genLeps):
        mother = first.mother(0) if first.numberOfMothers()>0 else None
        if mother is not None:
            mother_pdgId      = mother.pdgId()
            mother_ascend     = search.ascend(mother)
            grandmother       = mother_ascend.mother(0) if mother_ascend.numberOfMothers()>0 else None
            grandmother_pdgId = grandmother.pdgId() if grandmother is not None else 0
        else:
            mother_pdgId = 0
            grandmother_pdgId = 0 
        genLeps_dict[i_genLep]['mother_pdgId']      = mother_pdgId
        genLeps_dict[i_genLep]['grandmother_pdgId'] = grandmother_pdgId
    fill_vector_collection( event, "genLep", lep_all_varnames, genLeps_dict ) 

    # LHE Zs (Suman's example: https://github.com/Sumantifr/XtoYH/blob/master/Analysis/NTuplizer/plugins/NTuplizer_XYH.cc#L1973-L1992 )
    hepup = fwliteReader.products['lhe'].hepeup()
    lhe_particles = hepup.PUP
    for i_p, p in enumerate(lhe_particles):
        if hepup.IDUP[i_p]!=23: continue
        p4 = ROOT.TLorentzVector(lhe_particles[i_p][0], lhe_particles[i_p][1], lhe_particles[i_p][2], lhe_particles[i_p][3])
        event.LHE_genZ_pt  = p4.Pt()
        event.LHE_genZ_eta = p4.Eta()
        event.LHE_genZ_phi  = p4.Phi()
        event.LHE_genZ_mass = p4.M()
        break
    for i_p, p in enumerate(lhe_particles):
        if hepup.IDUP[i_p]!=25: continue
        p4 = ROOT.TLorentzVector(lhe_particles[i_p][0], lhe_particles[i_p][1], lhe_particles[i_p][2], lhe_particles[i_p][3])
        event.LHE_genH_pt  = p4.Pt()
        event.LHE_genH_eta = p4.Eta()
        event.LHE_genH_phi  = p4.Phi()
        event.LHE_genH_mass = p4.M()
        break
        #print "LHE Z", i_p, hepup.IDUP[i_p], p4.Pt(), "status", hepup.ISTUP[i_p]

    # generated Zs decaying to leptons
    genZs = [ (search.ascend(genZ), genZ) for genZ in filter( lambda p:abs(p.pdgId())==23 and search.isLast(p) and abs(p.daughter(0).pdgId()) in [11, 13, 15], gp) ]
    genZs.sort( key = lambda p: -p[1].pt() )
    genZs_dict = [ {var: getattr(last, var)() for var in boson_varnames} for first, last in genZs ]
    for i_genZ, (first, last) in enumerate(genZs):

        #if first.numberOfMothers()>1:
        #    print "Two mothers", first.mother(0).pdgId(), search.ascend(first.mother(0)).pdgId(), first.mother(1).pdgId(), search.ascend(first.mother(1)).pdgId(), search.ascend(first.mother(0)).eta(), search.ascend(first.mother(1)).eta()

        mother = first.mother(0) if first.numberOfMothers()>0 else None
        if mother is not None:
            mother_pdgId      = mother.pdgId()
            mother_ascend     = search.ascend(mother)
            grandmother       = mother_ascend.mother(0) if mother_ascend.numberOfMothers()>0 else None
            grandmother_pdgId = grandmother.pdgId() if grandmother is not None else 0
        else:
            mother_pdgId = 0
            grandmother_pdgId = 0 
        genZs_dict[i_genZ]['mother_pdgId']      = mother_pdgId
        genZs_dict[i_genZ]['grandmother_pdgId'] = grandmother_pdgId

        d1, d2 = last.daughter(0), last.daughter(1)
        if d1.pdgId()>0: 
            lm, lp = d1, d2
        else:
            lm, lp = d2, d1
        genZs_dict[i_genZ]['daughter_pdgId'] = lm.pdgId()
        genZs_dict[i_genZ]['l1_index'] = genLeps_from_bosons.index( last.daughter(0) )
        genZs_dict[i_genZ]['l2_index'] = genLeps_from_bosons.index( last.daughter(1) )
        genZs_dict[i_genZ]['cosThetaStar'] = cosThetaStar(last.mass(), last.pt(), last.eta(), last.phi(), lm.pt(), lm.eta(), lm.phi())
        
        #if genZs_dict[i_genZ]['mother_pdgId']==21:
        #    print genZs_dict[i_genZ]['eta']

    fill_vector_collection( event, "genZ", boson_all_varnames, genZs_dict ) 

    # generated W decaying to leptons
    genWs = [ (search.ascend(genW), genW) for genW in filter( lambda p:abs(p.pdgId())==24 and search.isLast(p) and abs(p.daughter(0).pdgId()) in [11, 12, 13, 14, 15, 16], gp)]
    genWs.sort( key = lambda p: -p[1].pt() )
    genWs_dict = [ {var: getattr(last, var)() for var in boson_varnames} for first, last in genWs ]
    for i_genW, (first, last) in enumerate(genWs):

        mother = first.mother(0) if first.numberOfMothers()>0 else None
        if mother is not None:
            mother_pdgId      = mother.pdgId()
            mother_ascend     = search.ascend(mother)
            grandmother       = mother_ascend.mother(0) if mother_ascend.numberOfMothers()>0 else None
            grandmother_pdgId = grandmother.pdgId() if grandmother is not None else 0
        else:
            mother_pdgId = 0
            grandmother_pdgId = 0 
        genWs_dict[i_genW]['mother_pdgId']      = mother_pdgId
        genWs_dict[i_genW]['grandmother_pdgId'] = grandmother_pdgId

        d1, d2 = last.daughter(0), last.daughter(1)
        if abs(d1.pdgId()) in [11, 13, 15]: 
            l, nu = d1, d2
        else:
            nu, l = d2, d1
        genWs_dict[i_genW]['daughter_pdgId'] = l.pdgId()
        genWs_dict[i_genW]['l1_index'] = genLeps_from_bosons.index( last.daughter(0) )
        genWs_dict[i_genW]['l2_index'] = genLeps_from_bosons.index( last.daughter(1) )
        genWs_dict[i_genW]['cosThetaStar'] = cosThetaStar(last.mass(), last.pt(), last.eta(), last.phi(), l.pt(), l.eta(), l.phi())
    fill_vector_collection( event, "genW", boson_all_varnames, genWs_dict ) 
    
#
#
### this def should only be called in debug mode
##        def printTopMothers():
##                genTopCheck = [ (search.ascend(l), l, search.ancestry( search.ascend(l) )) for l in filter( lambda p:abs(p.pdgId())==6,  gp) ]
##
##                print genTopCheck
##                genTopCheck.sort( key = lambda p: -p[1].pt() )
##                if len(genTopCheck) > 0:
##                    topPdg_ids = filter( lambda p:p!=2212, [abs(particle.pdgId()) for particle in genTopCheck[0][2]])
##                    print 'TOP 1'
##                    print topPdg_ids
##                    previous = genTopCheck[0][0].mother(1)
##                    print 'first', genTopCheck[0][0].pdgId()
##                    for i, pdg in enumerate(topPdg_ids):
##                        try:
##                            print 'mother', i, previous.pdgId()
##                            print 'mothers', i, [p.pdgId() for p in search.ancestry( previous )]
##                            previous = previous.mother(0)
##                        except Exception as val:
##                            print val
##                            break
##                if len(genTopCheck) > 1:
##                    topPdg_ids = filter( lambda p:p!=2212, [abs(particle.pdgId()) for particle in genTopCheck[1][2]])
##                    print 'TOP 2'
##                    print topPdg_ids
##                    previous = genTopCheck[1][0].mother(1)
##                    print 'first', genTopCheck[1][0].pdgId()
##                    for i, pdg in enumerate(topPdg_ids):
##                        try:
##                            print 'mother', i, previous.pdgId()
##                            print 'mothers', i, [p.pdgId() for p in search.ancestry( previous )]
##                            previous = previous.mother(0)
##                        except Exception as val:
##                            print val
##                            break
#
#
 # Gen photons: particle-level isolated gen photons
    genPhotons = [ ( search.ascend(l), l ) for l in filter( lambda p: abs( p.pdgId() ) == 22 and p.pt() > 20 and search.isLast(p) and p.status()==1, gp ) ]
    genPhotons.sort( key = lambda p: -p[1].pt() )
    genPhotons_dict = [ {var: getattr(genPhoton[1], var)() for var in boson_varnames} for genPhoton in genPhotons ]

    for i_genPhoton, (first, last) in enumerate(genPhotons):
        mother_pdgId = first.mother(0).pdgId() if first.numberOfMothers() > 0 else -999

        genPhotons_dict[i_genPhoton]['motherPdgId'] = mother_pdgId
        genPhotons_dict[i_genPhoton]['status']      = last.status()

        mother_ascend     = search.ascend( first.mother(0) )
        grandmother       = mother_ascend.mother(0) if first.mother(0).numberOfMothers() > 0 else None
        grandmother_pdgId = grandmother.pdgId() if grandmother else 0

        genPhotons_dict[i_genPhoton]['mother_pdgId']      = mother_pdgId
        genPhotons_dict[i_genPhoton]['grandmother_pdgId'] = grandmother_pdgId

        if abs(mother_pdgId) in [1,2,3,4,5,21,2212] and abs(grandmother_pdgId) in [1,2,3,4,5,21,2212]:
            genPhotons_dict[i_genPhoton]["isISR"] = 1 #also photons from gluons, as MG doesn't give you the right pdgId
        else:
            genPhotons_dict[i_genPhoton]["isISR"] = 0

        close_particles = filter( lambda p: p!=last and deltaR2( {'phi':last.phi(), 'eta':last.eta()}, {'phi':p.phi(), 'eta':p.eta()} ) < 0.16 , search.final_state_particles_no_neutrinos )
        genPhotons_dict[i_genPhoton]['relIso04'] = sum( [ p.pt() for p in close_particles ], 0 ) / last.pt()
        #GenPhoton['photonJetdR'] =  999
        #GenPhoton['photonLepdR'] =  999

    # store isolated photons
    fill_vector_collection( event, "genPhoton", gen_photon_varnames, filter( lambda g: g['relIso04']<0.3, genPhotons_dict) ) 

    # jets
    fwlite_genJets = fwliteReader.products['genJets']
    #ngenJets_preFilter = len(fwlite_genJets)
    #fwlite_genJets = filter( genJetId, fwlite_genJets )
    #ngenJets_postFilter = len(fwlite_genJets)

    # make dict
    allGenJets = map( lambda t:{var: getattr(t, var)() for var in jet_read_varnames}, filter( lambda j:j.pt()>30, fwlite_genJets) )

    #print allGenJets
    #print genLeps_dict
    ngenJets_preFilter = len(allGenJets)
    allGenJets, _ = cleanJetsAndLeptons( allGenJets, genLeps_dict )
    ngenJets_postFilter = len(allGenJets)
    #print "filtered:", ngenJets_preFilter-ngenJets_postFilter

    # filter genJets
    genJets = list( filter( lambda j:isGoodGenJet( j, max_jet_abseta = max_jet_abseta), allGenJets ) )

    ## cleaning of jets with isolated photons
    #genJets = list( filter( lambda j:min([999]+[deltaR2(j, p) for p in genPhotons_ ])>0.4**2, genJets))

    ## store minimum DR to jets
    #for genPhoton in genPhotons_:
    #    genPhoton['minJetDR'] =  min([999]+[deltaR(genPhoton, j) for j in genJets])

    # find b's from tops:
    #b_partons = [ b for b in filter( lambda p:abs(p.pdgId())==5 and p.numberOfMothers()==1 and abs(p.mother(0).pdgId())==6,  gp) ]
    b_partons = [ b for b in filter( lambda p:abs(p.pdgId())==5 and search.isLast(p),  gp) ]
    c_partons = [ c for c in filter( lambda p:abs(p.pdgId())==4 and search.isLast(p),  gp) ]
    
    #for b in b_partons:
    #    print  b.pt(), b.eta(), b.phi(), b.pdgId(), b.numberOfMothers(), b.mother(0).pdgId()
    #print

    # store if gen-jet is DR matched to a B parton
    for genJet in genJets:
        genJet['matchBParton'] = ( min([999]+[deltaR2(genJet, {'eta':b.eta(), 'phi':b.phi()}) for b in b_partons]) < 0.2**2 )
        if not genJet['matchBParton']:
            genJet['matchCParton'] = ( min([999]+[deltaR2(genJet, {'eta':c.eta(), 'phi':c.phi()}) for c in c_partons]) < 0.2**2 )
        else:
            genJet['matchCParton'] = False

    genJets = filter( lambda j: (min([999]+[deltaR2(j, l) for l in genLeps_dict if l['pt']>10 and l['status']==1 and abs(l['pdgId']) in [11,13]]) > 0.3**2 ), genJets )
    genJets.sort( key = lambda p:-p['pt'] )
    addIndex( genJets )
    fill_vector_collection( event, "genJet", jet_write_varnames, genJets)

    # gen b jets
    trueBjets = list( filter( lambda j: j['matchBParton'], genJets ) )
    trueNonBjets = list( filter( lambda j: not j['matchBParton'], genJets ) )

    # Mimick b reconstruction ( if the trailing b fails acceptance, we supplement with the leading non-b jet ) 
    genBj0, genBj1 = ( trueBjets + trueNonBjets + [None, None] )[:2]
    if genBj0: fill_vector( event, "genBj0", jet_write_varnames, genBj0) 
    if genBj1: fill_vector( event, "genBj1", jet_write_varnames, genBj1) 

    # VBS jets
    combs = list(itertools.combinations(range(len(genJets)),2))
    combs.sort(key=lambda comb:-genJets[comb[0]]['pt']-genJets[comb[1]]['pt'])
    for i_j1, i_j2 in combs:
        if genJets[i_j1]['eta']*genJets[i_j2]['eta']<0:
            vbsGenJet0 = genJets[i_j1]
            vbsGenJet1 = genJets[i_j2]
            event.genVBSj0_index, event.genVBSj1_index = i_j1, i_j2
            fill_vector( event, "genVBSj0", jet_write_varnames, vbsGenJet0)
            fill_vector( event, "genVBSj1", jet_write_varnames, vbsGenJet1)
            event.genVBS_dEta= vbsGenJet0['eta']-vbsGenJet1['eta']
            event.genVBS_dPhi= deltaPhi(vbsGenJet0['phi'],vbsGenJet1['phi'])
            event.genVBS_mjj = sqrt( 2*vbsGenJet0['pt']*vbsGenJet1['pt']*(cosh(event.genVBS_dEta)-cos(event.genVBS_dPhi)))
            break

    # Reco quantities
    if args.delphesEra is not None:
        #delphesReader.event.GetEntry(fwliteReader.position-1 ) # do this differently now ...  

        delphes_genJets = delphesReader.genJets()
        fill_vector_collection( event, "delphesGenJet", ['pt','eta', 'phi'], delphes_genJets)

        # add JEC info
        allRecoJets = delphesReader.jets()
        # add btag info
        for i_btagWP, btagWP in enumerate(btagWPs):
            count = 0
            for jet in allRecoJets:
                btag = ( jet["bTag"] & (2**i_btagWP) > 0 ) # Read b-tag bitmap
                jet["bTag_"+btagWP] = btag

        # read jets
        recoJets =  filter( lambda j: isGoodRecoJet(j, max_jet_abseta = max_jet_abseta), allRecoJets) 
        recoJets.sort( key = lambda p:-p['pt'] )
        addIndex( recoJets )

        # count b-tag multiplicities
        for i_btagWP, btagWP in enumerate(btagWPs):
            count = 0
            for jet in recoJets:
                if jet["bTag_"+btagWP]: count += 1
            setattr( event, "nBTag_"+btagWP, count )
            if btagWP == default_btagWP:
                setattr( event, "nBTag", count )

        # upgrade JEC are flavor dependent
        for jet in allRecoJets:
            #btag_ = jet ["bTag_"+default_btagWP]
            jet["matchGenBJet"] = min( [999]+[ deltaR( jet, trueBJet ) for trueBJet in trueBjets ] )<0.4
            #upgradeJECUncertainty.applyJECInfo( jet, flavor = 5 if btag else 0 )

        jet_combinatorics = []

        if args.combinatoricalBTags:
            from TMB.Tools.delphesCombinatoricalBTagEfficiency import getBJetProbabilities
            probabilities_2b, recoJets = getBJetProbabilities( 
                filter( lambda j: j["matchGenBJet"], recoJets ), 
                filter( lambda j: (not j["matchGenBJet"]), recoJets ),
                nbtag_min = 2)
            for bjet_indices, probability in probabilities_2b.iteritems(): 
                #event.combinatoricalBTagWeight2b = probability
                recoBJets    = [ jet for i_jet, jet in enumerate(recoJets) if i_jet in bjet_indices ]
                recoNonBJets = [ jet for i_jet, jet in enumerate(recoJets) if not (i_jet in bjet_indices) ]
                #print probability, bjet_indices, recoBJets
                jet_combinatorics.append( (probability, recoBJets, recoNonBJets) )
        else:
            # make reco b jets
            recoBJets    = filter( lambda j:     j['bTag_'+default_btagWP] and abs(j['eta'])<2.4 , recoJets )
            recoNonBJets = filter( lambda j:not (j['bTag_'+default_btagWP] and abs(j['eta'])<2.4), recoJets )

            jet_combinatorics.append( (1., recoBJets, recoNonBJets) )

        for probability, recoBJets, recoNonBJets in jet_combinatorics:
            if probability==0: continue

            if args.combinatoricalBTags:
                event.combinatoricalBTagWeight2b = probability
                event.nBTag = 2  
                for j in recoBJets:
                    j['bTag'] = True
                for j in recoNonBJets:
                    j['bTag'] = False
 
            recoBj0, recoBj1 = ( recoBJets + recoNonBJets + [None, None] )[:2]

            if recoBj0: fill_vector( event, "recoBj0", recoJet_varnames, recoBj0)
            if recoBj1: fill_vector( event, "recoBj1", recoJet_varnames, recoBj1) 

            # VBS jets
            combs = list(itertools.combinations(range(len(recoJets)),2))
            combs.sort(key=lambda comb:-recoJets[comb[0]]['pt']-recoJets[comb[1]]['pt'])
            for i_j1, i_j2 in combs:
                if recoJets[i_j1]['eta']*recoJets[i_j2]['eta']<0:
                    vbsRecoJet0 = recoJets[i_j1]
                    vbsRecoJet1 = recoJets[i_j2]
                    event.recoVBSj0_index, event.recoVBSj1_index = i_j1, i_j2
                    fill_vector( event, "recoVBSj0", recoJet_varnames, vbsRecoJet0)
                    fill_vector( event, "recoVBSj1", recoJet_varnames, vbsRecoJet1)
                    event.recoVBS_dEta= vbsRecoJet0['eta']-vbsRecoJet1['eta']
                    event.recoVBS_dPhi= deltaPhi(vbsRecoJet0['phi'],vbsRecoJet1['phi'])
                    event.recoVBS_mjj = sqrt( 2*vbsRecoJet0['pt']*vbsRecoJet1['pt']*(cosh(event.recoVBS_dEta)-cos(event.recoVBS_dPhi)))
                    break

            # read leptons
            allRecoLeps = delphesReader.muons() + delphesReader.electrons()
            allRecoLeps.sort( key = lambda p:-p['pt'] )
            recoLeps =  filter( isGoodRecoLepton, allRecoLeps )

            #delphesGenLeptons = filter( lambda p: abs(p['pdgId']) in [11,13] and p['status']==1, delphesReader.genParticles() )
            # gen-match leptons with delphes particles
            for recoLep in allRecoLeps:
                #recoLep['genMatched'] = any( deltaR( recoLep, genLep )<0.1 for genLep in delphesGenLeptons )
                recoLep['genMatched'] = any( deltaR( recoLep, genLep )<0.1 for genLep in genLeps_dict )
                #print recoLep, recoLep['genMatched'], recoLep['genMatched2']
                
                #print recoLep['genMatched'], [deltaR( recoLep, genLep ) for genLep in delphesGenLeptons], recoLep

            # Photons
            recoPhotons = filter( isGoodRecoPhoton, delphesReader.photons() )

            # Remove radiated photons in dR cone
            for recoPhoton in recoPhotons:
                recoPhoton['minLeptonDR'] = 999 
                recoPhoton['minLeptonPt'] = -1.
                dr_values = [deltaR(recoPhoton, l) for l in allRecoLeps]
                recoPhoton['minLeptonDR'] = min([999]+dr_values) 
                if len( dr_values )>0:
                    closest_lepton = dr_values.index( min(dr_values) ) 
                    recoPhoton['minLeptonPt'] = allRecoLeps[closest_lepton]['pt'] 
            recoPhotons = list(filter( lambda g: g['minLeptonDR']>0.4, recoPhotons))
            for recoPhoton in recoPhotons:
                recoPhoton['minJetDR'] =  min([999]+[deltaR(recoPhoton, j) for j in recoJets])
            recoPhotons = list(filter( lambda g: g['minJetDR']>0.4, recoPhotons))

            # cross-cleaning of reco-objects
            nrecoLeps_uncleaned = len( recoLeps )
            recoLeps = filter( lambda l: (min([999]+[deltaR2(l, j) for j in recoJets if j['pt']>30]) > 0.4**2 ), recoLeps )
            #logger.info( "Before photon cleaning: %i after: %i allRecoLeps: %i, recoLeps %i", nrecoLeps_uncleaned, len(recoLeps), len( allRecoLeps ), len( recoLeps ) )

            # give index to leptons
            addIndex( recoLeps )
        
            # lepton uncertainties
            event.reweight_id_mu = 1.
            event.reweight_id_ele = 1.
            for l in recoLeps:
                if abs(l['pdgId'])==11:
                    event.reweight_id_ele*=1.005
                elif abs(l['pdgId'])==13:
                    event.reweight_id_mu*=1.005
            # MET
            recoMet = delphesReader.met()[0]

            # reco-bjet/leading lepton association
            if len(recoLeps)>0 and recoBj0 and recoBj1:
                if vecSumPt( recoBj0, recoLeps[0], recoMet ) > vecSumPt( recoBj1, recoLeps[0], recoMet ):
                    event.recoBjLeadlep_index, event.recoBjLeadhad_index = recoBj0['index'], recoBj1['index']
                else:
                    event.recoBjLeadlep_index, event.recoBjLeadhad_index = recoBj1['index'], recoBj0['index']

            # Store
            fill_vector_collection( event, "recoLep",    recoLep_varnames, recoLeps )
            fill_vector_collection( event, "recoJet",    recoJet_varnames, recoJets )
            fill_vector_collection( event, "recoPhoton", recoPhoton_varnames, recoPhotons )

            event.recoMet_pt  = recoMet['pt']
            event.recoMet_phi = recoMet['phi']

            delphesGenMet = delphesReader.genMet()[0]
            event.delphesGenMet_pt  = delphesGenMet['pt']
            event.delphesGenMet_phi = delphesGenMet['phi']

            # search for reco Z in reco leptons
            (event.recoZ_mass, recoZ_l1_index, recoZ_l2_index) = closestOSDLMassToMZ(recoLeps)
            recoNonZ_indices = [ i for i in range(len(recoLeps)) if i not in [recoZ_l1_index, recoZ_l2_index] ]
            event.recoZ_l1_index    = recoLeps[recoZ_l1_index]['index'] if recoZ_l1_index>=0 else -1
            event.recoZ_l2_index    = recoLeps[recoZ_l2_index]['index'] if recoZ_l2_index>=0 else -1
            event.recoNonZ_l1_index = recoLeps[recoNonZ_indices[0]]['index'] if len(recoNonZ_indices)>0 else -1
            event.recoNonZ_l2_index = recoLeps[recoNonZ_indices[1]]['index'] if len(recoNonZ_indices)>1 else -1

            # Store Z information 
            if event.recoZ_mass>0:
                if recoLeps[event.recoZ_l1_index]['pdgId']*recoLeps[event.recoZ_l2_index]['pdgId']>0 or abs(recoLeps[event.recoZ_l1_index]['pdgId'])!=abs(recoLeps[event.recoZ_l2_index]['pdgId']): 
                    raise RuntimeError( "not a Z! Should not happen" )
                Z_l1 = ROOT.TLorentzVector()
                Z_l1.SetPtEtaPhiM(recoLeps[event.recoZ_l1_index]['pt'], recoLeps[event.recoZ_l1_index]['eta'], recoLeps[event.recoZ_l1_index]['phi'], 0 )
                Z_l2 = ROOT.TLorentzVector()
                Z_l2.SetPtEtaPhiM(recoLeps[event.recoZ_l2_index]['pt'], recoLeps[event.recoZ_l2_index]['eta'], recoLeps[event.recoZ_l2_index]['phi'], 0 )
                Z = Z_l1 + Z_l2
                event.recoZ_pt   = Z.Pt()
                event.recoZ_eta  = Z.Eta()
                event.recoZ_phi  = Z.Phi()
                event.recoZ_lldPhi = deltaPhi(recoLeps[event.recoZ_l1_index]['phi'], recoLeps[event.recoZ_l2_index]['phi'])
                event.recoZ_lldR   = deltaR(recoLeps[event.recoZ_l1_index], recoLeps[event.recoZ_l2_index])
                lm_index = event.recoZ_l1_index if recoLeps[event.recoZ_l1_index]['pdgId'] > 0 else event.recoZ_l2_index
                event.recoZ_cosThetaStar = cosThetaStar(event.recoZ_mass, event.recoZ_pt, event.recoZ_eta, event.recoZ_phi, recoLeps[lm_index]['pt'], recoLeps[lm_index]['eta'], recoLeps[lm_index]['phi'] )

                # reco-bjet/lepton association
                if event.recoNonZ_l1_index>=0 and recoBj0 and recoBj1:
                    if vecSumPt( recoBj0, recoLeps[event.recoNonZ_l1_index], recoMet ) > vecSumPt( recoBj1, recoLeps[event.recoNonZ_l1_index], recoMet ):
                        event.recoBjNonZlep_index, event.recoBjNonZhad_index = recoBj0['index'], recoBj1['index']
                    else:
                        event.recoBjNonZlep_index, event.recoBjNonZhad_index = recoBj1['index'], recoBj0['index']

            # H->bb
            H_vecP4 = None
            if len(recoBJets)>=2:
                addTLorentzVector( recoBJets[0] )
                addTLorentzVector( recoBJets[1] )
                event.H_dijet_mass = (recoBJets[0]['vecP4'] + recoBJets[1]['vecP4']).M()
                H_vecP4          = recoBJets[0]['vecP4'] + recoBJets[1]['vecP4']
                event.H_pt       = H_vecP4.Pt()
                event.H_j1_index = recoBJets[0]['index']
                event.H_j2_index = recoBJets[1]['index']


            # store WH information
            if len(recoLeps)>=1:
                WH_lepton       = recoLeps[0]
                addTLorentzVector( WH_lepton )
                random_no       = gRandom.Uniform(0,1)
                neutrino_vecP4  = VV_angles.neutrino_mom(WH_lepton['vecP4'], recoMet['pt'], recoMet['phi'], random_no)
                W_vecP4         = neutrino_vecP4 + WH_lepton['vecP4']
                event.WH_nu_pt, event.WH_nu_eta, event.WH_nu_phi, event.WH_nu_E = neutrino_vecP4.Pt(), neutrino_vecP4.Eta(), neutrino_vecP4.Phi(), neutrino_vecP4.E()
                event.WH_W_pt       = W_vecP4.Pt()
                event.WH_dPhiMetLep = abs(deltaPhi( recoMet['phi'], WH_lepton['phi'] ))
                event.WH_MT         = sqrt(2.*recoMet['pt']*WH_lepton['pt']*(1-cos(event.WH_dPhiMetLep)))
                if H_vecP4 is not None:
                    event.WH_Theta = VV_angles.getTheta(WH_lepton['vecP4'], neutrino_vecP4, H_vecP4)
                    event.WH_theta = VV_angles.gettheta(WH_lepton['vecP4'], neutrino_vecP4, H_vecP4)
                    event.WH_phi   = VV_angles.getphi(  WH_lepton['vecP4'], neutrino_vecP4, H_vecP4)

            # store ZH information
            if event.recoZ_pt>0:
                lepton1 = recoLeps[event.recoZ_l1_index]
                lepton2 = recoLeps[event.recoZ_l2_index]
                addTLorentzVector( lepton1 )
                addTLorentzVector( lepton2 )

                # first one should have negative charge 
                if lepton1['pdgId']<0:
                    lepton1, lepton2 = lepton2, lepton1

                event.ZH_l1_index, event.ZH_l2_index = lepton1['index'], lepton2['index']

                if H_vecP4 is not None:

                    event.ZH_Theta = VV_angles.getTheta(lepton1['vecP4'], lepton2['vecP4'], H_vecP4)
                    event.ZH_theta = VV_angles.gettheta(lepton1['vecP4'], lepton2['vecP4'], H_vecP4)
                    event.ZH_phi   = VV_angles.getphi(  lepton1['vecP4'], lepton2['vecP4'], H_vecP4)
            maker.fill()
        maker.event.init()

counter = 0
for reader in readers:
    reader.start()
maker.start()

# While loop with first reader, all others 'run' too.

#if args.delphesEra is not None:
#    readers[0].run( )
#    logger.warning( "Skip first GEN event to align with delphes ntuple. This is a workaround for a bug (?) in Delphes that makes it skipm the first event...")
#    if len(sample.files)!=1:
#        raise RuntimeError ("Delphes skip first event means it will only work on a single file!!") 

while readers[0].run( ):
    for reader in readers[1:]:
        reader.run()

    filler( maker.event )
    #if abs(map( lambda p: p.daughter(0).pdgId(), filter( lambda p: p.pdgId()==23 and p.numberOfDaughters()==2, fwliteReader.products['gp']))[0])==13: 
    #    maker.run()
    #    break
    #maker.run()
         
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
