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
from math                                import sqrt, cos, sin, pi, acos
import imp

#RootTools
from RootTools.core.standard             import *

# TMB
from TMB.Tools.user                   import skim_output_directory
from TMB.Tools.GenSearch              import GenSearch
from TMB.Tools.helpers                import deltaPhi, deltaR, deltaR2, cosThetaStar, closestOSDLMassToMZ, checkRootFile
from TMB.Tools.HyperPoly              import HyperPoly
from TMB.Tools.WeightInfo             import WeightInfo
#from TMB.Tools.DelphesProducer        import DelphesProducer
from TMB.Tools.objectSelection        import isGoodGenJet, isGoodGenLepton, isGoodGenPhoton, isGoodRecoMuon, isGoodRecoElectron, isGoodRecoLepton, isGoodRecoJet, isGoodRecoPhoton, genJetId
#from TMB.Tools.UpgradeJECUncertainty  import UpgradeJECUncertainty 
#
# Arguments
# 
import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--logLevel',           action='store',      default='INFO',          nargs='?', choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'], help="Log level for logging")
argParser.add_argument('--small',              action='store_true', help='Run only on a small subset of the data?')#, default = True)
#argParser.add_argument('--delphes',            action='store_true', help='Run Delphes?')
argParser.add_argument('--overwrite',          action='store',      nargs='?', choices = ['none', 'all', 'target'], default = 'none', help='Overwrite?')#, default = True)
argParser.add_argument('--targetDir',          action='store',      default='v5')
argParser.add_argument('--sample',             action='store',      default='fwlite_ttZ_ll_LO_scan', help="Name of the sample loaded from fwlite_benchmarks. Only if no inputFiles are specified")
argParser.add_argument('--inputFiles',         action='store',      nargs = '*', default=[])
#argParser.add_argument('--delphesEra',         action='store',      default='delphes_card_CMS', choices = ["RunII", "RunIICentral", "RunIInoDelphesIso", "RunIIPileUp", "PhaseII"], help="specify delphes era")
argParser.add_argument('--targetSampleName',   action='store',      default=None, help="Name of the sample in case inputFile are specified. Otherwise ignored")
argParser.add_argument('--nJobs',              action='store',      nargs='?', type=int, default=1,  help="Maximum number of simultaneous jobs.")
argParser.add_argument('--job',                action='store',      nargs='?', type=int, default=0,  help="Run only job i")
argParser.add_argument('--addReweights',       action='store_true',   help="Add reweights?")
#argParser.add_argument('--removeDelphesFiles', action='store_true',   help="remove Delphes file after postprocessing?")
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
    sample_file = "$CMSSW_BASE/python/TMB/Samples/gen.py"
    samples = imp.load_source( "samples", os.path.expandvars( sample_file ) )
    sample = getattr( samples, args.sample )
    logger.debug( 'Loaded sample %s with %i files.', sample.name, len(sample.files) )

maxEvents = -1
if args.small: 
    args.targetDir += "_small"
    maxEvents       = 1000 
    sample.files=sample.files[:1]

xsec = sample.xsec
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

    # weights for the ntuple
    rw_vector       = TreeVariable.fromString( "rw[w/F,"+",".join(w+'/F' for w in weightInfo.variables)+"]")
    rw_vector.nMax  = weightInfo.nid
    extra_variables.append(rw_vector)

    # coefficients for the weight parametrization
    param_vector      = TreeVariable.fromString( "p[C/F]" )
    param_vector.nMax = HyperPoly.get_ndof(weightInfo.nvar, args.interpolationOrder)
    hyperPoly         = HyperPoly( args.interpolationOrder )
    extra_variables.append(param_vector)
    extra_variables.append(TreeVariable.fromString( "chi2_ndof/F"))

    def interpret_weight(weight_id):
        str_s = weight_id.split('_')
        res={}
        for i in range(len(str_s)/2):
            res[str_s[2*i]] = float(str_s[2*i+1].replace('m','-').replace('p','.'))
        return res

# Run only job number "args.job" from total of "args.nJobs"
if args.nJobs>1:
    n_files_before = len(sample.files)
    sample = sample.split(args.nJobs)[args.job]
    n_files_after  = len(sample.files)
    logger.info( "Running job %i/%i over %i files from a total of %i.", args.job, args.nJobs, n_files_after, n_files_before)

products = {
    'lhe':{'type':'LHEEventProduct', 'label':("externalLHEProducer")},
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
jet_write_vars      = jet_read_vars+',matchBParton/I' 
jet_write_varnames  =  varnames( jet_write_vars )
variables += ["genJet[%s]"%jet_write_vars]
variables += ["genBj0_%s"%var for var in jet_write_vars.split(',')]
variables += ["genBj1_%s"%var for var in jet_write_vars.split(',')]
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

# to be stored for each boson
boson_varnames = [ 'pt', 'phi', 'eta', 'mass', 'status']
# Z vector from gen collection
boson_all_varnames = boson_varnames + ['cosThetaStar', 'daughter_pdgId','l1_index', 'l2_index']
variables     += ["genZ[pt/F,phi/F,eta/F,mass/F,status/I,cosThetaStar/F,daughter_pdgId/I,l1_index/I,l2_index/I]"]
variables     += ["genW[pt/F,phi/F,eta/F,mass/F,status/I,cosThetaStar/F,daughter_pdgId/I,l1_index/I,l2_index/I]"]
# Z vector from genleps
# gamma vector
#gen_photon_vars = "pt/F,phi/F,eta/F,mass/F,mother_pdgId/I,relIso04/F,minLeptonDR/F,minJetDR/F"
#variables     += ["genPhoton[%s]"%gen_photon_vars]
#gen_photon_varnames = varnames( gen_photon_vars )

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

#if args.delphes:
#    if args.delphesEra == 'RunII':
#        from TMB.Tools.DelphesReader          import DelphesReader
#        delphesCard = 'delphes_card_CMS'
#    elif args.delphesEra == 'RunIICentral':
#        from TMB.Tools.DelphesReader          import DelphesReader
#        delphesCard = 'delphes_card_CMS_Central'
#    elif args.delphesEra == 'RunIInoDelphesIso':
#        from TMB.Tools.DelphesReader          import DelphesReader
#        delphesCard = 'delphes_card_CMS_noLepIso'
#    elif args.delphesEra == 'RunIIPileUp':
#        from TMB.Tools.DelphesReader          import DelphesReader
#        delphesCard = 'delphes_card_CMS_PileUp'
#    elif args.delphesEra == 'PhaseII':
#        from TMB.Tools.DelphesReaderCMSHLLHC  import DelphesReader
#        delphesCard = 'CMS_PhaseII/CMS_PhaseII_200PU_v03'
readers = []

# FWLite reader if this is an EDM file
fwliteReader = sample.fwliteReader( products = products )
readers.append( fwliteReader )

## Delphes reader if we read Delphes
#if args.delphes:
#    delphes_file = os.path.join( output_directory, 'delphes', sample.name+'.root' )
#    if not os.path.exists( delphes_file ) or args.overwrite in ['all']:
#        logger.debug( "Reproducing delphes file %s", delphes_file)
#        delphesProducer = DelphesProducer( card = delphesCard )
#        delphesProducer.produce( sample.files, delphes_file)
#    delphesReader = DelphesReader( Sample.fromFiles( delphes_file, delphes_file, treeName = "Delphes" ) ) # RootTools version
#    readers.append( delphesReader )

def filler( event ):

    event.lumiweight1fb = lumiweight1fb

    event.run, event.lumi, event.evt = fwliteReader.evt
    if fwliteReader.position % 100==0: logger.info("At event %i/%i", fwliteReader.position, fwliteReader.nEvents)

    if args.addReweights:
        event.nrw = weightInfo.nid
        lhe_weights = fwliteReader.products['lhe'].weights()
        weights      = []
        param_points = []
        for weight in lhe_weights:
            # Store nominal weight (First position!) 
            if weight.id in ['rwgt_1','dummy']: event.rw_nominal = weight.wgt
            if not weight.id in weightInfo.id: continue
            pos = weightInfo.data[weight.id]
            event.rw_w[pos] = weight.wgt
            weights.append( weight.wgt )
            interpreted_weight = interpret_weight(weight.id) 
            for var in weightInfo.variables:
                getattr( event, "rw_"+var )[pos] = interpreted_weight[var]
            # weight data for interpolation
            if not hyperPoly.initialized: param_points.append( tuple(interpreted_weight[var] for var in weightInfo.variables) )

        # get list of values of ref point in specific order
        ref_point_coordinates = [weightInfo.ref_point_coordinates[var] for var in weightInfo.variables]

        # Initialize with Reference Point
        if not hyperPoly.initialized: hyperPoly.initialize( param_points, ref_point_coordinates )
        coeff = hyperPoly.get_parametrization( weights )

        # = HyperPoly(weight_data, args.interpolationOrder)
        event.np = hyperPoly.ndof
        event.chi2_ndof = hyperPoly.chi2_ndof(coeff, weights)
        #logger.debug( "chi2_ndof %f coeff %r", event.chi2_ndof, coeff )
        if event.chi2_ndof>10**-6: logger.warning( "chi2_ndof is large: %f", event.chi2_ndof )
        for n in xrange(hyperPoly.ndof):
            event.p_C[n] = coeff[n]

        # lumi weight / w0
        event.ref_lumiweight1fb = event.lumiweight1fb / coeff[0]

    # MET
    genMet = {'pt':fwliteReader.products['genMET'][0].pt(), 'phi':fwliteReader.products['genMET'][0].phi()}
    event.genMet_pt  = genMet['pt']
    event.genMet_phi = genMet['phi'] 

    # All gen particles
    gp      = fwliteReader.products['gp']

    # for searching
    search  = GenSearch( gp )

    # find heavy objects before they decay
    genTops = map( lambda t:{var: getattr(t, var)() for var in top_varnames}, filter( lambda p:abs(p.pdgId())==6 and search.isLast(p),  gp) )

    genTops.sort( key = lambda p:-p['pt'] )
    fill_vector_collection( event, "genTop", top_varnames, genTops ) 

    # generated leptons
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


    # generated Z's
    genZs = filter( lambda p:abs(p.pdgId())==23 and search.isLast(p) and abs(p.daughter(0).pdgId()) in [11, 13, 15], gp)
    genZs.sort( key = lambda p: -p.pt() )
    genZs_dict = [ {var: getattr(genZ, var)() for var in boson_varnames} for genZ in genZs ]
    for i_genZ, genZ in enumerate(genZs):
        d1, d2 = genZ.daughter(0), genZ.daughter(1)
        if d1.pdgId()>0: 
            lm, lp = d1, d2
        else:
            lm, lp = d2, d1
        genZs_dict[i_genZ]['daughter_pdgId'] = lm.pdgId()
        genZs_dict[i_genZ]['l1_index'] = genLeps_from_bosons.index( genZ.daughter(0) )
        genZs_dict[i_genZ]['l2_index'] = genLeps_from_bosons.index( genZ.daughter(1) )
        genZs_dict[i_genZ]['cosThetaStar'] = cosThetaStar(genZ.mass(), genZ.pt(), genZ.eta(), genZ.phi(), lm.pt(), lm.eta(), lm.phi())
    fill_vector_collection( event, "genZ", boson_all_varnames, genZs_dict ) 

    # generated W's
    genWs = filter( lambda p:abs(p.pdgId())==24 and search.isLast(p) and abs(p.daughter(0).pdgId()) in [11, 12, 13, 14, 15, 16], gp)
    genWs.sort( key = lambda p: -p.pt() )
    genWs_dict = [ {var: getattr(genW, var)() for var in boson_varnames} for genW in genWs ]
    for i_genW, genW in enumerate(genWs):
        d1, d2 = genW.daughter(0), genW.daughter(1)
        if abs(d1.pdgId()) in [11, 113, 15]: 
            l, nu = d1, d2
        else:
            nu, l = d2, d1
        genWs_dict[i_genW]['daughter_pdgId'] = l.pdgId()
        genWs_dict[i_genW]['l1_index'] = genLeps_from_bosons.index( genW.daughter(0) )
        genWs_dict[i_genW]['l2_index'] = genLeps_from_bosons.index( genW.daughter(1) )
        genWs_dict[i_genW]['cosThetaStar'] = cosThetaStar(genW.mass(), genW.pt(), genW.eta(), genW.phi(), l.pt(), l.eta(), l.phi())
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
#
##    for ttgamma and tt events, categorize events:
##        a1) ttgamma vertex, gamma from t/g, isolated, must be in ttgamma sample   = ttgamma signal (photonSignal == 1)
##            -> gamma isolated
##            -> gamma with only tops/gluons in ancestry
##        a2) ttgamma, gamma from ISR, must be in ttgamma sample                    = tt ISR bg (photonISR == 1)
##            -> gamma isolated
##            -> gamma with no top in ancestry
##            -> gamma with only gluons and light quarks in ancestry
##            -> direct gamma mother != gluon!!! (this is a1)
##        b) tt + isolated gamma from W, l, tau, must be in tt sample               = ttgamma (W,l,tau) bg (photonLep == 1)
##            -> gamma isolated
##            -> gamma with direct abs mother being 11, 13, 15 or 24 
##        c1) tt + non isolated gamma from anything including mesons                = tt bg (ttBg == 1)
##            -> gamma non isolated or meson in ancestry
##        c2) tt + isolated gamma from bottom quark (from t decay) or jets (from W) = tt bottom bg (ttBg == 1)
##            -> gamma isolated from b or j (from t/W decay)
##            ATTENTION: gammas from bottoms with off-shell tops where the
##                       top decay is done by pythia are currently labeled as ISR!
##                       However this case is currently not simulated in any sample 
##        d) tt + gamma fake                                                        = ttgamma fake (photonFake == 1)
##            -> everything else does not contain a photon
##            -> if it still passes selection: it is a photon fake
#
#        genPhotonsSignalCheck = [ (search.ascend(l), l, search.ancestry( search.ascend(l) )) for l in filter( lambda p:abs(p.pdgId())==22 and p.pt()>10 and abs(p.eta())<2.5 and search.isLast(p) and p.status()==1, gp) ]
#        genPhotonsSignalCheck.sort( key = lambda p: -p[1].pt() )
#
#        photonSignal = 0    #a1
#        photonISR = 0       #a2
#        photonLep = 0       #b
#        ttBg = 0            #c1
#        photonJets = 0      #c2
#        photonFake = 0      #d
#
#
#        if len( genPhotonsSignalCheck ) > 0:     
#            # check hardest photon with pT>13 and abs(eta)<2.5
#            first, last, ancestry = genPhotonsSignalCheck[0]
#            # get abs pdgIDs of ancestry
#            pdg_ids = filter( lambda p:p!=2212, [abs(particle.pdgId()) for particle in ancestry])
#            # check if particles are close by
#            close_particles = filter( lambda p: p!=last and p.pt()>5 and deltaR2( {'phi':last.phi(), 'eta':last.eta()}, {'phi':p.phi(), 'eta':p.eta()} )<0.2**2 , search.final_state_particles_no_neutrinos )
#        
#    #        print 'mothers pdg', pdg_ids
#    #        print 'close', [p.pdgId() for p in close_particles]
#    #        print 'first mother pdg', first.mother(0).pdgId()
#
#    #        if len(pdg_ids) < 999:
#    #            previous = first.mother(0)
#    #            for i, pdg in enumerate(pdg_ids):
#    #                try:
#    #                    print 'mother', i, previous.pdgId()
#    #                    previous = previous.mother(0)
#    #                except Exception as val:
#    #                    print val
#    #                    break
#
#            # deside the categories
#            if max(pdg_ids)>100 or len(close_particles) != 0:
#                #photon with meson in ancestry or non isolated -> cat c1)
#                ttBg = 1
#            elif abs(first.mother(0).pdgId()) in [ 11, 13, 15, 24 ]:
#                #isolated photon with W, l or tau direct mother -> cat b)
#                photonLep = 1
#    #        elif all( [ p in [ 6, 21 ] for p in pdg_ids ] ) or abs(first.mother(0).pdgId()) == 21: # not working for photons from top, as the gluons can come from light quarks
#            elif abs(first.mother(0).pdgId()) in [ 6, 21 ]:
#                #isolated photon with photon from top or gluon -> cat a1)
#                photonSignal = 1
#    #            printTopMothers()
#            elif all( [ p in [ 1, 2, 3, 4, 5, 21 ] for p in pdg_ids ] ):
#                #isolated photon with photon ancestry only containing light quarks or gluons (ISR) -> cat a1)
#                photonISR = 1
#            else:
#                #isolated gammas from bottoms originating from the top decay or jets from W -> cat c2) 
#                photonJets = 1
#
#        else:
#            # if events with photonFake == 1 pass selection: fake gamma -> cat d)
#            photonFake = 1
#
#        # if all flags are 0, it is an isolated gamma from a process I havn't thought of!
#        # should not be there! - check!
#        event.signalPhoton = photonSignal
#        event.isrPhoton    = photonISR
#        event.lepPhoton    = photonLep
#        event.nonIsoPhoton = ttBg
#        event.jetPhoton    = photonJets
#        event.fakePhoton   = photonFake
#
#        # gen photons: particle-level isolated gen photons
#        genPhotons = [ (search.ascend(l), l) for l in filter( lambda p:abs(p.pdgId())==22 and p.pt()>15 and search.isLast(p) and p.status()==1, gp) ]
#        genPhotons.sort( key = lambda p: -p[1].pt() )
#        genPhotons_   = []
#
#        for first, last in genPhotons[:100]: 
#            mother_pdgId = first.mother(0).pdgId() if first.numberOfMothers()>0 else 0
#            genPhoton_ = {var: getattr(last, var)() for var in boson_read_varnames}
#            # kinematic photon selection
#            if not isGoodGenPhoton( genPhoton_): continue
#            genPhoton_['mother_pdgId'] = mother_pdgId
#            genPhoton_['status']      = last.status()
#            
#            close_particles = filter( lambda p: p!=last and deltaR2( {'phi':last.phi(), 'eta':last.eta()}, {'phi':p.phi(), 'eta':p.eta()} )<0.4**2 , search.final_state_particles_no_neutrinos )
#            genPhoton_['relIso04'] = sum( [p.pt() for p in close_particles], 0) / last.pt()
#            # require isolation
#            if genPhoton_['relIso04']<0.4:
#                genPhotons_.append( genPhoton_ )
#        
#        # genLeptons: prompt gen-leptons 
#        genLeptons = [ (search.ascend(l), l) for l in filter( lambda p:abs(p.pdgId()) in [11, 13] and search.isLast(p) and p.pt()>=0 and p.status()==1,  gp) ]
#        promptGenLeps    = []
#        allGenLeps    = []
#        for first, last in genLeptons:
#            mother = first.mother(0) if first.numberOfMothers()>0 else None
#            if mother is not None:
#                mother_pdgId      = mother.pdgId()
#                mother_ascend     = search.ascend(mother)
#                grandmother       = mother_ascend.mother(0) if mother.numberOfMothers()>0 else None
#                grandmother_pdgId = grandmother.pdgId() if grandmother is not None else 0
#            else:
#                mother_pdgId = 0
#                grandmother_pdgId = 0 
#            genLep = {var: getattr(last, var)() for var in lep_varnames}
#            genLep['mother_pdgId']      = mother_pdgId
#            genLep['grandmother_pdgId'] = grandmother_pdgId
#            allGenLeps.append( genLep )
#            if abs(genLep['mother_pdgId']) in [ 11, 13, 15, 23, 24, 25 ]:
#                promptGenLeps.append(genLep )
#
#        # filter gen leptons
#        promptGenLeps =  list( filter( lambda l:isGoodGenLepton( l ), promptGenLeps ) )
#        promptGenLeps.sort( key = lambda p:-p['pt'] )
#        addIndex( promptGenLeps )
#
#        ## removing photons in dR cone leptons (radiation photons)
#        for genPhoton in genPhotons_:
#            genPhoton['minLeptonDR'] =  min([999]+[deltaR(genPhoton, l) for l in allGenLeps])
#        genPhotons_ = list(filter( lambda g: g['minLeptonDR']>0.4, genPhotons_))
#        addIndex( genPhotons_ )
#
    # jets
    fwlite_genJets = filter( genJetId, fwliteReader.products['genJets'] )
    genJets = map( lambda t:{var: getattr(t, var)() for var in jet_read_varnames}, filter( lambda j:j.pt()>30, fwlite_genJets) )
    # filter genJets
    genJets = list( filter( lambda j:isGoodGenJet( j ), genJets ) )
    ## cleaning of jets with isolated photons
    #genJets = list( filter( lambda j:min([999]+[deltaR2(j, p) for p in genPhotons_ ])>0.4**2, genJets))

    ## store minimum DR to jets
    #for genPhoton in genPhotons_:
    #    genPhoton['minJetDR'] =  min([999]+[deltaR(genPhoton, j) for j in genJets])

    # find b's from tops:
    b_partons = [ b for b in filter( lambda p:abs(p.pdgId())==5 and p.numberOfMothers()==1 and abs(p.mother(0).pdgId())==6,  gp) ]

    # store if gen-jet is DR matched to a B parton
    for genJet in genJets:
        genJet['matchBParton'] = ( min([999]+[deltaR2(genJet, {'eta':b.eta(), 'phi':b.phi()}) for b in b_partons]) < 0.2**2 )

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
#
#        # reco-bjet/leading lepton association
#        if len(promptGenLeps)>0 and genBj0 and genBj1:
#            if vecSumPt( genBj0, promptGenLeps[0], genMet ) > vecSumPt( genBj1, promptGenLeps[0], genMet ):
#                event.genBjLeadlep_index, event.genBjLeadhad_index = genBj0['index'], genBj1['index']
#            else:
#                event.genBjLeadlep_index, event.genBjLeadhad_index = genBj1['index'], genBj0['index']
#
#        # find Z in genLep
#        (event.genLepZ_mass, genLepZ_l1_index, genLepZ_l2_index) = closestOSDLMassToMZ(promptGenLeps)
#        genLepNonZ_indices = [ i for i in range(len(promptGenLeps)) if i not in [genLepZ_l1_index, genLepZ_l2_index] ]
#        event.genLepZ_l1_index    = promptGenLeps[genLepZ_l1_index]['index'] if genLepZ_l1_index>=0 else -1
#        event.genLepZ_l2_index    = promptGenLeps[genLepZ_l2_index]['index'] if genLepZ_l2_index>=0 else -1
#        event.genLepNonZ_l1_index = promptGenLeps[genLepNonZ_indices[0]]['index'] if len(genLepNonZ_indices)>0 else -1
#        event.genLepNonZ_l2_index = promptGenLeps[genLepNonZ_indices[1]]['index'] if len(genLepNonZ_indices)>1 else -1
#        # store genLepZ stuff
#        if event.genLepZ_mass>0:
#            genLepZ_l1 = ROOT.TLorentzVector()
#            genLepZ_l1.SetPtEtaPhiM(promptGenLeps[event.genLepZ_l1_index]['pt'], promptGenLeps[event.genLepZ_l1_index]['eta'], promptGenLeps[event.genLepZ_l1_index]['phi'], 0 )
#            genLepZ_l2 = ROOT.TLorentzVector()
#            genLepZ_l2.SetPtEtaPhiM(promptGenLeps[event.genLepZ_l2_index]['pt'], promptGenLeps[event.genLepZ_l2_index]['eta'], promptGenLeps[event.genLepZ_l2_index]['phi'], 0 )
#            genLepZ = genLepZ_l1 + genLepZ_l2
#            event.genLepZ_pt   = genLepZ.Pt()
#            event.genLepZ_eta  = genLepZ.Eta()
#            event.genLepZ_phi  = genLepZ.Phi()
#            event.genLepZ_lldPhi = deltaPhi(promptGenLeps[event.genLepZ_l1_index]['phi'], promptGenLeps[event.genLepZ_l2_index]['phi'])
#            event.genLepZ_lldR   = deltaR(promptGenLeps[event.genLepZ_l1_index], promptGenLeps[event.genLepZ_l2_index])
#            genLepMinus_index = event.genLepZ_l1_index if promptGenLeps[event.genLepZ_l1_index]['pdgId'] > 0 else event.genLepZ_l2_index
#            event.genLepZ_cosThetaStar = cosThetaStar(event.genLepZ_mass, event.genLepZ_pt, event.genLepZ_eta, event.genLepZ_phi, promptGenLeps[genLepMinus_index]['pt'], promptGenLeps[genLepMinus_index]['eta'], promptGenLeps[genLepMinus_index]['phi'] )
#
#        # reco-bjet/nonZ lepton association
#        if event.genLepNonZ_l1_index>=0 and genBj0 and genBj1:
#            if vecSumPt( genBj0, promptGenLeps[event.genLepNonZ_l1_index], genMet ) > vecSumPt( genBj1, promptGenLeps[event.genLepNonZ_l1_index], genMet ):
#                event.genBjNonZlep_index, event.genBjNonZhad_index = genBj0['index'], genBj1['index']
#            else:
#                event.genBjNonZlep_index, event.genBjNonZhad_index = genBj1['index'], genBj0['index']
#
#        #for jet in genJets:
#        #    print jet['isMuon'], jet['isElectron'], jet['isPhoton'], min([999]+[deltaR2(jet, l) for l in promptGenLeps if l['pt']>10]), jet
#
#        # jet/lepton disambiguation -> remove jets, because gen-jets cluster all leptons
#        #if args.logLevel == 'DEBUG':
#        #    for jet in filter( lambda j: not (min([999]+[deltaR2(j, l) for l in promptGenLeps if l['pt']>10]) > 0.3**2 ), genJets ):
#        #        logger.debug( "Filtered gen %f jet %r lep %r", sqrt((min([999]+[deltaR2(jet, l) for l in promptGenLeps if l['pt']>10]))), jet, [ (l['eta'], jet['pt']/l['pt']) for l in promptGenLeps] )
#        #        assert False, ""
#
#
#        fill_vector_collection( event, "genPhoton", gen_photon_varnames, genPhotons_ ) 
#        fill_vector_collection( event, "genLep", lep_all_varnames, promptGenLeps)

#    # Reco quantities
#    if args.delphes:
#        #delphesReader.event.GetEntry(fwliteReader.position-1 ) # do this differently now ...  
#
#        delphes_genJets = delphesReader.genJets()
#        fill_vector_collection( event, "delphesGenJet", ['pt','eta', 'phi'], delphes_genJets)
#
#        # add JEC info
#        allRecoJets = delphesReader.jets()
#        # add btag info
#        for i_btagWP, btagWP in enumerate(btagWPs):
#            count = 0
#            for jet in allRecoJets:
#                btag = ( jet["bTag"] & (2**i_btagWP) > 0 ) # Read b-tag bitmap
#                jet["bTag_"+btagWP] = btag
#
#        # read jets
#        recoJets =  filter( isGoodRecoJet, allRecoJets) 
#        recoJets.sort( key = lambda p:-p['pt'] )
#        addIndex( recoJets )
#
#        # count b-tag multiplicities
#        for i_btagWP, btagWP in enumerate(btagWPs):
#            count = 0
#            for jet in recoJets:
#                if jet["bTag_"+btagWP]: count += 1
#            setattr( event, "nBTag_"+btagWP, count )
#            if btagWP == default_btagWP:
#                setattr( event, "nBTag", count )
#
#        # upgrade JEC are flavor dependent
#        for jet in allRecoJets:
#            btag_ = jet ["bTag_"+default_btagWP]
#            #upgradeJECUncertainty.applyJECInfo( jet, flavor = 5 if btag else 0 )
#        
#        # count JEC varied jet multiplicities
#        event.nrecoJets_JEC_up = len( filter( lambda j: isGoodRecoJet(j, pt_var = 'pt_JEC_up'), allRecoJets ) ) 
#        event.nrecoJets_JEC_down = len( filter( lambda j: isGoodRecoJet(j, pt_var = 'pt_JEC_down'), allRecoJets ) ) 
#
#        event.nBTag_JEC_up   = len( filter( lambda j: j["bTag_"+default_btagWP] and isGoodRecoJet(j, pt_var = 'pt_JEC_up'), allRecoJets ) )
#        event.nBTag_JEC_down = len( filter( lambda j: j["bTag_"+default_btagWP] and isGoodRecoJet(j, pt_var = 'pt_JEC_down'), allRecoJets ) )
#        
#        # read jets
#        recoJets =  filter( isGoodRecoJet, allRecoJets ) 
#        recoJets.sort( key = lambda p:-p['pt'] )
#        addIndex( recoJets )
#
#        # make reco b jets
#        recoBJets    = filter( lambda j:j['bTag_'+default_btagWP], recoJets )
#        recoNonBJets = filter( lambda j:not j['bTag_'+default_btagWP], recoJets )
#        recoBj0, recoBj1 = ( recoBJets + recoNonBJets + [None, None] )[:2]
#        if recoBj0: fill_vector( event, "recoBj0", recoJet_varnames, recoBj0)
#        if recoBj1: fill_vector( event, "recoBj1", recoJet_varnames, recoBj1) 
#
#        # add b-tag reweights
#        event.reweight_BTag_B = getBTagSF_1a(default_btagWP, 'B', recoBJets, recoNonBJets )
#        event.reweight_BTag_L = getBTagSF_1a(default_btagWP, 'L', recoBJets, recoNonBJets )
#
#        # read leptons
#        allRecoLeps = delphesReader.muons() + delphesReader.electrons()
#        allRecoLeps.sort( key = lambda p:-p['pt'] )
#        recoLeps =  filter( isGoodRecoLepton, allRecoLeps )
#
#        delphesGenLeptons = filter( lambda p: abs(p['pdgId']) in [11,13] and p['status']==1, delphesReader.genParticles() )
#        # gen-match leptons with delphes particles
#        for recoLep in allRecoLeps:
#            recoLep['genMatched'] = any( deltaR( recoLep, genLep )<0.1 for genLep in delphesGenLeptons )
#            
#            #print recoLep['genMatched'], [deltaR( recoLep, genLep ) for genLep in delphesGenLeptons], recoLep
#
#        # Photons
#        recoPhotons = filter( isGoodRecoPhoton, delphesReader.photons() )
#
#        # Remove radiated photons in dR cone
#        for recoPhoton in recoPhotons:
#            recoPhoton['minLeptonDR'] = 999 
#            recoPhoton['minLeptonPt'] = -1.
#            dr_values = [deltaR(recoPhoton, l) for l in allRecoLeps]
#            recoPhoton['minLeptonDR'] = min([999]+dr_values) 
#            if len( dr_values )>0:
#                closest_lepton = dr_values.index( min(dr_values) ) 
#                recoPhoton['minLeptonPt'] = allRecoLeps[closest_lepton]['pt'] 
#        recoPhotons = list(filter( lambda g: g['minLeptonDR']>0.4, recoPhotons))
#        for recoPhoton in recoPhotons:
#            recoPhoton['minJetDR'] =  min([999]+[deltaR(recoPhoton, j) for j in recoJets])
#        recoPhotons = list(filter( lambda g: g['minJetDR']>0.4, recoPhotons))
#
#        # cross-cleaning of reco-objects
#        nrecoLeps_uncleaned = len( recoLeps )
#        recoLeps = filter( lambda l: (min([999]+[deltaR2(l, j) for j in recoJets if j['pt']>30]) > 0.3**2 ), recoLeps )
#        #logger.info( "Before photon cleaning: %i after: %i allRecoLeps: %i, recoLeps %i", nrecoLeps_uncleaned, len(recoLeps), len( allRecoLeps ), len( recoLeps ) )
#
#        # give index to leptons
#        addIndex( recoLeps )
#    
#        # lepton uncertainties
#        event.reweight_id_mu = 1.
#        event.reweight_id_ele = 1.
#        for l in recoLeps:
#            if abs(l['pdgId'])==11:
#                event.reweight_id_ele*=1.005
#            elif abs(l['pdgId'])==13:
#                event.reweight_id_mu*=1.005
#        # MET
#        recoMet = delphesReader.met()[0]
#
#        # reco-bjet/leading lepton association
#        if len(recoLeps)>0 and recoBj0 and recoBj1:
#            if vecSumPt( recoBj0, recoLeps[0], recoMet ) > vecSumPt( recoBj1, recoLeps[0], recoMet ):
#                event.recoBjLeadlep_index, event.recoBjLeadhad_index = recoBj0['index'], recoBj1['index']
#            else:
#                event.recoBjLeadlep_index, event.recoBjLeadhad_index = recoBj1['index'], recoBj0['index']
#
#        # Photons
#        for recoPhoton in recoPhotons:
#            recoPhoton['genIndex'] = -1
#            minDR = 999
#            for index, genPhoton in enumerate(genPhotons_):
#                dr = deltaR( recoPhoton, genPhoton )
#                if dr<0.4 and dr<minDR:
#                    minDR = dr
#                    recoPhoton['genIndex'] = index 
#        # Store
#        fill_vector_collection( event, "recoLep",    recoLep_varnames, recoLeps )
#        fill_vector_collection( event, "recoJet",    recoJet_varnames, recoJets )
#        fill_vector_collection( event, "recoPhoton", recoPhoton_varnames, recoPhotons )
#
#        event.recoMet_pt  = recoMet['pt']
#        event.recoMet_phi = recoMet['phi']
#
#        # search for reco Z in reco leptons
#        (event.recoZ_mass, recoZ_l1_index, recoZ_l2_index) = closestOSDLMassToMZ(recoLeps)
#        recoNonZ_indices = [ i for i in range(len(recoLeps)) if i not in [recoZ_l1_index, recoZ_l2_index] ]
#        event.recoZ_l1_index    = recoLeps[recoZ_l1_index]['index'] if recoZ_l1_index>=0 else -1
#        event.recoZ_l2_index    = recoLeps[recoZ_l2_index]['index'] if recoZ_l2_index>=0 else -1
#        event.recoNonZ_l1_index = recoLeps[recoNonZ_indices[0]]['index'] if len(recoNonZ_indices)>0 else -1
#        event.recoNonZ_l2_index = recoLeps[recoNonZ_indices[1]]['index'] if len(recoNonZ_indices)>1 else -1
#
#        # Store Z information 
#        if event.recoZ_mass>=0:
#            if recoLeps[event.recoZ_l1_index]['pdgId']*recoLeps[event.recoZ_l2_index]['pdgId']>0 or abs(recoLeps[event.recoZ_l1_index]['pdgId'])!=abs(recoLeps[event.recoZ_l2_index]['pdgId']): 
#                raise RuntimeError( "not a Z! Should not happen" )
#            Z_l1 = ROOT.TLorentzVector()
#            Z_l1.SetPtEtaPhiM(recoLeps[event.recoZ_l1_index]['pt'], recoLeps[event.recoZ_l1_index]['eta'], recoLeps[event.recoZ_l1_index]['phi'], 0 )
#            Z_l2 = ROOT.TLorentzVector()
#            Z_l2.SetPtEtaPhiM(recoLeps[event.recoZ_l2_index]['pt'], recoLeps[event.recoZ_l2_index]['eta'], recoLeps[event.recoZ_l2_index]['phi'], 0 )
#            Z = Z_l1 + Z_l2
#            event.recoZ_pt   = Z.Pt()
#            event.recoZ_eta  = Z.Eta()
#            event.recoZ_phi  = Z.Phi()
#            event.recoZ_lldPhi = deltaPhi(recoLeps[event.recoZ_l1_index]['phi'], recoLeps[event.recoZ_l2_index]['phi'])
#            event.recoZ_lldR   = deltaR(recoLeps[event.recoZ_l1_index], recoLeps[event.recoZ_l2_index])
#            lm_index = event.recoZ_l1_index if recoLeps[event.recoZ_l1_index]['pdgId'] > 0 else event.recoZ_l2_index
#            event.recoZ_cosThetaStar = cosThetaStar(event.recoZ_mass, event.recoZ_pt, event.recoZ_eta, event.recoZ_phi, recoLeps[lm_index]['pt'], recoLeps[lm_index]['eta'], recoLeps[lm_index]['phi'] )
#
#            # reco-bjet/lepton association
#            if event.recoNonZ_l1_index>=0 and recoBj0 and recoBj1:
#                if vecSumPt( recoBj0, recoLeps[event.recoNonZ_l1_index], recoMet ) > vecSumPt( recoBj1, recoLeps[event.recoNonZ_l1_index], recoMet ):
#                    event.recoBjNonZlep_index, event.recoBjNonZhad_index = recoBj0['index'], recoBj1['index']
#                else:
#                    event.recoBjNonZlep_index, event.recoBjNonZhad_index = recoBj1['index'], recoBj0['index']

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
    sequence  = [ filler ],
    variables = [ TreeVariable.fromString(x) for x in variables ] + extra_variables,
    treeName = "Events"
    )

tmp_dir.cd()

counter = 0
for reader in readers:
    reader.start()
maker.start()

# While loop with first reader, all others 'run' too.
while readers[0].run( ):
    for reader in readers[1:]:
        reader.run()
    #if abs(map( lambda p: p.daughter(0).pdgId(), filter( lambda p: p.pdgId()==23 and p.numberOfDaughters()==2, fwliteReader.products['gp']))[0])==13: 
    #    maker.run()
    #    break
    maker.run()
         
    counter += 1
    if counter == maxEvents:  break

logger.info( "Done with running over %i events.", readers[0].nEvents )

output_file.cd()
maker.tree.Write()
output_file.Close()

logger.info( "Written output file %s", output_filename )

##cleanup delphes file:
#if os.path.exists( output_filename ) and args.delphes and args.removeDelphesFiles:
#    os.remove( delphes_file )
#    logger.info( "Removing Delphes file %s", delphes_file )
