# WeightInfo

from    Analysis.Tools.WeightInfo   import WeightInfo
import  Analysis.Tools.syncer       as syncer
from    Analysis.Tools.helpers      import deltaPhi, deltaR, getObjDict

from    RootTools.core.standard     import *

from    TMB.Tools.helpers           import getCollection, mZ
from    TMB.Tools.user              import plot_directory
import  TMB.Tools.VV_angles         as VV_angles
from    TMB.Tools.delphesCutInterpreter import cutInterpreter

import  ROOT
import  os 
from    math                        import sqrt, sin, cos, sinh, cosh, copysign, pi

# Arguments
import argparse
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--logLevel',           action='store',      default='INFO',          nargs='?', choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'], help="Log level for logging")
argParser.add_argument('--plot_directory',     action='store',      default='delphes')
argParser.add_argument('--selection',          action='store',      default='singlelep-WHJet')
argParser.add_argument('--signal',             action='store',      default='WH', choices = ['WH', 'ZH'])
argParser.add_argument('--small',                                   action='store_true',     help='Run only on a small subset of the data?')
args = argParser.parse_args()

# Logger
import TMB.Tools.logger as logger
import RootTools.core.logger as logger_rt
logger    = logger.get_logger(   args.logLevel, logFile = None)
logger_rt = logger_rt.get_logger(args.logLevel, logFile = None)

maxN = -1
plot_directory = os.path.join(plot_directory, args.plot_directory,  args.signal )
if args.small: 
    plot_directory += "_small"
    maxN = 10000

import  TMB.Samples.pp_gen_v10 as samples

sample = getattr( samples, args.signal )

sample.weightInfo = WeightInfo(sample.reweight_pkl)
sample.weightInfo.set_order(2)
sample.setSelectionString( cutInterpreter.cutString(args.selection) )

# read_variables
jetVars          = ['pt/F', 'eta/F', 'phi/F', 'bTag/F', 'bTagPhys/I']
jetVarNames      = [x.split('/')[0] for x in jetVars]
lepVars          = ['pt/F','eta/F','phi/F','pdgId/I','isolationVar/F', 'isolationVarRhoCorr/F']
lepVarNames      = [x.split('/')[0] for x in lepVars]
read_variables = [\
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
    "evt/l", "run/I", "lumi/I", "np/I", "nweight/I", 
]
read_variables += [VectorTreeVariable.fromString( "p[C/F]", nMax=200 )]
read_variables += [VectorTreeVariable.fromString( "weight[base/F]", nMax=200 )]

def addTLorentzVector( p_dict ):
    ''' add a TLorentz 4D Vector for further calculations
    '''
    p_dict['vecP4'] = ROOT.TLorentzVector( p_dict['pt']*cos(p_dict['phi']), p_dict['pt']*sin(p_dict['phi']),  p_dict['pt']*sinh(p_dict['eta']), p_dict['pt']*cosh(p_dict['eta']) )

sequence = []

from TMB.Tools.objectSelection import isBJet
def makeJets( event, sample ):
    event.jets     = [getObjDict(event, 'recoJet_', jetVarNames, i) for i in range(int(event.nrecoJet))]
    for p in event.jets:
        #addTransverseVector( p )
        addTLorentzVector( p )
    event.bJets    = filter(lambda j:j['bTag']>=1 and abs(j['eta'])<=2.4    , event.jets)

sequence.append( makeJets )

gRandom = ROOT.TRandom3()
def makeLeptonic( event, sample ):

    # Read leptons, do not yet filter
    event.all_leps = getCollection( event, 'recoLep', ['pt', 'eta', 'phi', 'pdgId', 'isolationVar', 'isolationVarRhoCorr'], 'nrecoLep' )
    # Add extra vectors
    for p in event.all_leps:
        #addTransverseVector( p )
        addTLorentzVector( p )

    # Sort
    event.leps = sorted( event.all_leps, key=lambda k: -k['pt'] )

    # Cross-cleaning: remove leptons that overlap with a jet within 0.4
    #before_cleaning = len(event.leps)
    event.leps = list(filter( lambda l: min( [ deltaR(l, j) for j in event.jets ] + [999] ) > 0.4 , event.leps ))
    #if sample.name == signal.name:
    #    print 
    #    print_leps ( "reco", event.leps ) 
    if sample.name == "WH":
    #print "Jet cleaning removed %i leps" %( before_cleaning-len(event.leps))
        if len(event.leps)>0:
            event.lepton = event.leps[0]
            random_no      = gRandom.Uniform(0,1)
            event.neutrino_vecP4 = VV_angles.neutrino_mom(event.lepton['vecP4'], event.recoMet_pt, event.recoMet_phi, random_no)
            event.V_vecP4  = event.neutrino_vecP4 + event.lepton['vecP4']

            #event.neutrino_vecP4_genMet = VV_angles.neutrino_mom(event.lepton['vecP4'], event.genMet_pt, event.genMet_phi, random_no)
            #event.V_vecP4_genMet        = event.neutrino_vecP4_genMet + event.lepton['vecP4']  

            #print "rand", random_no, "lep M", event.lepton['vecP4'].M(), "Pt", event.lepton['pt'], "phi", event.lepton['phi'], "eta", event.lepton['eta']
            #print event.recoMet_pt, event.neutrino_vecP4.Pt(), event.genMet_pt
            event.V_pt          = event.V_vecP4.Pt()
            event.has_highPt_V  = event.V_pt > 150
            event.dPhiMetLep    = abs(deltaPhi( event.recoMet_phi, event.lepton['phi'] ))
            event.MT            = sqrt(2.*event.recoMet_pt*event.lepton['pt']*(1-cos(event.dPhiMetLep)))
        else:
            event.lepton            = None
            event.neutrino_vecP4    = None
            event.V_vecP4           = None
            event.V_pt              = float('nan')
            event.has_highPt_V      = False
            event.dPhiMetLep        = float('nan')
            event.MT                = float('nan')

    elif sample.name == "ZH":

        #"recoZ_pt/F", "recoZ_eta/F", "recoZ_phi/F", "recoZ_mass/F", "recoZ_cosThetaStar/F", "recoZ_lldPhi/F", "recoZ_lldR/F", "recoZ_l1_index/I", "recoZ_l2_index/I",
        if event.recoZ_pt>0:
            event.V_pt          =   event.recoZ_pt
            event.V_vecP4       =   event.all_leps[event.recoZ_l1_index]['vecP4'] + event.all_leps[event.recoZ_l2_index]['vecP4']
            event.lepton1       =   event.all_leps[event.recoZ_l1_index]
            event.lepton2       =   event.all_leps[event.recoZ_l2_index]

            # first one should have negative charge 
            if event.lepton1['pdgId']<0:
                event.lepton1, event.lepton2 = event.lepton2, event.lepton1

            event.has_highPt_V  =   event.recoZ_pt > 75
        else:
            event.V_pt          =   float('nan')
            event.V_vecP4       =   False
            event.lepton1       =   None
            event.lepton2       =   None
            event.has_highPt_V  =   False

sequence.append( makeLeptonic )

def makeH( event, sample ):
    if len(event.bJets)>=2:
        event.dijet_mass = (event.bJets[0]['vecP4'] + event.bJets[1]['vecP4']).M()
        event.H_vecP4 = event.bJets[0]['vecP4'] + event.bJets[1]['vecP4']
        event.H_pt    = event.H_vecP4.Pt()
    else:
        event.dijet_mass = float('nan')
        event.H_vecP4 = None
        event.H_pt    = float('nan')
    if event.dijet_mass>90 and event.dijet_mass<150:
        event.has_H = 1
    else:
        event.has_H = 0

sequence.append( makeH )
 
def make_VV_angles( event, sample ):
    event.Theta = None 
    event.theta = None 
    event.phi   = None 
    if sample.name == "WH":
        if event.lepton is None: return
        v = [event.lepton['vecP4'], event.neutrino_vecP4 ]
        event.Theta = VV_angles.getTheta(event.lepton['vecP4'], event.neutrino_vecP4, event.H_vecP4)
        event.theta = VV_angles.gettheta(event.lepton['vecP4'], event.neutrino_vecP4, event.H_vecP4)
        event.phi   = VV_angles.getphi(  event.lepton['vecP4'], event.neutrino_vecP4, event.H_vecP4)
    elif sample.name == "ZH":
        if event.lepton1 is None or event.lepton2 is None: return
        event.Theta = VV_angles.getTheta(event.lepton1['vecP4'], event.lepton2['vecP4'], event.H_vecP4)
        event.theta = VV_angles.gettheta(event.lepton1['vecP4'], event.lepton2['vecP4'], event.H_vecP4)
        event.phi   = VV_angles.getphi(  event.lepton1['vecP4'],   event.lepton2['vecP4'], event.H_vecP4)

    event.pos_phi = event.neg_phi = float('nan')
    event.sin2thetaSin2Theta = sin(2*event.theta)*sin(2*event.Theta)

    if event.sin2thetaSin2Theta>0:
        event.pos_phi = event.phi
    else:
        event.neg_phi = event.phi

sequence.append( make_VV_angles )

# selection bools
if sample.name == 'WH':
    def makeSelection( event, sample):
        event.selection_realNu      = (event.neutrino_vecP4 is not None) and (event.neutrino_vecP4.E()>0.1)
        event.selection             = event.has_highPt_V and event.has_H and event.dPhiMetLep < 2
        event.selection_noH         = event.has_highPt_V and event.dPhiMetLep < 2 and event.H_vecP4 is not None
        event.selection_noHighPtV   = event.has_H and event.dPhiMetLep < 2
        event.selection_noPhiMetLep = event.has_highPt_V and event.has_H
elif sample.name == 'ZH':
    def makeSelection( event, sample):
        event.selection     = event.has_highPt_V and event.has_H
        event.selection_noH = event.has_highPt_V
        event.selection_noHighPtV = event.has_H

sequence.append( makeSelection )

#rw_str_sm = sample.weightInfo.get_weight_string()
#rw_str_1  = sample.weightInfo.get_weight_string(cHj3=1)
#rw_str_2  = sample.weightInfo.get_weight_string(cHj3=2)

rw_sm = sample.weightInfo.get_weight_func()
rw_cHj3_1  = sample.weightInfo.get_weight_func(cHj3=1)
rw_cHj3_2  = sample.weightInfo.get_weight_func(cHj3=2)

rw_cHW_1  = sample.weightInfo.get_weight_func(cHW=1)
rw_cHW_2  = sample.weightInfo.get_weight_func(cHW=2)
rw_cHWtil_1  = sample.weightInfo.get_weight_func(cHWtil=1)
rw_cHWtil_2  = sample.weightInfo.get_weight_func(cHWtil=2)

h_SM        = ROOT.TH1F("SM", "SM", 20, -pi, pi)
h_cHW_1     = ROOT.TH1F("cHW_1", "cHW_1", 20, -pi, pi)
h_cHW_2     = ROOT.TH1F("cHW_2", "cHW_2", 20, -pi, pi)
h_cHWtil_1  = ROOT.TH1F("cHWtil_1", "cHWtil_1", 20, -pi, pi)
h_cHWtil_2  = ROOT.TH1F("cHWtil_2", "cHWtil_2", 20, -pi, pi)

r = sample.treeReader(variables=read_variables)
r.start()

counter = 0
while r.run():

    for f in sequence:
        f( r.event, sample )


    #print "GEN pt W", r.event.genW_pt[0]#, r.event.p_C[0], r.event.p_C[2], r.event.p_C[33], r.event.p_C[2]/r.event.p_C[0], r.event.p_C[33]/r.event.p_C[0]
    #print "SM", rw_sm(r.event, sample), "cHj3=1", rw_cHj3_1(r.event, sample), "cHj3=2", rw_cHj3_2(r.event, sample)
    #print "SM", r.event.weight_base[0], "cHj3=1", r.event.weight_base[2], "cHj3=2", r.event.weight_base[33] 
    #print
    #print "SM", rw_sm(r.event, sample), "cHW=1", rw_cHW_1(r.event, sample), "cHW=2", rw_cHW_2(r.event, sample)
    #print "SM", r.event.weight_base[0], "cHW=1", r.event.weight_base[7], "cHW=2", r.event.weight_base[98]
    #print 
    #print "SM", rw_sm(r.event, sample), "cHWtil=1", rw_cHWtil_1(r.event, sample), "cHWtil=2", rw_cHWtil_2(r.event, sample)
    #print "SM", r.event.weight_base[0], "cHWtil=1", r.event.weight_base[8], "cHWtil=2", r.event.weight_base[108]
    #print

    counter+=1
    #if not r.event.selection: continue
    if r.event.phi is None: continue

    if sample.name == 'WH' and not r.event.selection_realNu: continue
    #print rw_sm( r.event, sample ), r.event.sin2thetaSin2Theta, copysign( rw_sm( r.event, sample ),  r.event.sin2thetaSin2Theta )

    h_SM.Fill       ( r.event.phi, copysign( rw_sm( r.event, sample ),  r.event.sin2thetaSin2Theta )) 
    h_cHW_1.Fill    ( r.event.phi, copysign( rw_cHW_1( r.event, sample ), r.event.sin2thetaSin2Theta )) 
    h_cHW_2.Fill    ( r.event.phi, copysign( rw_cHW_2( r.event, sample ), r.event.sin2thetaSin2Theta )) 
    h_cHWtil_1.Fill ( r.event.phi, copysign( rw_cHWtil_1( r.event, sample ), r.event.sin2thetaSin2Theta )) 
    h_cHWtil_2.Fill ( r.event.phi, copysign( rw_cHWtil_2( r.event, sample ), r.event.sin2thetaSin2Theta )) 

    if maxN>0 and counter>=maxN: break
   
h_SM        .SetLineColor(ROOT.kBlack)
h_cHW_1     .SetLineColor(ROOT.kBlue)
h_cHW_2     .SetLineColor(ROOT.kRed) 
h_cHWtil_1  .SetLineColor(ROOT.kGreen)
h_cHWtil_2  .SetLineColor(ROOT.kOrange)

h_SM.Scale(-1)
for h in [ h_cHW_2, h_cHW_1, h_cHWtil_2, h_cHWtil_1]:
    h.Add(h_SM)

c1 = ROOT.TCanvas()
same = ""
for h in [ h_cHW_2, h_cHW_1, h_cHWtil_2, h_cHWtil_1]:
    h.Draw("hist"+same)
    same    =   "same" 

c1.Print(os.path.join( plot_directory, "phi_shape_%s_maxN_%i.png"%(sample.name, maxN)) )

syncer.sync()


#lumi  = 137
#lumi_weight = lambda event, sample: lumi*event.lumiweight1fb
#
#eft_configs = [
#    (ROOT.kBlack, {}, "SM"),
#    (ROOT.kGreen+1, {'cHW':1}, "c_{W}=1"),
#    (ROOT.kOrange-1, {'cHj3':1}, "c_{Hq3}=1"),
#    (ROOT.kOrange-2, {'cHj3':.7}, "c_{Hq3}=.7"),
#    ]
#
#def get_eft_reweight( eft, weightInfo_):
#    func1 = sample.weightInfo.get_weight_func(**eft)
#    func2 = sample.weightInfo.get_weight_func()

