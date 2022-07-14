''' Class for reading Delphes files.
    Based on Suchitas class, now moved to RootTools reader
'''

# Standard imports
import ROOT

# Delphes Reader RootTools
from RootTools.core.DelphesReaderBase import DelphesReaderBase

class DelphesReader( DelphesReaderBase ): # version RootTools reader
#class DelphesReader: # version Suchtia-reader

#    def __init__( self, filename ):
#        ''' Return a Suchita-reader
#        '''
#        self.reader =  ROOT.DelphesRecoClass(filename)

    # Read a vector collection from the Delphes reader
#    def read_collection( self, collection, variables ):
#        ''' read delphes collection and rename leaves'''
#        nColl   = getattr( self.reader, collection+"_size" )
#        buffers = {var_old: getattr( self.reader, collection+'_'+var_old) for var_old, var_new in variables}
#        return [{var_new:buffers[var_old][i] for var_old, var_new in variables} for i in range(nColl)]

    def muons( self ):
        res = self.read_collection( 'Muon', 
            [   ('PT', 'pt'), ( 'Eta', 'eta'), ('Phi', 'phi'), 
                ('Charge', 'charge'), ('IsolationVar', 'isolationVar'), ('IsolationVarRhoCorr', 'isolationVarRhoCorr'),  
                ('SumPtCharged', 'sumPtCharged'),  ('SumPtNeutral', 'sumPtNeutral'), ('SumPtChargedPU', 'sumPtChargedPU'),  ('SumPt', 'sumPt') 
            ])
        for r in res:
            r['pdgId'] = -13*r['charge']
            r['ehadOverEem'] = float('nan')
        return res

    def electrons( self ):
        res = self.read_collection( 'Electron', 
            [   ('PT', 'pt'), ( 'Eta', 'eta'), ('Phi', 'phi'), 
                ('Charge', 'charge'), ('IsolationVar', 'isolationVar'), ('IsolationVarRhoCorr', 'isolationVarRhoCorr'),  
                ('SumPtCharged', 'sumPtCharged'),  ('SumPtNeutral', 'sumPtNeutral'), ('SumPtChargedPU', 'sumPtChargedPU'),  ('SumPt', 'sumPt'),
                ('EhadOverEem','ehadOverEem')
            ])
        for r in res:
            r['pdgId'] = -11*r['charge']
        return res

    def jets( self ):
        return self.read_collection( 'Jet', 
            [   ('PT', 'pt'), ( 'Eta', 'eta'), ('Phi', 'phi'),
                ('BTag', 'bTag'), ( 'BTagPhys', 'bTagPhys'), ('Flavor', 'flavor'),
                ('NCharged', 'nCharged'), ('NNeutrals', 'nNeutrals'), 
            ])

    def genJets( self ):
        return self.read_collection( 'GenJet', 
            [   ('PT', 'pt'), ( 'Eta', 'eta'), ('Phi', 'phi'),
            ])

    def photons( self ):
        return self.read_collection( 'Photon', 
            [   ('PT', 'pt'), ( 'Eta', 'eta'), ('Phi', 'phi'),
                ('IsolationVar', 'isolationVar'), ('IsolationVarRhoCorr', 'isolationVarRhoCorr'),  
                ('SumPtCharged', 'sumPtCharged'),  ('SumPtNeutral', 'sumPtNeutral'), ('SumPtChargedPU', 'sumPtChargedPU'),  ('SumPt', 'sumPt'),
                ('EhadOverEem','ehadOverEem') 
            ])

    def met( self ):
        return self.read_collection( 'MissingET', [('MET', 'pt'), ('Phi', 'phi')] )
    def genMet( self ):
        return self.read_collection( 'GenMissingET', [('MET', 'pt'), ('Phi', 'phi')] )

#    def genParticles( self ):
#        return self.read_collection( 'Particle', 
#            [   ('PT', 'pt'), ( 'Eta', 'eta'), ('Phi', 'phi'), ('Status', 'status'), ('PID', 'pdgId')
#            ])

#    def getEntry( self, entry ):
#        # Set the Suchi-Reader 
#        self.reader.GetEntry( entry )

#root [3] Delphes->GetListOfLeaves()->ls()
#OBJ: TObjArray  TObjArray   An array of objects : 0
# OBJ: TLeafElement  Event_  Event_ : 0 at: 0x58a61d0
# OBJ: TLeafElement  Event.fUniqueID fUniqueID[Event_] : 0 at: 0x58a2740
# OBJ: TLeafElement  Event.fBits fBits[Event_] : 0 at: 0x58a68a0
# OBJ: TLeafElement  Event.Number    Number[Event_] : 0 at: 0x58a6f70
# OBJ: TLeafElement  Event.ReadTime  ReadTime[Event_] : 0 at: 0x58a7670
# OBJ: TLeafElement  Event.ProcTime  ProcTime[Event_] : 0 at: 0x58a7da0
# OBJ: TLeafElement  Event.ProcessID ProcessID[Event_] : 0 at: 0x58a84d0
# OBJ: TLeafElement  Event.MPI   MPI[Event_] : 0 at: 0x58a8bd0
# OBJ: TLeafElement  Event.Weight    Weight[Event_] : 0 at: 0x58a92a0
# OBJ: TLeafElement  Event.Scale Scale[Event_] : 0 at: 0x58a9970
# OBJ: TLeafElement  Event.AlphaQED  AlphaQED[Event_] : 0 at: 0x58aa070
# OBJ: TLeafElement  Event.AlphaQCD  AlphaQCD[Event_] : 0 at: 0x58aa7a0
# OBJ: TLeafElement  Event.ID1   ID1[Event_] : 0 at: 0x58aaea0
# OBJ: TLeafElement  Event.ID2   ID2[Event_] : 0 at: 0x58ab570
# OBJ: TLeafElement  Event.X1    X1[Event_] : 0 at: 0x58abc40
# OBJ: TLeafElement  Event.X2    X2[Event_] : 0 at: 0x58ac310
# OBJ: TLeafElement  Event.ScalePDF  ScalePDF[Event_] : 0 at: 0x58aca10
# OBJ: TLeafElement  Event.PDF1  PDF1[Event_] : 0 at: 0x58ad110
# OBJ: TLeafElement  Event.PDF2  PDF2[Event_] : 0 at: 0x58ad7e0
# OBJ: TLeafI    Event_size  Event_size : 0 at: 0x58adde0
# OBJ: TLeafElement  Weight_ Weight_ : 0 at: 0x58b0500
# OBJ: TLeafElement  Weight.fUniqueID    fUniqueID[Weight_] : 0 at: 0x58b0420
# OBJ: TLeafElement  Weight.fBits    fBits[Weight_] : 0 at: 0x58b0bd0
# OBJ: TLeafElement  Weight.Weight   Weight[Weight_] : 0 at: 0x58b12a0
# OBJ: TLeafI    Weight_size Weight_size : 0 at: 0x58b18a0
# OBJ: TLeafElement  Particle_   Particle_ : 0 at: 0x58b2670
# OBJ: TLeafElement  Particle.fUniqueID  fUniqueID[Particle_] : 0 at: 0x58b2590
# OBJ: TLeafElement  Particle.fBits  fBits[Particle_] : 0 at: 0x58b2d70
# OBJ: TLeafElement  Particle.PID    PID[Particle_] : 0 at: 0x58b3470
# OBJ: TLeafElement  Particle.Status Status[Particle_] : 0 at: 0x58b3b70
# OBJ: TLeafElement  Particle.IsPU   IsPU[Particle_] : 0 at: 0x58b4270
# OBJ: TLeafElement  Particle.M1 M1[Particle_] : 0 at: 0x58b4940
# OBJ: TLeafElement  Particle.M2 M2[Particle_] : 0 at: 0x58b5010
# OBJ: TLeafElement  Particle.D1 D1[Particle_] : 0 at: 0x58b56e0
# OBJ: TLeafElement  Particle.D2 D2[Particle_] : 0 at: 0x58b5db0
# OBJ: TLeafElement  Particle.Charge Charge[Particle_] : 0 at: 0x58b64b0
# OBJ: TLeafElement  Particle.Mass   Mass[Particle_] : 0 at: 0x58b6bb0
# OBJ: TLeafElement  Particle.E  E[Particle_] : 0 at: 0x58b7280
# OBJ: TLeafElement  Particle.Px Px[Particle_] : 0 at: 0x58b7950
# OBJ: TLeafElement  Particle.Py Py[Particle_] : 0 at: 0x58b8020
# OBJ: TLeafElement  Particle.Pz Pz[Particle_] : 0 at: 0x58b86f0
# OBJ: TLeafElement  Particle.P  P[Particle_] : 0 at: 0x58b8dc0
# OBJ: TLeafElement  Particle.PT PT[Particle_] : 0 at: 0x58b9490
# OBJ: TLeafElement  Particle.Eta    Eta[Particle_] : 0 at: 0x58b9b60
# OBJ: TLeafElement  Particle.Phi    Phi[Particle_] : 0 at: 0x58ba230
# OBJ: TLeafElement  Particle.Rapidity   Rapidity[Particle_] : 0 at: 0x58ba960
# OBJ: TLeafElement  Particle.CtgTheta   CtgTheta[Particle_] : 0 at: 0x58bb0f0
# OBJ: TLeafElement  Particle.D0 D0[Particle_] : 0 at: 0x58bb820
# OBJ: TLeafElement  Particle.DZ DZ[Particle_] : 0 at: 0x58bbef0
# OBJ: TLeafElement  Particle.T  T[Particle_] : 0 at: 0x58bc5c0
# OBJ: TLeafElement  Particle.X  X[Particle_] : 0 at: 0x58bcc90
# OBJ: TLeafElement  Particle.Y  Y[Particle_] : 0 at: 0x58bd360
# OBJ: TLeafElement  Particle.Z  Z[Particle_] : 0 at: 0x58bda30
# OBJ: TLeafI    Particle_size   Particle_size : 0 at: 0x58be030
# OBJ: TLeafElement  GenJet_ GenJet_ : 0 at: 0x58bee80
# OBJ: TLeafElement  GenJet.fUniqueID    fUniqueID[GenJet_] : 0 at: 0x58beda0
# OBJ: TLeafElement  GenJet.fBits    fBits[GenJet_] : 0 at: 0x58bf550
# OBJ: TLeafElement  GenJet.PT   PT[GenJet_] : 0 at: 0x58bfbf0
# OBJ: TLeafElement  GenJet.Eta  Eta[GenJet_] : 0 at: 0x58c02c0
# OBJ: TLeafElement  GenJet.Phi  Phi[GenJet_] : 0 at: 0x58c0990
# OBJ: TLeafElement  GenJet.T    T[GenJet_] : 0 at: 0x58c1060
# OBJ: TLeafElement  GenJet.Mass Mass[GenJet_] : 0 at: 0x58c1730
# OBJ: TLeafElement  GenJet.DeltaEta DeltaEta[GenJet_] : 0 at: 0x58c1e30
# OBJ: TLeafElement  GenJet.DeltaPhi DeltaPhi[GenJet_] : 0 at: 0x58c2560
# OBJ: TLeafElement  GenJet.Flavor   Flavor[GenJet_] : 0 at: 0x58c2c60
# OBJ: TLeafElement  GenJet.FlavorAlgo   FlavorAlgo[GenJet_] : 0 at: 0x58c3390
# OBJ: TLeafElement  GenJet.FlavorPhys   FlavorPhys[GenJet_] : 0 at: 0x58c3b20
# OBJ: TLeafElement  GenJet.BTag BTag[GenJet_] : 0 at: 0x58c4250
# OBJ: TLeafElement  GenJet.BTagAlgo BTagAlgo[GenJet_] : 0 at: 0x58c4950
# OBJ: TLeafElement  GenJet.BTagPhys BTagPhys[GenJet_] : 0 at: 0x58c5080
# OBJ: TLeafElement  GenJet.TauTag   TauTag[GenJet_] : 0 at: 0x58c5780
# OBJ: TLeafElement  GenJet.TauWeight    TauWeight[GenJet_] : 0 at: 0x58c5eb0
# OBJ: TLeafElement  GenJet.Charge   Charge[GenJet_] : 0 at: 0x58c65e0
# OBJ: TLeafElement  GenJet.EhadOverEem  EhadOverEem[GenJet_] : 0 at: 0x58c6d10
# OBJ: TLeafElement  GenJet.NCharged NCharged[GenJet_] : 0 at: 0x58c7470
# OBJ: TLeafElement  GenJet.NNeutrals    NNeutrals[GenJet_] : 0 at: 0x58c7bd0
# OBJ: TLeafElement  GenJet.Beta Beta[GenJet_] : 0 at: 0x58c8300
# OBJ: TLeafElement  GenJet.BetaStar BetaStar[GenJet_] : 0 at: 0x58c8a00
# OBJ: TLeafElement  GenJet.MeanSqDeltaR MeanSqDeltaR[GenJet_] : 0 at: 0x58c9160
# OBJ: TLeafElement  GenJet.PTD  PTD[GenJet_] : 0 at: 0x58c9890
# OBJ: TLeafElement  GenJet.FracPt   FracPt[GenJet_] : 0 at: 0x58c9f90
# OBJ: TLeafElement  GenJet.Tau  Tau[GenJet_] : 0 at: 0x58ca660
# OBJ: TLeafElement  GenJet.SoftDroppedJet   SoftDroppedJet[GenJet_] : 0 at: 0x58cad90
# OBJ: TLeafElement  GenJet.SoftDroppedSubJet1   SoftDroppedSubJet1[GenJet_] : 0 at: 0x58cb520
# OBJ: TLeafElement  GenJet.SoftDroppedSubJet2   SoftDroppedSubJet2[GenJet_] : 0 at: 0x58cbcb0
# OBJ: TLeafElement  GenJet.TrimmedP4    TrimmedP4[GenJet_] : 0 at: 0x58cc440
# OBJ: TLeafElement  GenJet.PrunedP4 PrunedP4[GenJet_] : 0 at: 0x58ccbd0
# OBJ: TLeafElement  GenJet.SoftDroppedP4    SoftDroppedP4[GenJet_] : 0 at: 0x58cd330
# OBJ: TLeafElement  GenJet.NSubJetsTrimmed  NSubJetsTrimmed[GenJet_] : 0 at: 0x58cdac0
# OBJ: TLeafElement  GenJet.NSubJetsPruned   NSubJetsPruned[GenJet_] : 0 at: 0x58ce250
# OBJ: TLeafElement  GenJet.NSubJetsSoftDropped  NSubJetsSoftDropped[GenJet_] : 0 at: 0x58ce9e0
# OBJ: TLeafElement  GenJet.ExclYmerge23 ExclYmerge23[GenJet_] : 0 at: 0x58cf170
# OBJ: TLeafElement  GenJet.ExclYmerge34 ExclYmerge34[GenJet_] : 0 at: 0x58cf900
# OBJ: TLeafElement  GenJet.ExclYmerge45 ExclYmerge45[GenJet_] : 0 at: 0x58d0090
# OBJ: TLeafElement  GenJet.ExclYmerge56 ExclYmerge56[GenJet_] : 0 at: 0x58d0820
# OBJ: TLeafElement  GenJet.Constituents Constituents[GenJet_] : 0 at: 0x58d0fb0
# OBJ: TLeafElement  GenJet.Particles    Particles[GenJet_] : 0 at: 0x58d1740
# OBJ: TLeafElement  GenJet.Area Area[GenJet_] : 0 at: 0x58d1e70
# OBJ: TLeafI    GenJet_size GenJet_size : 0 at: 0x58d2470
# OBJ: TLeafElement  GenMissingET_   GenMissingET_ : 0 at: 0x58d31f0
# OBJ: TLeafElement  GenMissingET.fUniqueID  fUniqueID[GenMissingET_] : 0 at: 0x58d3110
# OBJ: TLeafElement  GenMissingET.fBits  fBits[GenMissingET_] : 0 at: 0x58d3920
# OBJ: TLeafElement  GenMissingET.MET    MET[GenMissingET_] : 0 at: 0x58d40b0
# OBJ: TLeafElement  GenMissingET.Eta    Eta[GenMissingET_] : 0 at: 0x58d4840
# OBJ: TLeafElement  GenMissingET.Phi    Phi[GenMissingET_] : 0 at: 0x58d4fd0
# OBJ: TLeafI    GenMissingET_size   GenMissingET_size : 0 at: 0x58d5690
# OBJ: TLeafElement  Jet_    Jet_ : 0 at: 0x58d6480
# OBJ: TLeafElement  Jet.fUniqueID   fUniqueID[Jet_] : 0 at: 0x58d6400
# OBJ: TLeafElement  Jet.fBits   fBits[Jet_] : 0 at: 0x58d6b50
# OBJ: TLeafElement  Jet.PT  PT[Jet_] : 0 at: 0x58d7220
# OBJ: TLeafElement  Jet.Eta Eta[Jet_] : 0 at: 0x58d78f0
# OBJ: TLeafElement  Jet.Phi Phi[Jet_] : 0 at: 0x58d7fc0
# OBJ: TLeafElement  Jet.T   T[Jet_] : 0 at: 0x58d8690
# OBJ: TLeafElement  Jet.Mass    Mass[Jet_] : 0 at: 0x58d8d60
# OBJ: TLeafElement  Jet.DeltaEta    DeltaEta[Jet_] : 0 at: 0x58d9430
# OBJ: TLeafElement  Jet.DeltaPhi    DeltaPhi[Jet_] : 0 at: 0x58d9b00
# OBJ: TLeafElement  Jet.Flavor  Flavor[Jet_] : 0 at: 0x58da1d0
# OBJ: TLeafElement  Jet.FlavorAlgo  FlavorAlgo[Jet_] : 0 at: 0x58da8d0
# OBJ: TLeafElement  Jet.FlavorPhys  FlavorPhys[Jet_] : 0 at: 0x58db000
# OBJ: TLeafElement  Jet.BTag    BTag[Jet_] : 0 at: 0x58db700
# OBJ: TLeafElement  Jet.BTagAlgo    BTagAlgo[Jet_] : 0 at: 0x58dbdd0
# OBJ: TLeafElement  Jet.BTagPhys    BTagPhys[Jet_] : 0 at: 0x58dc4a0
# OBJ: TLeafElement  Jet.TauTag  TauTag[Jet_] : 0 at: 0x58dcb70
# OBJ: TLeafElement  Jet.TauWeight   TauWeight[Jet_] : 0 at: 0x58dd240
# OBJ: TLeafElement  Jet.Charge  Charge[Jet_] : 0 at: 0x58dd910
# OBJ: TLeafElement  Jet.EhadOverEem EhadOverEem[Jet_] : 0 at: 0x58de010
# OBJ: TLeafElement  Jet.NCharged    NCharged[Jet_] : 0 at: 0x58de710
# OBJ: TLeafElement  Jet.NNeutrals   NNeutrals[Jet_] : 0 at: 0x58dede0
# OBJ: TLeafElement  Jet.Beta    Beta[Jet_] : 0 at: 0x58df4b0
# OBJ: TLeafElement  Jet.BetaStar    BetaStar[Jet_] : 0 at: 0x58dfb80
# OBJ: TLeafElement  Jet.MeanSqDeltaR    MeanSqDeltaR[Jet_] : 0 at: 0x58e02b0
# OBJ: TLeafElement  Jet.PTD PTD[Jet_] : 0 at: 0x58e09e0
# OBJ: TLeafElement  Jet.FracPt  FracPt[Jet_] : 0 at: 0x58e10b0
# OBJ: TLeafElement  Jet.Tau Tau[Jet_] : 0 at: 0x58e1780
# OBJ: TLeafElement  Jet.SoftDroppedJet  SoftDroppedJet[Jet_] : 0 at: 0x58e1eb0
# OBJ: TLeafElement  Jet.SoftDroppedSubJet1  SoftDroppedSubJet1[Jet_] : 0 at: 0x58e2640
# OBJ: TLeafElement  Jet.SoftDroppedSubJet2  SoftDroppedSubJet2[Jet_] : 0 at: 0x58e2dd0
# OBJ: TLeafElement  Jet.TrimmedP4   TrimmedP4[Jet_] : 0 at: 0x58e3530
# OBJ: TLeafElement  Jet.PrunedP4    PrunedP4[Jet_] : 0 at: 0x58e3c00
# OBJ: TLeafElement  Jet.SoftDroppedP4   SoftDroppedP4[Jet_] : 0 at: 0x58e4330
# OBJ: TLeafElement  Jet.NSubJetsTrimmed NSubJetsTrimmed[Jet_] : 0 at: 0x58e4ac0
# OBJ: TLeafElement  Jet.NSubJetsPruned  NSubJetsPruned[Jet_] : 0 at: 0x58e5250
# OBJ: TLeafElement  Jet.NSubJetsSoftDropped NSubJetsSoftDropped[Jet_] : 0 at: 0x58e59e0
# OBJ: TLeafElement  Jet.ExclYmerge23    ExclYmerge23[Jet_] : 0 at: 0x58e6170
# OBJ: TLeafElement  Jet.ExclYmerge34    ExclYmerge34[Jet_] : 0 at: 0x58e6900
# OBJ: TLeafElement  Jet.ExclYmerge45    ExclYmerge45[Jet_] : 0 at: 0x58e7090
# OBJ: TLeafElement  Jet.ExclYmerge56    ExclYmerge56[Jet_] : 0 at: 0x58e7820
# OBJ: TLeafElement  Jet.Constituents    Constituents[Jet_] : 0 at: 0x58e7fb0
# OBJ: TLeafElement  Jet.Particles   Particles[Jet_] : 0 at: 0x58e86e0
# OBJ: TLeafElement  Jet.Area    Area[Jet_] : 0 at: 0x58e8db0
# OBJ: TLeafI    Jet_size    Jet_size : 0 at: 0x58e93b0
# OBJ: TLeafElement  Electron_   Electron_ : 0 at: 0x58ea130
# OBJ: TLeafElement  Electron.fUniqueID  fUniqueID[Electron_] : 0 at: 0x58ea050
# OBJ: TLeafElement  Electron.fBits  fBits[Electron_] : 0 at: 0x58ea830
# OBJ: TLeafElement  Electron.PT PT[Electron_] : 0 at: 0x58eaf30
# OBJ: TLeafElement  Electron.Eta    Eta[Electron_] : 0 at: 0x58eb600
# OBJ: TLeafElement  Electron.Phi    Phi[Electron_] : 0 at: 0x58ebcd0
# OBJ: TLeafElement  Electron.T  T[Electron_] : 0 at: 0x58ec3a0
# OBJ: TLeafElement  Electron.Charge Charge[Electron_] : 0 at: 0x58ecaa0
# OBJ: TLeafElement  Electron.EhadOverEem    EhadOverEem[Electron_] : 0 at: 0x58ed200
# OBJ: TLeafElement  Electron.Particle   Particle[Electron_] : 0 at: 0x58ed990
# OBJ: TLeafElement  Electron.IsolationVar   IsolationVar[Electron_] : 0 at: 0x58ee120
# OBJ: TLeafElement  Electron.IsolationVarRhoCorr    IsolationVarRhoCorr[Electron_] : 0 at: 0x58ee8b0
# OBJ: TLeafElement  Electron.SumPtCharged   SumPtCharged[Electron_] : 0 at: 0x58ef040
# OBJ: TLeafElement  Electron.SumPtNeutral   SumPtNeutral[Electron_] : 0 at: 0x58ef7d0
# OBJ: TLeafElement  Electron.SumPtChargedPU SumPtChargedPU[Electron_] : 0 at: 0x58eff60
# OBJ: TLeafElement  Electron.SumPt  SumPt[Electron_] : 0 at: 0x58f06c0
# OBJ: TLeafI    Electron_size   Electron_size : 0 at: 0x58f0cf0
# OBJ: TLeafElement  Photon_ Photon_ : 0 at: 0x58f1a70
# OBJ: TLeafElement  Photon.fUniqueID    fUniqueID[Photon_] : 0 at: 0x58f1990
# OBJ: TLeafElement  Photon.fBits    fBits[Photon_] : 0 at: 0x58f2140
# OBJ: TLeafElement  Photon.PT   PT[Photon_] : 0 at: 0x58f2810
# OBJ: TLeafElement  Photon.Eta  Eta[Photon_] : 0 at: 0x58f2ee0
# OBJ: TLeafElement  Photon.Phi  Phi[Photon_] : 0 at: 0x58f35b0
# OBJ: TLeafElement  Photon.E    E[Photon_] : 0 at: 0x58f3c80
# OBJ: TLeafElement  Photon.T    T[Photon_] : 0 at: 0x58f4350
# OBJ: TLeafElement  Photon.EhadOverEem  EhadOverEem[Photon_] : 0 at: 0x58f4a80
# OBJ: TLeafElement  Photon.Particles    Particles[Photon_] : 0 at: 0x58f5210
# OBJ: TLeafElement  Photon.IsolationVar IsolationVar[Photon_] : 0 at: 0x58f59a0
# OBJ: TLeafElement  Photon.IsolationVarRhoCorr  IsolationVarRhoCorr[Photon_] : 0 at: 0x58f6130
# OBJ: TLeafElement  Photon.SumPtCharged SumPtCharged[Photon_] : 0 at: 0x58f68c0
# OBJ: TLeafElement  Photon.SumPtNeutral SumPtNeutral[Photon_] : 0 at: 0x58f7050
# OBJ: TLeafElement  Photon.SumPtChargedPU   SumPtChargedPU[Photon_] : 0 at: 0x58f77e0
# OBJ: TLeafElement  Photon.SumPt    SumPt[Photon_] : 0 at: 0x58f7f10
# OBJ: TLeafElement  Photon.Status   Status[Photon_] : 0 at: 0x58f85e0
# OBJ: TLeafI    Photon_size Photon_size : 0 at: 0x3f30fc0
# OBJ: TLeafElement  Muon_   Muon_ : 0 at: 0x3f31ce0
# OBJ: TLeafElement  Muon.fUniqueID  fUniqueID[Muon_] : 0 at: 0x3f31c30
# OBJ: TLeafElement  Muon.fBits  fBits[Muon_] : 0 at: 0x3f323b0
# OBJ: TLeafElement  Muon.PT PT[Muon_] : 0 at: 0x3f32a80
# OBJ: TLeafElement  Muon.Eta    Eta[Muon_] : 0 at: 0x3f33150
# OBJ: TLeafElement  Muon.Phi    Phi[Muon_] : 0 at: 0x3f33820
# OBJ: TLeafElement  Muon.T  T[Muon_] : 0 at: 0x3f33ef0
# OBJ: TLeafElement  Muon.Charge Charge[Muon_] : 0 at: 0x3f345c0
# OBJ: TLeafElement  Muon.Particle   Particle[Muon_] : 0 at: 0x3f34c90
# OBJ: TLeafElement  Muon.IsolationVar   IsolationVar[Muon_] : 0 at: 0x3f353c0
# OBJ: TLeafElement  Muon.IsolationVarRhoCorr    IsolationVarRhoCorr[Muon_] : 0 at: 0x3f35b50
# OBJ: TLeafElement  Muon.SumPtCharged   SumPtCharged[Muon_] : 0 at: 0x3f362e0
# OBJ: TLeafElement  Muon.SumPtNeutral   SumPtNeutral[Muon_] : 0 at: 0x3f36a70
# OBJ: TLeafElement  Muon.SumPtChargedPU SumPtChargedPU[Muon_] : 0 at: 0x5904d00
# OBJ: TLeafElement  Muon.SumPt  SumPt[Muon_] : 0 at: 0x5905430
# OBJ: TLeafI    Muon_size   Muon_size : 0 at: 0x5905a30
# OBJ: TLeafElement  MissingET_  MissingET_ : 0 at: 0x59067b0
# OBJ: TLeafElement  MissingET.fUniqueID fUniqueID[MissingET_] : 0 at: 0x59066d0
# OBJ: TLeafElement  MissingET.fBits fBits[MissingET_] : 0 at: 0x5906eb0
# OBJ: TLeafElement  MissingET.MET   MET[MissingET_] : 0 at: 0x59075b0
# OBJ: TLeafElement  MissingET.Eta   Eta[MissingET_] : 0 at: 0x5907c80
# OBJ: TLeafElement  MissingET.Phi   Phi[MissingET_] : 0 at: 0x5908350
# OBJ: TLeafI    MissingET_size  MissingET_size : 0 at: 0x5908980
# OBJ: TLeafElement  ScalarHT_   ScalarHT_ : 0 at: 0x5909700
# OBJ: TLeafElement  ScalarHT.fUniqueID  fUniqueID[ScalarHT_] : 0 at: 0x5909620
# OBJ: TLeafElement  ScalarHT.fBits  fBits[ScalarHT_] : 0 at: 0x5909e00
# OBJ: TLeafElement  ScalarHT.HT HT[ScalarHT_] : 0 at: 0x590a500
# OBJ: TLeafI    ScalarHT_size   ScalarHT_size : 0 at: 0x590ab00
# OBJ: TLeafElement  Rho_    Rho_ : 0 at: 0x590b7c0
# OBJ: TLeafElement  Rho.fUniqueID   fUniqueID[Rho_] : 0 at: 0x590b740
# OBJ: TLeafElement  Rho.fBits   fBits[Rho_] : 0 at: 0x590be90
# OBJ: TLeafElement  Rho.Rho Rho[Rho_] : 0 at: 0x590c560
# OBJ: TLeafElement  Rho.Edges   Edges[Rho_] : 0 at: 0x590cc30
# OBJ: TLeafI    Rho_size    Rho_size : 0 at: 0x590d230
# OBJ: TLeafElement  Vertex_ Vertex_ : 0 at: 0x590dfc0
# OBJ: TLeafElement  Vertex.fUniqueID    fUniqueID[Vertex_] : 0 at: 0x590dee0
# OBJ: TLeafElement  Vertex.fBits    fBits[Vertex_] : 0 at: 0x590e690
# OBJ: TLeafElement  Vertex.T    T[Vertex_] : 0 at: 0x590ed60
# OBJ: TLeafElement  Vertex.X    X[Vertex_] : 0 at: 0x590f430
# OBJ: TLeafElement  Vertex.Y    Y[Vertex_] : 0 at: 0x590fb00
# OBJ: TLeafElement  Vertex.Z    Z[Vertex_] : 0 at: 0x59101d0
# OBJ: TLeafElement  Vertex.ErrorT   ErrorT[Vertex_] : 0 at: 0x59108a0
# OBJ: TLeafElement  Vertex.ErrorX   ErrorX[Vertex_] : 0 at: 0x5910f70
# OBJ: TLeafElement  Vertex.ErrorY   ErrorY[Vertex_] : 0 at: 0x5911640
# OBJ: TLeafElement  Vertex.ErrorZ   ErrorZ[Vertex_] : 0 at: 0x5911d10
# OBJ: TLeafElement  Vertex.Index    Index[Vertex_] : 0 at: 0x59123e0
# OBJ: TLeafElement  Vertex.NDF  NDF[Vertex_] : 0 at: 0x5912ab0
# OBJ: TLeafElement  Vertex.Sigma    Sigma[Vertex_] : 0 at: 0x5913180
# OBJ: TLeafElement  Vertex.SumPT2   SumPT2[Vertex_] : 0 at: 0x5913850
# OBJ: TLeafElement  Vertex.GenSumPT2    GenSumPT2[Vertex_] : 0 at: 0x5913f80
# OBJ: TLeafElement  Vertex.GenDeltaZ    GenDeltaZ[Vertex_] : 0 at: 0x5914710
# OBJ: TLeafElement  Vertex.BTVSumPT2    BTVSumPT2[Vertex_] : 0 at: 0x5914ea0
# OBJ: TLeafElement  Vertex.Constituents Constituents[Vertex_] : 0 at: 0x5915630
# OBJ: TLeafI    Vertex_size Vertex_size : 0 at: 0x5915c90

