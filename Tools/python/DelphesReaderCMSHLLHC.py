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
        res = self.read_collection( 'MuonLoose', 
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
        return self.read_collection( 'JetPUPPI', 
            [   ('PT', 'pt'), ( 'Eta', 'eta'), ('Phi', 'phi'),
                ('BTag', 'bTag'), ( 'BTagPhys', 'bTagPhys'), ('Flavor', 'flavor'),
                ('NCharged', 'nCharged'), ('NNeutrals', 'nNeutrals'), 
            ])

    def photons( self ):
        return self.read_collection( 'PhotonTight', 
            [   ('PT', 'pt'), ( 'Eta', 'eta'), ('Phi', 'phi'),
                ('IsolationVar', 'isolationVar'), ('IsolationVarRhoCorr', 'isolationVarRhoCorr'),  
                ('SumPtCharged', 'sumPtCharged'),  ('SumPtNeutral', 'sumPtNeutral'), ('SumPtChargedPU', 'sumPtChargedPU'),  ('SumPt', 'sumPt'),
                ('EhadOverEem','ehadOverEem') 
            ])

    def genJets( self ):
        return self.read_collection( 'GenJet', 
            [   ('PT', 'pt'), ( 'Eta', 'eta'), ('Phi', 'phi'),
            ])

    def met( self ):
        return self.read_collection( 'PuppiMissingET', [('MET', 'pt'), ('Phi', 'phi')] )

#    def getEntry( self, entry ):
#        # Set the Suchi-Reader 
#        self.reader.GetEntry( entry )

# OBJ: TLeafElement  JetPUPPI_   JetPUPPI_ : 0 at: 0x5296ed0
# OBJ: TLeafElement  JetPUPPI.fUniqueID  fUniqueID[JetPUPPI_] : 0 at: 0x5296df0
# OBJ: TLeafElement  JetPUPPI.fBits  fBits[JetPUPPI_] : 0 at: 0x52976b0
# OBJ: TLeafElement  JetPUPPI.PT PT[JetPUPPI_] : 0 at: 0x5297e90
# OBJ: TLeafElement  JetPUPPI.Eta    Eta[JetPUPPI_] : 0 at: 0x5298640
# OBJ: TLeafElement  JetPUPPI.Phi    Phi[JetPUPPI_] : 0 at: 0x5298df0
# OBJ: TLeafElement  JetPUPPI.T  T[JetPUPPI_] : 0 at: 0x52995a0
# OBJ: TLeafElement  JetPUPPI.Mass   Mass[JetPUPPI_] : 0 at: 0x5299d70
# OBJ: TLeafElement  JetPUPPI.DeltaEta   DeltaEta[JetPUPPI_] : 0 at: 0x529a5a0
# OBJ: TLeafElement  JetPUPPI.DeltaPhi   DeltaPhi[JetPUPPI_] : 0 at: 0x529ae10
# OBJ: TLeafElement  JetPUPPI.Flavor Flavor[JetPUPPI_] : 0 at: 0x529b670
# OBJ: TLeafElement  JetPUPPI.FlavorAlgo FlavorAlgo[JetPUPPI_] : 0 at: 0x529bed0
# OBJ: TLeafElement  JetPUPPI.FlavorPhys FlavorPhys[JetPUPPI_] : 0 at: 0x529c740
# OBJ: TLeafElement  JetPUPPI.BTag   BTag[JetPUPPI_] : 0 at: 0x529cf70
# OBJ: TLeafElement  JetPUPPI.BTagAlgo   BTagAlgo[JetPUPPI_] : 0 at: 0x529d7a0
# OBJ: TLeafElement  JetPUPPI.BTagPhys   BTagPhys[JetPUPPI_] : 0 at: 0x529e010
# OBJ: TLeafElement  JetPUPPI.TauTag TauTag[JetPUPPI_] : 0 at: 0x529e870
# OBJ: TLeafElement  JetPUPPI.TauWeight  TauWeight[JetPUPPI_] : 0 at: 0x529f0d0
# OBJ: TLeafElement  JetPUPPI.Charge Charge[JetPUPPI_] : 0 at: 0x529f930
# OBJ: TLeafElement  JetPUPPI.EhadOverEem    EhadOverEem[JetPUPPI_] : 0 at: 0x52a0190
# OBJ: TLeafElement  JetPUPPI.NCharged   NCharged[JetPUPPI_] : 0 at: 0x52a0a00
# OBJ: TLeafElement  JetPUPPI.NNeutrals  NNeutrals[JetPUPPI_] : 0 at: 0x52a1270
# OBJ: TLeafElement  JetPUPPI.Beta   Beta[JetPUPPI_] : 0 at: 0x52a1aa0
# OBJ: TLeafElement  JetPUPPI.BetaStar   BetaStar[JetPUPPI_] : 0 at: 0x52a22d0
# OBJ: TLeafElement  JetPUPPI.MeanSqDeltaR   MeanSqDeltaR[JetPUPPI_] : 0 at: 0x52a2b40
# OBJ: TLeafElement  JetPUPPI.PTD    PTD[JetPUPPI_] : 0 at: 0x523be00
# OBJ: TLeafElement  JetPUPPI.FracPt FracPt[JetPUPPI_] : 0 at: 0x523c610
# OBJ: TLeafElement  JetPUPPI.Tau    Tau[JetPUPPI_] : 0 at: 0x523ce30
# OBJ: TLeafElement  JetPUPPI.SoftDroppedJet SoftDroppedJet[JetPUPPI_] : 0 at: 0x523d640
# OBJ: TLeafElement  JetPUPPI.SoftDroppedSubJet1 SoftDroppedSubJet1[JetPUPPI_] : 0 at: 0x523def0
# OBJ: TLeafElement  JetPUPPI.SoftDroppedSubJet2 SoftDroppedSubJet2[JetPUPPI_] : 0 at: 0x523e7a0
# OBJ: TLeafElement  JetPUPPI.TrimmedP4  TrimmedP4[JetPUPPI_] : 0 at: 0x523f050
# OBJ: TLeafElement  JetPUPPI.PrunedP4   PrunedP4[JetPUPPI_] : 0 at: 0x523fa10
# OBJ: TLeafElement  JetPUPPI.SoftDroppedP4  SoftDroppedP4[JetPUPPI_] : 0 at: 0x52403d0
# OBJ: TLeafElement  JetPUPPI.NSubJetsTrimmed    NSubJetsTrimmed[JetPUPPI_] : 0 at: 0x5240d90
# OBJ: TLeafElement  JetPUPPI.NSubJetsPruned NSubJetsPruned[JetPUPPI_] : 0 at: 0x5241600
# OBJ: TLeafElement  JetPUPPI.NSubJetsSoftDropped    NSubJetsSoftDropped[JetPUPPI_] : 0 at: 0x5241e70
# OBJ: TLeafElement  JetPUPPI.ExclYmerge23   ExclYmerge23[JetPUPPI_] : 0 at: 0x52426e0
# OBJ: TLeafElement  JetPUPPI.ExclYmerge34   ExclYmerge34[JetPUPPI_] : 0 at: 0x5242f50
# OBJ: TLeafElement  JetPUPPI.ExclYmerge45   ExclYmerge45[JetPUPPI_] : 0 at: 0x52437c0
# OBJ: TLeafElement  JetPUPPI.ExclYmerge56   ExclYmerge56[JetPUPPI_] : 0 at: 0x5244030
# OBJ: TLeafElement  JetPUPPI.Constituents   Constituents[JetPUPPI_] : 0 at: 0x52448a0
# OBJ: TLeafElement  JetPUPPI.Particles  Particles[JetPUPPI_] : 0 at: 0x5245150
# OBJ: TLeafElement  JetPUPPI.Area   Area[JetPUPPI_] : 0 at: 0x52459f0

