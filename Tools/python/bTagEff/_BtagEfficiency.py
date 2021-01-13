#''' Implementation of b-tagging reweighting
#'''
#
## TODO:
## add read pt borders
## add scalefactor calculation
## check method 1a
#
#from operator import mul
#
## Logging
#import logging
#logger = logging.getLogger(__name__)
#
## binning in pt
## taken from https://twiki.cern.ch/twiki/pub/CMS/YR2018Systematics/btag_syst_Phase2.txt
#ptinf = [ 20, 30, 50,  70, 100, 140, 200, 300,  600, 1000, -1 ]
##ptsup = [ 30, 50, 70, 100, 140, 200, 300, 600, 1000, 3000 ]
#
#syst_B_loose   = [ 0.015, 0.01, 0.01, 0.01, 0.01, 0.01, 0.012, 0.018, 0.021, 0.042 ]
#syst_B_medium  = [ 0.017, 0.01, 0.01, 0.01, 0.01, 0.01, 0.016, 0.018, 0.023, 0.046 ]
#syst_B_tight   = [ 0.018, 0.01, 0.01, 0.01, 0.01, 0.01, 0.016, 0.018, 0.029, 0.058 ]
#syst_B_default = syst_B_medium
#
#syst_L_loose   = [ 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05 ]
#syst_L_medium  = [ 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10 ]
#syst_L_tight   = [ 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15, 0.15 ]
#syst_L_default = syst_L_medium
#
#syst_B = [ (pt, ptinf[i+1], syst_B_default[i]) for i, pt in enumerate(ptinf[:-1]) ]
#syst_L = [ (pt, ptinf[i+1], syst_L_default[i]) for i, pt in enumerate(ptinf[:-1]) ]
#
#class BtagEfficiency:
#
#    def getBTagSF_1a( self, var, bJets, nonBJets ):
#
#        if var not in self.btagWeightNames:
#            raise ValueError( "Don't know what to do with b-tag variation %s" %var )
#
#        # wrong calculation!?!
#        # return reduce( mul, [ j['jetSF'][var] for j in bJets ] + [ 1-j['jetSF'][var] for j in nonBJets ], 1 )
#        return reduce( mul, [ j['jetSF'][var] for j in bJets ] + [ 1-j['jetSF'][var] for j in nonBJets ], 1 )
#
#
#    def __init__( self ):
#
#        # All btag weight names per jet
#        self.btagWeightNames = ['SF', 'SF_b_Down', 'SF_b_Up', 'SF_l_Down', 'SF_l_Up']
#
#    def addJetTagEffToJet( self, j, btagWP ):
#        # BTag SF Not implemented below 20 GeV
#        if j['pt'] < ptinf[0]: 
#            j['jetSF'] =  { sf:1 for sf in self.btagWeightNames }
#        
#        j['jetSF'] = {}
#
#        if j['bTag_'+btagWP]:
#            #SF for b
#            j['jetSF']['SF']         = self.getB_SF( j['pt'] )
#            j['jetSF']['SF_b_Down']  = self.getB_SF( j['pt_JEC_down'] )
#            j['jetSF']['SF_b_Up']    = self.getB_SF( j['pt_JEC_up'] )
#            j['jetSF']['SF_l_Down']  = 1.
#            j['jetSF']['SF_l_Up']    = 1.
#        else:
#            #SF for light flavours
#            j['jetSF']['SF']         = self.getL_SF( j['pt'] )
#            j['jetSF']['SF_b_Down']  = 1.
#            j['jetSF']['SF_b_Up']    = 1.
#            j['jetSF']['SF_l_Down']  = self.getL_SF( j['pt_JEC_down'] )
#            j['jetSF']['SF_l_Up']    = self.getL_SF( j['pt_JEC_up'] )
#
#
#    def getB_SF( self, pt ):
#        ''' Get SF for b jet
#        '''
#        for bin in syst_B:
#            if pt>=bin[0] and (pt<bin[1] or bin[1]<0):
#                return 1 + bin[2]                
#
#        logger.debug( "No SF for pt %f" %pt )
#        return 1
#
#    def getL_SF( self, pt ):
#        ''' Get SF for light jet
#        '''
#        for bin in syst_L:
#            if pt>=bin[0] and (pt<bin[1] or bin[1]<0):
#                return 1 + bin[2]                
#
#        logger.debug( "No SF for pt %f" %pt )
#        return 1
#
#
