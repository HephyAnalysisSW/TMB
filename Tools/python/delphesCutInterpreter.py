''' Class to interpret string based cuts
'''

import logging
logger = logging.getLogger(__name__)

jetSelection    = "nJetGood"
bJetSelectionM  = "nBTag"

special_cuts = {
    "photon" :   "(recoPhoton_pt[0]>20)",
    "singlelep": "(nrecoLep>=1)&&(recoLep_pt[0]>25)&&Sum$(recoLep_pt>25)<=1",
    "dilep" :    "(nrecoLep>=2)&&(recoLep_pt[0]>25)&&Sum$(recoLep_pt>15)>=2&&Sum$(recoLep_pt>25)<=2",
    "WHJet" :    "Sum$(recoJet_pt>25&&abs(recoJet_eta)<2.4&&recoJet_bTag>=1)>=2&&Sum$(recoJet_pt>25&&abs(recoJet_eta)<2.4)<=3",
    "realW" :    "WH_nu_E>0",
    #"protectWH" :    "sqrt(acos(cos(H_phi-WH_W_phi)**2+(H_eta-WH_W_eta)**2)<5.5",
    "ZHJet" :    "Sum$(recoJet_pt>25&&abs(recoJet_eta)<2.4&&recoJet_bTag>=1)>=2&&Sum$(recoJet_pt>25&&abs(recoJet_eta)<2.4)<=3",
    "trilep":    "(nrecoLep>=3)&&(recoLep_pt[0]>40&&recoLep_pt[1]>20&&recoLep_pt[2]>10)",
    "onZ"   :    "abs(recoZ_mass-91.2)<10",
    "onH"   :    "(H_dijet_mass>90&&H_dijet_mass<150)",
    "LepGBB":    "(sqrt(acos(cos(recoLep_phi[0]-recoPhoton_phi[0]))**2+(recoLep_eta[0]-recoPhoton_eta[0])**2)>2.5)",
    "odd":       "(evt%2==1)",
    "even":      "(evt%2==0)",
  }

continous_variables = [ ("ptG", "recoPhoton_pt[0]"), ("met", "recoMet_pt"), ("ptZ", "recoZ_pt"), ("ptW", "WH_W_pt")]
discrete_variables  = [ ("njet", "nrecoJet"), ("btag", "nBTag") ]

class cutInterpreter:
    ''' Translate var100to200-var2p etc.
    '''

    @staticmethod
    def translate_cut_to_string( string ):

        if string.startswith("multiIso"):
            str_ = mIsoWP[ string.replace('multiIso','') ]
            return "l1_mIsoWP>%i&&l2_mIsoWP>%i" % (str_, str_)
        elif string.startswith("relIso"):
           iso = float( string.replace('relIso','') )
           raise ValueError("We do not want to use relIso for our analysis anymore!")
           return "l1_relIso03<%3.2f&&l2_relIso03<%3.2f"%( iso, iso )
        elif string.startswith("miniIso"):
           iso = float( string.replace('miniIso','') )
           return "l1_miniRelIso<%3.2f&&l2_miniRelIso<%3.2f"%( iso, iso )
        # special cuts
        if string in special_cuts.keys(): return special_cuts[string]

        # continous Variables and discrete variables with "To"
        for var, tree_var in continous_variables + discrete_variables:
            isDiscrete = (var, tree_var) in discrete_variables
            if string.startswith( var ):
                num_str = string[len( var ):].replace("to","To").split("To")
                # don't do discrete variables without "To"
                if isDiscrete and len(num_str)<=1:
                    continue
                upper = None
                lower = None
                if len(num_str)==2:
                    lower, upper = num_str
                elif len(num_str)==1:
                    lower = num_str[0]
                else:
                    raise ValueError( "Can't interpret string %s" % string )
                res_string = []
                if lower: res_string.append( tree_var+">="+lower )
                leq = "<=" if isDiscrete else "<"
                if upper: res_string.append( tree_var+leq+upper )
                return "&&".join( res_string )

        # discrete Variables
        for var, tree_var in discrete_variables:
            logger.debug("Reading discrete cut %s as %s"%(var, tree_var))
            if string.startswith( var ):
                # Omit discrete variables, done above
                if string[len( var ):].replace("to","To").count("To"):
                    continue
                    #raise NotImplementedError( "Can't interpret string with 'to' for discrete variable: %s. You just volunteered." % string )

                num_str = string[len( var ):]
                # logger.debug("Num string is %s"%(num_str))
                # var1p -> tree_var >= 1
                if num_str[-1] == 'p' and len(num_str)==2:
                    # logger.debug("Using cut string %s"%(tree_var+">="+num_str[0]))
                    return tree_var+">="+num_str[0]
                # var123->tree_var==1||tree_var==2||tree_var==3
                else:
                    vls = [ tree_var+"=="+c for c in num_str ]
                    if len(vls)==1:
                      # logger.debug("Using cut string %s"%vls[0])
                      return vls[0]
                    else:
                      # logger.debug("Using cut string %s"%'('+'||'.join(vls)+')')
                      return '('+'||'.join(vls)+')'
        raise ValueError( "Can't interpret string %s. All cuts %s" % (string,  ", ".join( [ c[0] for c in continous_variables + discrete_variables] +  special_cuts.keys() ) ) )

    @staticmethod
    def cutString( cut, select = [""], ignore = []):
        ''' Cutstring syntax: cut1-cut2-cut3
        '''
        if cut is None or cut=="":
            return "(1)"
        cuts = cut.split('-')
        # require selected
        cuts = filter( lambda c: any( sel in c for sel in select ), cuts )
        # ignore
        cuts = filter( lambda c: not any( ign in c for ign in ignore ), cuts )

        cutString = "&&".join( map( cutInterpreter.translate_cut_to_string, cuts ) )

        return cutString
    
    @staticmethod
    def cutList ( cut, select = [""], ignore = []):
        ''' Cutstring syntax: cut1-cut2-cut3
        '''
        cuts = cut.split('-')
        # require selected
        cuts = filter( lambda c: any( sel in c for sel in select ), cuts )
        # ignore
        cuts = filter( lambda c: not any( ign in c for ign in ignore ), cuts )
        return [ cutInterpreter.translate_cut_to_string(cut) for cut in cuts ] 
        #return  "&&".join( map( cutInterpreter.translate_cut_to_string, cuts ) )

if __name__ == "__main__":
    print cutInterpreter.cutString("njet2-btag0p-multiIsoVT-relIso0.12-looseLeptonVeto-mll20-onZ-met80-metSig5-dPhiJet0-dPhiJet1-mt2ll100")
    print
    print cutInterpreter.cutList("njet2-btag0p-multiIsoVT-relIso0.12-looseLeptonVeto-mll20-onZ-met80-metSig5-dPhiJet0-dPhiJet1-mt2ll100")
