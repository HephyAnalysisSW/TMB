# Standard Imports
import textwrap
from Analysis.Tools.helpers import getVarValue, getObjDict

photonIdCutBasedBitmap = {                     'loose':1, 'medium':2, 'tight':4 }  # NanoAOD Version ID bitmap, 2^(0:loose, 1:medium, 2:tight)
photonIdCutBased       = { 'fail':0,           'loose':1, 'medium':2, 'tight':3 }  # NanoAOD Version
electronIdCutBased     = { 'fail':0, 'veto':1, 'loose':2, 'medium':3, 'tight':4 }  # NanoAOD Version

muonPfIsoId            = { 'PFIsoVeryLoose':1, 'PFIsoLoose':2, 'PFIsoMedium':3, 'PFIsoTight':4, 'PFIsoVeryTight':5, 'PFIsoVeryVeryTight':6 }

jetIdNamingList        = [ "tightLepVeto", "tight", "loose" ]

# see https://cms-nanoaod-integration.web.cern.ch/integration/master-102X/mc80X_doc.html
# or  https://cms-nanoaod-integration.web.cern.ch/integration/master-102X/mc102X_doc.html
# Attention: only for nanoAOD v94x or higher (in 80x, only 2 bits are used)
vidNestedWPBitMapNamingList = [ 'MinPtCut', 'GsfEleSCEtaMultiRangeCut', 'GsfEleDEtaInSeedCut', 'GsfEleDPhiInCut', 'GsfEleFull5x5SigmaIEtaIEtaCut', 'GsfEleHadronicOverEMEnergyScaledCut', 'GsfEleEInverseMinusPInverseCut', 'GsfEleRelPFIsoScaledCut', 'GsfEleConversionVetoCut', 'GsfEleMissingHitsCut' ][::-1]
vidNestedWPBitMap           = { 'fail':0, 'veto':1, 'loose':2, 'medium':3, 'tight':4 }  # Bitwise (Electron vidNestedWPBitMap ID flags (3 bits per cut), '000'=0 is fail, '001'=1 is veto, '010'=2 is loose, '011'=3 is medium, '100'=4 is tight)

# the PhoAnyPFIsoWithEACut is used twice in the root info, thus I call one PhoAnyPFIsoWithEACut2 (they don't have the same cuts)
# seems like the list is implemented reversed, thus the [::-1] at the end (values make more sense in that way)
vidNestedWPBitMapNamingListPhoton = [ 'MinPtCut', 'PhoSCEtaMultiRangeCut', 'PhoSingleTowerHadOverEmCut', 'PhoFull5x5SigmaIEtaIEtaCut', 'PhoAnyPFIsoWithEACut', 'PhoAnyPFIsoWithEAAndQuadScalingCut', 'PhoAnyPFIsoWithEACut2' ][::-1] # HoverE (2), SIEIE (3), chgIso (4), NeuIso (5), PhoIso (6) cuts



# Attention: only for nanoAOD v94x or higher (in 80x, only 2 bits are used)
def jetIdBitMapToDict( val ):
    # convert int of vidNestedWPBitMap ( e.g. val = 6 ) to bitmap ( e.g. "110"), then to list ( e.g. [ 1, 1, 0 ] )
    # Jet ID flags bit1 is loose (always false in 2017 since it does not exist), bit2 is tight, bit3 is tightLepVeto : 0 at: 0x5f1a030
    # in 80x no tightLepVeto exists, only 2 bits, bit3 is set to 0
    # create dictionary
    # however some entries of jets don't have a jedId entry, no clue why (very limited number of events e.g. in DY M50 Fall2017)
    try: idList = map( lambda x: int(x), "{0:03b}".format( val ) )
    except: idList = [ 0, 0, 0 ]
    return dict( zip( jetIdNamingList, idList ) )

def photonIdBitMapToDict( val ):
    idList = [ int( x, 2 ) for x in textwrap.wrap( "{0:014b}".format( val ) , 2) ]
    return dict( zip( vidNestedWPBitMapNamingListPhoton, idList ) )

# Attention: only for nanoAOD v94x or higher (in 80x, only 2 bits are used)
def vidNestedWPBitMapToDict( val ):
    # convert int of vidNestedWPBitMap ( e.g. val = 611099940 ) to bitmap ( e.g. "100100011011001010010100100100")
    # split vidBitmap string (containing 3 bits per cut) in parts of 3 bits ( e.g. ["100","100","011","011","001","010","010","100","100","100"] )
    # convert 3 bits to int ( e.g. [4, 4, 3, 3, 1, 2, 2, 4, 4, 4])
    # create dictionary
    idList = [ int( x, 2 ) for x in textwrap.wrap( "{0:030b}".format( val ) , 3) ] #use 2 for nanoAOD version 80x
    return dict( zip( vidNestedWPBitMapNamingList, idList ) )

# General Selection Functions
def deltaRCleaning( cleaningParticles, otherParticles, dRCut = 0.4 ):

    from Analysis.Tools.helpers  import deltaR

    res = []
    for part in cleaningParticles:
        clean = True
        for otherParticle in otherParticles:
            if deltaR( otherParticle, part ) < dRCut:
                clean = False
                break
        if clean:
            res.append( part )
    res.sort( key = lambda p: -p['pt'] )
    return res

def getParticles( c, collVars, coll ):
    from Analysis.Tools.helpers import getVarValue, getObjDict
    return [ getObjDict( c, coll+'_', collVars, i ) for i in range(int(getVarValue(c, 'n'+coll))) ]

# Reco b-Jet Filter
def isBJet( j, tagger='DeepCSV', year=2016 ):
    if tagger == 'CSVv2':
        if year == 2016:
            # https://twiki.cern.ch/twikix/bin/viewauth/CMS/BtagRecommendation80XReReco
            return j['btagCSVV2'] > 0.8484 
        elif year == 2017:
            # https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation94X
            return j['btagCSVV2'] > 0.8838 
        elif year == 2018:
            # UPDATE WHEN AVAILABLE
            return j['btagCSVV2'] > 0.8838 
        else:
            raise (NotImplementedError, "Don't know what cut to use for year %s"%year)
    elif tagger == 'DeepCSV':
        if year == 2016:
            # https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation2016Legacy
            return j['btagDeepB'] > 0.6321
        elif year == 2017:
            # https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation94X
            return j['btagDeepB'] > 0.4941
        elif year == 2018:
            # https://twiki.cern.ch/twiki/bin/viewauth/CMS/BtagRecommendation102X
            return j['btagDeepB'] > 0.4184
        else:
            raise (NotImplementedError, "Don't know what cut to use for year %s"%year)

def vertexSelector( l, barrelECdiff=False ):
    if abs(l['pdgId']) == 11: absEta = abs(l["eta"] + l["deltaEtaSC"])   # eta supercluster
    else:                     absEta = abs(l["eta"])                     # eta
    EC = absEta > 1.479 and barrelECdiff # only if difference for EndCaps
    if abs(l["dxy"]) > 0.05 + 0.05*EC: return False
    if abs(l["dz"])  > 0.1  + 0.1*EC:  return False
    return True

def removekey(d, key):
    r = dict(d)
    del r[key]
    return r

def photonVIDSelector( p, idVal, removedCuts=[] ):
    # https://twiki.cern.ch/twiki/bin/view/CMS/CutBasedPhotonIdentificationRun2

    vidDict = photonIdBitMapToDict( p["vidNestedWPBitmap"] )
    if not removedCuts:
        return all( [ cut >= idVal for cut in vidDict.values() ] )

#    # all cuts in removedCuts
#    if not list( set(vidNestedWPBitMapNamingListPhoton) - set(removedCuts) ):
#        return True

    # sieie is implemented according to the twiki
    if ("sieie"          in removedCuts):
        vidDict = removekey( vidDict, "PhoFull5x5SigmaIEtaIEtaCut" )
    # hoe is implemented according to the twiki
    if ("hoe"            in removedCuts):
        vidDict = removekey( vidDict, "PhoSingleTowerHadOverEmCut" )
    # PhoAnyPFIsoWithEACut seems to be chg-iso cut, pretty sure
    if ("pfRelIso03_chg" in removedCuts):
        vidDict = removekey( vidDict, "PhoAnyPFIsoWithEACut" )
    # PhoAnyPFIsoWithEACut2 seems to be very close to all-iso cut, but not everywhere, sometimes it's the neutral-iso cut
    # no clue what PhoAnyPFIsoWithEAAndQuadScalingCut is, I think its the all-iso cut with additional pt**2 dependence, hard to check
    if ("pfRelIso03_all" in removedCuts):
        vidDict = removekey( vidDict, "PhoAnyPFIsoWithEAAndQuadScalingCut" )
        vidDict = removekey( vidDict, "PhoAnyPFIsoWithEACut2" )
    # no clue what the scEtaMultiRange cut really is, probably out of scope for this analysis as we only use SC barrel photons
    if ("scEtaMultiRange" in removedCuts):
        vidDict = removekey( vidDict, "PhoSCEtaMultiRangeCut" )
    # no clue what the minPt cut really is
    if ("minPt"          in removedCuts):
        vidDict = removekey( vidDict, "MinPtCut" )
        
    return all( [ cut >= idVal for cut in vidDict.values() ] )

def electronVIDSelector( l, idVal, removedCuts=[] ):

    vidDict    = vidNestedWPBitMapToDict( l['vidNestedWPBitmap'] )
    if not removedCuts:
        return all( [ cut >= idVal for cut in vidDict.values() ] )

    if ("pt"             in removedCuts):
        vidDict = removekey( vidDict, "MinPtCut" )
    if ("sieie"          in removedCuts):
        vidDict = removekey( vidDict, "GsfEleFull5x5SigmaIEtaIEtaCut" )
    if ("hoe"            in removedCuts):
        vidDict = removekey( vidDict, "GsfEleHadronicOverEMEnergyScaledCut" )
    if ("pfRelIso03_all" in removedCuts):
        vidDict = removekey( vidDict, "GsfEleRelPFIsoScaledCut" )
    if ("SCEta" in removedCuts):
        vidDict = removekey( vidDict, "GsfEleSCEtaMultiRangeCut" )
    if ("dEtaSeed" in removedCuts):
        vidDict = removekey( vidDict, "GsfEleDEtaInSeedCut" )
    if ("dPhiInCut" in removedCuts):
        vidDict = removekey( vidDict, "GsfEleDPhiInCut" )
    if ("EinvMinusPinv" in removedCuts):
        vidDict = removekey( vidDict, "GsfEleEInverseMinusPInverseCut" )
    if ("convVeto" in removedCuts):
        vidDict = removekey( vidDict, "GsfEleConversionVetoCut" )
    if ("lostHits" in removedCuts):
        vidDict = removekey( vidDict, "GsfEleMissingHitsCut" )

    return all( [ cut >= idVal for cut in vidDict.values() ] )

# nanoAOD 94x or higher
def triggerEmulatorSelector( l, wp ):
    cutDict    = vidNestedWPBitMapToDict( l['vidNestedWPBitmap'] )
    wpCutValue = vidNestedWPBitMap[ wp ]     
    return all( val >= wpCutValue for val in cutDict.values() )

#def triggerEmulatorSelector( l, wp ):
#    # still missing: dPhiSC, dEtaSC (switch to triggerEmulatorCuts in nanoAOD 94x or higher)
#    ECSc = abs(l["eta"] + l["deltaEtaSC"]) > 1.479
#    if l["sieie"]            >= (0.011+0.019*ECSc): return False
#    if l["eInvMinusPInv"]    <= -0.05:              return False
#    if l["eInvMinusPInv"]    >= (0.01-0.005*ECSc):  return False
#    if l["hoe"]              >= (0.10-0.03*ECSc):   return False
#    return True

def barrelEndcapVeto( p ):
    if abs(p['pdgId']) == 11: absEta = abs(p["eta"] + p["deltaEtaSC"])   # eta supercluster
    else:                     absEta = abs(p["eta"])                     # eta
    return ( absEta > 1.566 or absEta < 1.4442 )

# Reco Selectors
def jetSelector( year ):
    # https://twiki.cern.ch/twiki/bin/viewauth/CMS/JetIDhttps://twiki.cern.ch/twiki/bin/viewauth/CMS/JetID
    if year == 2016:
        # According to AN-2017/197
        # jetID cuts, pT and eta cuts
        def func(j, removedCuts=[], ptVar="pt"):
#            if not j["cleanmask"]:                           return False # too much cleaning
            if not "pt" in removedCuts:
                if j[ptVar]      < 30:                      return False
            if not "eta" in removedCuts:
                if abs(j["eta"]) >= 2.4:                     return False
            if not jetIdBitMapToDict( j["jetId"] )["loose"]: return False
            return True
        return func

    elif year == 2017 or year == 2018:
        # jetID cuts, pT and eta cuts
        def func(j, removedCuts=[], ptVar="pt"):
#            if not j["cleanmask"]:                           return False
            if not "pt" in removedCuts:
                if j[ptVar]      < 30:                      return False
            if not "eta" in removedCuts:
                if abs(j["eta"]) >= 2.4:                     return False
            if not jetIdBitMapToDict( j["jetId"] )["tight"]: return False
            return True
        return func

    else:
        raise (NotImplementedError, "Don't know what cut to use for year %s"%year)

#muonRelIsoCut = 0.12
muonRelIsoCut      = 0.15
muonRelIsoCutLoose = 0.20
muonRelIsoCutVeto  = 0.25
def muonSelector( lepton_selection ):
    # According to AN-2017/197
    if lepton_selection == 'tight':
        def func(l, removedCuts=[], ptVar="pt"):
            if not "pt" in removedCuts:
                if l[ptVar]          <  30:                        return False
            if not "eta" in removedCuts:
                if abs(l["eta"])    > 2.4:                        return False
            if not l["tightId"]:                                  return False
#            if not "vertex" in removedCuts:
#                if not vertexSelector(l):                         return False #sync with FNAL
            if not "pfRelIso04_all" in removedCuts:
                if l['pfRelIso04_all']  > muonRelIsoCut:          return False
#            if not "sip3d" in removedCuts:
#                if l["sip3d"]           > 4:                      return False #sync with FNAL
#            if l['pfIsoId']        <  muonPfIsoId['PFIsoMedium']: return False
            return True
        return func

    elif lepton_selection == 'medium':
        def func(l, leading=False, removedCuts=[], ptVar="pt"):
            if not "pt" in removedCuts:
                if leading:
                    if l[ptVar]         <= 25:            return False
                else:
                    if l[ptVar]         <= 15:            return False
            if not "eta" in removedCuts:
                if abs(l["eta"])    > 2.4:               return False
            if not l["mediumId"]:                        return False
            if not "vertex" in removedCuts:
                if not vertexSelector(l):                return False
            if not "pfRelIso04_all" in removedCuts:
                if l['pfRelIso04_all']  > muonRelIsoCut: return False
            if not "sip3d" in removedCuts:
                if l["sip3d"]           > 4:             return False
#            if l['pfIsoId']        <  muonPfIsoId['PFIsoMedium']: return False
            return True
        return func

    elif lepton_selection == 'veto2l':
        # muon loose requirement
        def func(l, removedCuts=[], ptVar="pt"):
            if not "pt" in removedCuts:
                if l[ptVar]          < 15:                            return False
            if not "eta" in removedCuts:
                if abs(l["eta"])    > 2.4:                           return False
            if not "ID" in removedCuts:
                if not l["isPFcand"]:                                return False
                if not ( l["isGlobal"] or l["isTracker"] ):          return False
            if not "vertex" in removedCuts:
                if not vertexSelector(l):                            return False
            if not "pfRelIso04_all" in removedCuts:
                if l['pfRelIso04_all']  > 0.4:                       return False
            if not "sip3d" in removedCuts:
                if l["sip3d"]           > 4:                         return False
#            if l['pfIsoId']        <  muonPfIsoId['PFIsoVeryLoose']: return False
            return True
        return func

    elif lepton_selection == 'veto':
        # muon loose requirement
        def func(l, removedCuts=[], ptVar="pt"):
            if not "pt" in removedCuts:
                if l[ptVar]          < 15:                            return False
            if not "eta" in removedCuts:
                if abs(l["eta"])    > 2.4:                           return False
            if not "ID" in removedCuts:
                if not l["isPFcand"]:                                return False
                if not ( l["isGlobal"] or l["isTracker"] ):          return False
#            if not "vertex" in removedCuts:
#                if not vertexSelector(l):                            return False #sync with FNAL
            if not "pfRelIso04_all" in removedCuts:
                if l['pfRelIso04_all']  > muonRelIsoCutVeto:                      return False
#            if not "sip3d" in removedCuts:
#                if l["sip3d"]           > 4:                         return False #sync with FNAL
#            if l['pfIsoId']        <  muonPfIsoId['PFIsoVeryLoose']: return False
            return True
        return func

    else:
        raise (NotImplementedError, "Don't know what cut to use for muons")

def getElectronIsoCutV2( pt, eta, id ):
    if id == "tight":
        if abs(eta) <= 1.479: return 0.0287+0.506/pt
        else:                 return 0.0445+0.963/pt
    elif id == "medium":
        if abs(eta) <= 1.479: return 0.0478+0.506/pt
        else:                 return 0.0658+0.963/pt
    elif id == "loose":
        if abs(eta) <= 1.479: return 0.112+0.506/pt
        else:                 return 0.108+0.963/pt
    elif id == "veto":
        if abs(eta) <= 1.479: return 0.198+0.506/pt
        else:                 return 0.203+0.963/pt
    else:
        raise (NotImplementedError, "Don't know what cut to use for electrons")


# electrons 
def eleSelector( lepton_selection ):
    idVar = "cutBased"
    # According to AN-2017/197
    if lepton_selection == 'tight':
        def func(l, removedCuts=[], ptVar="pt"):
            if not barrelEndcapVeto(l):                      return False
            if not "pt" in removedCuts:
                if l[ptVar]          < 35:                    return False
            if not "eta" in removedCuts:
                if abs(l["eta"])    > 2.4:                   return False
            if not electronVIDSelector( l, vidNestedWPBitMap["tight"], removedCuts=removedCuts ): return False
#            if not "pfRelIso03_all" in removedCuts:
#                if l['pfRelIso03_all']  > 0.12:              return False
#            if not "sip3d" in removedCuts:
#                if l["sip3d"]           > 4:                 return False
            if not "vertex" in removedCuts:
                if not vertexSelector(l, barrelECdiff=True):                    return False
#            if l[idVar] < electronIdCutBased['tight']:       return False
#            if int(l["lostHits"])  != 0:                 return False # in the cutbased id? Ghent!
#            if not l["convVeto"]:                        return False # in the cutbased id? Ghent!
            return True
        return func

    elif lepton_selection == 'tight2l':
        def func(l, leading=False, removedCuts=[], ptVar="pt"):
            if not barrelEndcapVeto(l):                  return False
            if not "pt" in removedCuts:
                if leading:
                    if l[ptVar]         <= 25:            return False
                else:
                    if l[ptVar]         <= 15:            return False
            if not "eta" in removedCuts:
                if abs(l["eta"])   >= 2.4:               return False
            if not electronVIDSelector( l, vidNestedWPBitMap["tight"], removedCuts=removedCuts ): return False
            if not "pfRelIso03_all" in removedCuts:
                if l['pfRelIso03_all']  > 0.12:          return False
            if not "sip3d" in removedCuts:
                if l["sip3d"]           > 4:             return False
            if not "vertex" in removedCuts:
                if not vertexSelector(l):                return False
#            if l[idVar] < electronIdCutBased['tight']:   return False
#            if int(l["lostHits"])  != 0:                 return False # in the cutbased id? Ghent!
#            if not l["convVeto"]:                        return False # in the cutbased id? Ghent!
            return True
        return func

    elif lepton_selection == 'medium':
        def func(l, leading=False, removedCuts=[], ptVar="pt"):
            if not barrelEndcapVeto(l):                  return False
            if not "pt" in removedCuts:
                if leading:
                    if l[ptVar]         <= 25:            return False
                else:
                    if l[ptVar]         <= 15:            return False
            if not "eta" in removedCuts:
                if abs(l["eta"])   >= 2.4:               return False
#            if not electronVIDSelector( l, vidNestedWPBitMap["medium"], removedCuts=removedCuts ): return False
            if not "pfRelIso03_all" in removedCuts:
                if l['pfRelIso03_all']  > 0.12:          return False
            if not "sip3d" in removedCuts:
                if l["sip3d"]           > 4:             return False
            if not "vertex" in removedCuts:
                if not vertexSelector(l):                return False
            if l[idVar] < electronIdCutBased['medium']:  return False
#            if int(l["lostHits"])  != 0:                 return False # in the cutbased id? Ghent!
#            if not l["convVeto"]:                        return False # in the cutbased id? Ghent!
            return True
        return func

    elif lepton_selection == 'veto':
        def func(l, removedCuts=[], ptVar="pt"):
            if not barrelEndcapVeto(l):                return False
            if not "pt" in removedCuts:
                if l[ptVar]          < 15:              return False
            if not "eta" in removedCuts:
                if abs(l["eta"])    > 2.4:             return False
            if not electronVIDSelector( l, vidNestedWPBitMap["veto"], removedCuts=removedCuts ): return False
#            if not "pfRelIso03_all" in removedCuts:
#                if l['pfRelIso03_all']  > 0.4:         return False
#            if not "sip3d" in removedCuts:
#                if l["sip3d"]           > 4:           return False
            if not "vertex" in removedCuts:
                if not vertexSelector(l, barrelECdiff=True):              return False
#            if l[idVar] < electronIdCutBased['veto']:  return False
            return True
        return func

    else:
        raise (NotImplementedError, "Don't know what cut to use for electrons")

def tauSelector( lepton_selection ):
    # dummy function
        def func(l):
            return True
        return func

def photonSelector( selection, version ):
    # According to AN-2017/197
    if version == 2016:
        idVar    = "cutBased"       
        photonId = photonIdCutBased 
    else:
        idVar    =  "cutBasedBitmap"
        photonId =  photonIdCutBasedBitmap

    if selection == "mva":
        def func(g, removedCuts=[], ptVar="pt"):
            if not "pt" in removedCuts:
                if g[ptVar]       <= 20:        return False
#                if not g["isScEtaEB"]:        return False # Supercluster Barrel only
            if not "eta" in removedCuts:
                if abs(g["eta"]) >= 1.4442:    return False # Barrel only
            if not "pixelSeed" in removedCuts:
                if g["pixelSeed"]:             return False
#            if not "electronVeto" in removedCuts:
#                if not g["electronVeto"]:      return False
            if not "ID" in removedCuts:
                if not g["mvaID_WP90"]:        return False
            return True
        return func

    if selection == "medium":
        def func(g, removedCuts=[], ptVar="pt"):
            if not "pt" in removedCuts:
                if g[ptVar]       <= 20:                                                  return False
#                if not g["isScEtaEB"]:                                                   return False # Supercluster Barrel only
            if not "eta" in removedCuts:
                if abs(g["eta"]) >= 1.4442:                                              return False # Barrel only
            if not "pixelSeed" in removedCuts:
                if g["pixelSeed"]:                                                       return False
#            if not "electronVeto" in removedCuts:
#                if not g["electronVeto"]:                                                return False
#            if g[idVar]          <  photonId[selection]:                                 return False
            if not "ID" in removedCuts:
                if not photonVIDSelector( g, photonId[selection], removedCuts=removedCuts ): return False
            return True
        return func

    elif selection == "loose":
        def func(g, removedCuts=[], ptVar="pt"):
            if not "pt" in removedCuts:
                if g[ptVar]       <= 20:                                             return False
#                if not g["isScEtaEB"]:                                                   return False # Supercluster Barrel only
            if not "eta" in removedCuts:
                if abs(g["eta"]) >= 1.479:                                          return False # Barrel only
            if not "pixelSeed" in removedCuts:
                if g["pixelSeed"]:                                                  return False
#            if not "electronVeto" in removedCuts:
#                if not g["electronVeto"]:                                           return False
#            if g[idVar]          <  photonId[selection]:                            return False
            if not "ID" in removedCuts:
                if not photonVIDSelector( g, photonId[selection], removedCuts=removedCuts ): return False
            return True
        return func

    else:
        raise (NotImplementedError, "Don't know what cut to use for photons")


# Gen Selectors
def genJetSelector( selection=None ):

    if selection == 'ATLASUnfolding':
        def func(j):
            if j["pt"]       < 25:  return False
            if abs(j["eta"]) > 2.5: return False
            return True
        return func

    elif selection == 'CMSUnfolding':
        def func(j):
            if j["pt"]       < 30:  return False
            if abs(j["eta"]) > 2.4: return False
            return True
        return func

    else:
        # According to AN-2017/197
        def func(j):
            if j["pt"]       <= 25:  return False
            if abs(j["eta"]) >= 2.5: return False
            return True
    return func

def genLeptonSelector( selection=None ):

    if selection == 'ATLASUnfolding':
        def func(l):
            if abs(l["pdgId"]) not in [11,13]: return False
            if l["pt"]       < 25:  return False
            if abs(l["eta"]) > 2.5: return False
            return True
        return func

    elif selection == 'CMSUnfolding':
        def func(l):
            if abs(l["pdgId"]) not in [11,13]: return False
            if abs(l["pdgId"]) == 11:
                if l["pt"]       < 35:         return False
                if abs(l["eta"]) > 2.4:        return False
            elif abs(l["pdgId"]) == 13:
                if l["pt"]       < 30:         return False
                if abs(l["eta"]) > 2.4:        return False
            return True
        return func

    else:
        # According to AN-2017/197
        def func(l):
            if l["pt"]       <= 25:  return False
            if abs(l["eta"]) >= 2.5: return False
            return True
    return func

def genPhotonSelector( photon_selection=None ):
    # According to AN-2017/197
    if photon_selection == 'overlapTTGamma':
        # Remove events from ttbar sample, keep ttgamma events
        def func(g):
            if g["pt"]       < 10:  return False
            if abs(g["eta"]) > 5.0: return False
            return True
        return func

    elif photon_selection == 'overlapZWGamma':
        # Remove events from DY and W+jets sample, keep Zgamma and Wgamma events
        def func(g):
            if g["pt"]       < 15:  return False
            if abs(g["eta"]) > 2.6: return False
            return True
        return func

    elif photon_selection == 'overlapZGamma':
        # Remove events from DY, keep Zgamma lowGPt
        def func(g):
            if g["pt"]       < 9:  return False
            if abs(g["eta"]) > 2.6: return False
            return True
        return func

    elif photon_selection == 'overlapSingleTopTch':
        # Remove events from single top t-channel sample, keep single top + photon events
        def func(g):
            if g["pt"]       < 10:  return False
            if abs(g["eta"]) > 2.6: return False
            return True
        return func

    elif photon_selection == 'overlapGJets':
        # Remove events from QCD, keep gamma+jets events
        def func(g):
            if g["pt"]       < 25:  return False
            if abs(g["eta"]) > 2.5: return False
            return True
        return func

    elif photon_selection == 'ATLASUnfolding':
        def func(g):
            if g["pt"]       < 20:   return False
            if abs(g["eta"]) > 2.37: return False
            return True
        return func

    elif photon_selection == 'CMSUnfolding':
        def func(g):
            if g["pt"]       < 20:     return False
            if abs(g["eta"]) > 1.4442: return False
            return True
        return func

    else:
        # general gen-photon selection
        def func(g):
            if g["pt"]       < 13:  return False
            if abs(g["eta"]) > 5.0: return False
            return True
        return func

# Gen Particle Filter
def filterGenElectrons( genParts, status=None ):
    if   status == 'first': stat = [23]
    elif status == 'last':  stat = [1]
    else:                   stat = [1,23]
    electrons = list( filter( lambda l: abs(l['pdgId']) == 11 and l['status'] in stat, genParts ) )
    return electrons

def filterGenMuons( genParts, status=None ):
    if   status == 'first': stat = [23]
    elif status == 'last':  stat = [1]
    else:                   stat = [1,23]
    muons = list( filter( lambda l: abs(l['pdgId']) == 13 and l['status'] in stat, genParts ) )
    return muons

def filterGenTaus( genParts, status=None ):
    if   status == 'first': stat = [23]
    elif status == 'last':  stat = [1]
    else:                   stat = [1,23]
    taus = list( filter( lambda l: abs(l['pdgId']) == 15 and l['status'] in stat, genParts ) )
    return taus

def filterGenPhotons( genParts, status=None ):
    if   status == 'first': stat = [23]
    elif status == 'last':  stat = [1]
    if status == "all":
        photons = list( filter( lambda l: abs(l['pdgId']) == 22 and l['status'] >= 0, genParts ) )
    else:
        photons = list( filter( lambda l: abs(l['pdgId']) == 22 and l['status'] in stat, genParts ) )
    return photons

def filterGenTops( genParts ):
    tops = list( filter( lambda l: abs(l['pdgId']) == 6 and l['status'] == 62, genParts ) )
    return tops

def filterGenBJets( genJets ):
    bjets = list( filter( lambda j: abs(j['partonFlavour']) == 5, genJets ) )
    return bjets

# Pythia status flags:
#  1: outgoing particle (e.g. photon from decay products of top)
# 23: e.g. intermediate photon (mother of photon with status 1)
# 44: ISR photons (outgoing, but initially produced and then nothing changed (MG photons?)
# 52: photons from e.g. top, probably also initially produced (MG photons?)
# 51: photons from e.g. top, probably also initially produced (MG photons?)

# 62: outgoing subprocess particle with primordial kT included, last gencopy before it decays
# 23: outgoing
# 71: copied partons to collect into contiguous colour singlet (jets)
#  1: stage of event generation inside PYTHIA
# 22: intermediate (intended to have preserved mass) (tops)


# for 4 muon default plot script

muonVars_data = ['pt','eta','phi','pdgId','mediumId','miniPFRelIso_all','pfRelIso03_all','sip3d','dxy','dz','charge']
muonVars = muonVars_data + []

electronVars_data = ['pt','eta','phi','pdgId','cutBased','miniPFRelIso_all','pfRelIso03_all','sip3d','lostHits','convVeto','dxy','dz','charge','deltaEtaSC','mvaFall17V2noIso_WP80', 'vidNestedWPBitmap']
electronVars = electronVars_data + []

def alwaysTrue(*args, **kwargs):
  return True

def alwaysFalse(*args, **kwargs):
  return False

def getMuons(c, collVars=muonVars):
    return [getObjDict(c, 'Muon_', collVars, i) for i in range(int(getVarValue(c, 'nMuon')))]

def getElectrons(c, collVars=electronVars):
    return [getObjDict(c, 'Electron_', collVars, i) for i in range(int(getVarValue(c, 'nElectron')))]

def getGoodMuons(c, collVars=muonVars, mu_selector = alwaysFalse):
    return [l for l in getMuons(c, collVars) if mu_selector(l)]

def getGoodElectrons(c, collVars=electronVars, ele_selector = alwaysFalse):
    return [l for l in getElectrons(c, collVars) if ele_selector(l)]

