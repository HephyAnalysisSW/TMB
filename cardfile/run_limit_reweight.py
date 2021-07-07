import ROOT
import os
import argparse
from RootTools.core.Sample import Sample
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--logLevel',       action='store',      default='INFO',         nargs='?', choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'],             help="Log level for logging")
argParser.add_argument("--signal",         action='store',      default='dipoles',       nargs='?', choices=["dipoles", "currents"], help="which signal scan?")
argParser.add_argument("--model",          action='store',      default='dim6top_LO',   nargs='?', choices=["dim6top_LO", "ewkDM"], help="which signal model?")
argParser.add_argument("--only",           action='store',      default=None,           nargs='?',                                                                                           help="pick only one signal point?")
argParser.add_argument("--scale",          action='store',      default=1.0,            type=float,nargs='?',                                                                                help="scaling all yields")
argParser.add_argument("--overwrite",      action='store_true', default = False,        help="Overwrite existing output files")
argParser.add_argument("--useXSec",        action='store_true', help="Use the x-sec information?")
argParser.add_argument("--useShape",       action='store_true', help="Use the shape information?")
argParser.add_argument("--statOnly",       action='store_true', help="Use only statistical uncertainty?")
argParser.add_argument("--controlRegion",  action='store',      default='', choices = ['', 'nbtag0-njet3p', 'nbtag1p-njet02', 'nbtag1p-njet2', 'nbtag0-njet02', 'nbtag0-njet0p', 'nbtag0-njet1p', 'nbtag0-njet2p'], help="Use any CRs cut?")
argParser.add_argument("--includeCR",      action='store_true', help="Do simultaneous SR and CR fit")
argParser.add_argument("--inclusiveRegions", action='store_true', help="Use inclusive signal regions?")
argParser.add_argument("--calcNuisances",  action='store_true', help="Extract the nuisances and store them in text files?")
argParser.add_argument("--unblind",        action='store_true', help="Unblind? Currently also correlated with controlRegion option for safety.")
argParser.add_argument("--WZtoPowheg",     action='store_true', help="Use reweighting from WZ amc@NLO sample to powheg?")
argParser.add_argument("--expected",       action='store_true', help="Run expected NLL (=blinded, no pseudo-data needed)?")
argParser.add_argument("--merged",         action='store_true', help="Run expected NLL (=blinded, no pseudo-data needed)?")
argParser.add_argument("--year",           action='store',      type=int, default=2016, choices = [ 2016, 2017, 20167 ], help='Which year?')


args = argParser.parse_args()


# Logging
import Analysis.Tools.logger as logger
logger = logger.get_logger(args.logLevel, logFile = None )
import RootTools.core.logger as logger_rt
logger_rt = logger_rt.get_logger(args.logLevel, logFile = None )

# MVA_HISTOS #
filename = '/mnt/hephy/cms/rosmarie.schoefbeck/www/tWZ/plots/analysisPlots/ttZ_dy_noData/Run2016/all/dilepM-onZ1-minDLmass12-njet5p-btag1p/ttz_dy_TTZ.root'
f = ROOT.TFile.Open(filename)
canvas = f.Get(f.GetListOfKeys().At(0).GetName())
#number of bins 
nbins = 50 
#signal
sig = canvas.GetListOfPrimitives().At(1)
#backround(s)
bkg = [] 
bkg.append( canvas.GetListOfPrimitives().At(2) ) 

from math                               import sqrt
from copy                               import deepcopy

from TMB.Tools.user                     import combineReleaseLocation, results_directory, plot_directory
from Ananlysis.Tools.CardFileWriter     import CardFileWriter

year = int(args.year)
limitDir = '/mnt/hephy/cms/rosmarie.schoefbeck/cardfiles'
    
uncertainties = {}

def wrapper(s):
    
    logger.info("Now working on %s", s.name)
    xSecScale = 1
    c = cardFileWriter.cardFileWriter()
    c.releaseLocation = combineReleaseLocation
    
    cardFileName = os.path.join(limitDir, s.name+'.txt')
    if not os.path.exists(cardFileName) or overWrite:
        counter=0
        c.reset()
        c.setPrecision(3)
        postfix = '_%s'%args.year
        c.addUncertainty('PU',                  'lnN') # correlated
        c.addUncertainty('JEC',                 'lnN') # correlated
        c.addUncertainty('JER',                 'lnN') # correlated
        c.addUncertainty('btag_heavy'+postfix,  'lnN') # uncorrelated, wait for offical recommendation
        c.addUncertainty('btag_light'+postfix,  'lnN') # uncorrelated, wait for offical recommendation
        c.addUncertainty('trigger'+postfix,     'lnN') # uncorrelated, statistics dominated
        c.addUncertainty('leptonSFSyst',        'lnN') # correlated
        c.addUncertainty('leptonTracking',      'lnN') # correlated
        c.addUncertainty('eleSFStat'+postfix,   'lnN') # uncorrelated
        c.addUncertainty('muSFStat'+postfix,    'lnN') # uncorrelated
        c.addUncertainty('scale',               'lnN') # correlated.
        c.addUncertainty('scale_sig',           'lnN') # correlated.
        c.addUncertainty('PDF',                 'lnN') # correlated.
        c.addUncertainty('PartonShower',        'lnN') # correlated.
        c.addUncertainty('nonprompt',           'lnN') # correlated?!
        c.addUncertainty('WZ_xsec',             'lnN') # correlated.
        c.addUncertainty('WZ_bb',               'lnN') # correlated
        c.addUncertainty('WZ_powheg',           'lnN') # correlated
        c.addUncertainty('WZ_njet',             'lnN') # correlated
        c.addUncertainty('ZZ_xsec',             'lnN') # correlated.
        c.addUncertainty('XG_xsec',             'lnN') # correlated.
        c.addUncertainty('rare',                'lnN') # correlated.
        c.addUncertainty('ttX',                 'lnN') # correlated.
        c.addUncertainty('Lumi'+postfix,        'lnN')

        uncList = ['PU', 'JEC', 'btag_heavy', 'btag_light', 'leptonSFSyst', 'trigger']
        for unc in uncList:
            uncertainties[unc] = []
        
        signal      = sig

        for b in nbins:
            #signalevents= sig.GetBinContent(b)
            totalBackground = u_float(0)

            niceName = sig.name  

            binname = 'Bin'+str(counter)
            logger.info("Working on %s", binname)
            counter += 1
            c.addBin(binname, [e.name for e in bkg], niceName)
            #c.addBin(binname, 'nonPromptDD', niceName)

            for e in bkg:  
                name = e.name
                expected = e.GetBinContent(b) 
                logger.info("Adding expectation %s for process %s", expected.val, name)
                c.specifyExpectation(binname, name, expected.val if expected.val > 0.01 else 0.01)

                totalBackground += expected

                if not args.statOnly:
                    # uncertainties
                    pu          = 1 + e.PUSystematic( r, channel, setup).val            if expected.val>0.01 else 1.1
                    jec         = 1 + e.JECSystematic( r, channel, setup).val           if expected.val>0.01 else 1.1
                    jer         = 1 + e.JERSystematic( r, channel, setup).val           if expected.val>0.01 else 1.1
                    btag_heavy  = 1 + e.btaggingSFbSystematic(r, channel, setup).val    if expected.val>0.01 else 1.1
                    btag_light  = 1 + e.btaggingSFlSystematic(r, channel, setup).val    if expected.val>0.01 else 1.1
                    trigger     = 1 + e.triggerSystematic(r, channel, setup).val        if expected.val>0.01 else 1.1
                    leptonSFSyst= 1 + e.leptonSFSystematic(r, channel, setup).val       if expected.val>0.01 else 1.1
                    leptonReco  = 1 + e.leptonTrackingSystematic(r, channel, setup).val if expected.val>0.01 else 1.1
                    eleSFStat   = 1 + e.eleSFSystematic(r, channel, setup).val          if expected.val>0.01 else 1.1
                    muSFStat    = 1 + e.muSFSystematic(r, channel, setup).val           if expected.val>0.01 else 1.1

                    c.specifyUncertainty('PU',          binname, name, 1 + e.PUSystematic( r, channel, setup).val)

                    if not name.count('nonprompt'):
                        c.specifyUncertainty('JEC',                 binname, name, jec)
                        c.specifyUncertainty('JER',                 binname, name, jer)
                        c.specifyUncertainty('btag_heavy'+postfix,  binname, name, btag_heavy)
                        c.specifyUncertainty('btag_light'+postfix,  binname, name, btag_light)
                        c.specifyUncertainty('trigger'+postfix,     binname, name, trigger)
                        c.specifyUncertainty('leptonSFSyst',        binname, name, leptonSFSyst)
                        c.specifyUncertainty('leptonTracking',      binname, name, leptonReco)
                        c.specifyUncertainty('eleSFStat'+postfix,   binname, name, eleSFStat)
                        c.specifyUncertainty('muSFStat'+postfix,    binname, name, muSFStat)
                        c.specifyUncertainty('scale',               binname, name, 1.01) 
                        c.specifyUncertainty('PDF',                 binname, name, 1.01)
                        if year == 2016:
                            c.specifyUncertainty('Lumi'+postfix,        binname, name, 1.025 )
                        else:
                            c.specifyUncertainty('Lumi'+postfix,        binname, name, 1.023 )

                    if name.count('ZZ'):    c.specifyUncertainty('ZZ_xsec',     binname, name, 1.10)
                    if name.count('XG'):    c.specifyUncertainty('XG_xsec',     binname, name, 1.20)
                    if name.count('WZ') and not name.count("WZZ") and not name.count("WWZ"):
                        c.specifyUncertainty('WZ_xsec',     binname, name, 1.10)
                        njetUnc = e.highNJetSystematic(r, channel, setup).val
                        if njetUnc>0:
                            c.specifyUncertainty('WZ_njet',     binname, name, 1+njetUnc)
                        if setup == setup3l:
                            c.specifyUncertainty('WZ_bb',     binname, name, 1.16)
                        c.specifyUncertainty('WZ_powheg',     binname, name, WZ_powheg)
                    
                    if name.count('nonprompt'):    c.specifyUncertainty('nonprompt',   binname, name, 1.30)
                    if name.count('rare'):    c.specifyUncertainty('rare',        binname, name, 1.50)
                    if name.count('TTX'):     c.specifyUncertainty('ttX',         binname, name, 1.11)


                #MC bkg stat (some condition to neglect the smaller ones?)
                uname = 'Stat_'+binname+'_'+name+postfix
                c.addUncertainty(uname, 'lnN')
                if expected.val > 0:
                    c.specifyUncertainty(uname, binname, name, 1 + expected.sigma/expected.val )
                else:
                    c.specifyUncertainty(uname, binname, name, 1.01 )
            
            uname = 'Stat_'+binname+'_nonprompt'+postfix
            c.addUncertainty(uname, 'lnN')
            
#            if setup.nLeptons == 3 and setupNP:
#                nonprompt   = FakeEstimate(name="nonPromptDD_%s"%args.year, sample=setup.samples["Data"], setup=setupNP, cacheDir=setup.defaultCacheDir())
#                np = nonprompt.cachedEstimate(r, channel, setupNP)
#                if np.val < 0.01:
#                    np = u_float(0.01,0.)
#                c.specifyExpectation(binname, 'nonPromptDD', np.val ) 
#                c.specifyUncertainty(uname,   binname, "nonPromptDD", 1 + np.sigma/np.val )
#                c.specifyUncertainty('nonprompt',   binname, "nonPromptDD", 1.30)
#            else:
#                np = u_float(0)
#                c.specifyExpectation(binname, 'nonPromptDD', np.val)
            
#            if args.expected:
#                sig = signal.cachedEstimate(r, channel, setup)
#                obs = totalBackground + sig + np
#            elif args.unblind or (setup == setup3l_CR) or (setup == setup4l_CR):
#                obs = observation.cachedObservation(r, channel, setup)
#            else:
#                obs = observation.cachedEstimate(r, channel, setup)
            obs = observation.cachedEstimate(r, channel, setup)
            c.specifyObservation(binname, int(round(obs.val,0)) )

#            if args.useShape:
#                logger.info("Using 2D reweighting method for shapes")
#                if args.model == "dim6top_LO":
#                    source_gen = dim6top_central
#                elif args.model == "ewkDM":
#                    source_gen = ewkDM_central
#
#                signalReweighting = SignalReweighting( source_sample = source_gen, target_sample = s, cacheDir = reweightCache)
#                f = signalReweighting.cachedReweightingFunc( setup.genSelection )
#                sig = signal.reweight2D(r, channel, setup, f)
#            else:
#                sig = signal.cachedEstimate(r, channel, setup)
            sig = signal.cachedEstimate(r, channel, setup)

            xSecMod = 1
            if args.useXSec:
                xSecMod = xsec.val/xsec_central.val
            
            logger.info("x-sec is multiplied by %s",xSecMod)
            
            signalName = 'signal'
            #c.specifyExpectation(binname, 'ttZ', 0) # this is just a fake signal for combine
            c.specifyExpectation(binname, signalName, sig.val * xSecScale * xSecMod ) # this is the real signal 
            logger.info('Adding signal %s'%(sig.val * xSecScale * xSecMod))
            
            if sig.val>0:
                if year == 2016:
                    c.specifyUncertainty('Lumi'+postfix, binname, signalName, 1.025 )
                else:
                    c.specifyUncertainty('Lumi'+postfix, binname, signalName, 1.023 )
                if not args.statOnly:
                    # uncertainties
                    pu          = 1 + signal.PUSystematic( r, channel, setup).val
                    jec         = 1 + signal.JECSystematic( r, channel, setup).val
                    jer         = 1 + signal.JERSystematic( r, channel, setup).val
                    btag_heavy  = 1 + signal.btaggingSFbSystematic(r, channel, setup).val
                    btag_light  = 1 + signal.btaggingSFlSystematic(r, channel, setup).val
                    trigger     = 1 + signal.triggerSystematic(r, channel, setup).val
                    leptonSFSyst= 1 + signal.leptonSFSystematic(r, channel, setup).val
                    leptonReco  = 1 + signal.leptonTrackingSystematic(r, channel, setup).val
                    eleSFStat   = 1 + signal.eleSFSystematic(r, channel, setup).val
                    muSFStat    = 1 + signal.muSFSystematic(r, channel, setup).val 

                    if sig.sigma/sig.val < 0.05:
                        uncertainties['PU']         += [pu]
                        uncertainties['JEC']        += [jec]
                        uncertainties['btag_heavy'] += [btag_heavy]
                        uncertainties['btag_light'] += [btag_light]
                        uncertainties['trigger']    += [trigger]
                        uncertainties['leptonSFSyst']   += [leptonSFSyst]

                    c.specifyUncertainty('PU',                  binname, signalName, pu)
                    c.specifyUncertainty('JEC',                 binname, signalName, jec)
                    c.specifyUncertainty('JER',                 binname, signalName, jer)
                    c.specifyUncertainty('btag_heavy'+postfix,  binname, signalName, btag_heavy)
                    c.specifyUncertainty('btag_light'+postfix,  binname, signalName, btag_light)
                    c.specifyUncertainty('trigger'+postfix,     binname, signalName, trigger)
                    c.specifyUncertainty('leptonSFSyst',        binname, signalName, leptonSFSyst)
                    c.specifyUncertainty('leptonTracking',      binname, signalName, leptonReco)
                    c.specifyUncertainty('eleSFStat'+postfix,   binname, signalName, eleSFStat)
                    c.specifyUncertainty('muSFStat'+postfix,    binname, signalName, muSFStat)
                    # This doesn't get the right uncertainty in CRs. However, signal doesn't matter there anyway.
                    if setup in [setup3l, setup4l]:
                        if args.inclusiveRegions:
                            c.specifyUncertainty('scale_sig',   binname, signalName, 1.02)
                            c.specifyUncertainty('PDF',         binname, signalName, 1.04)
                            c.specifyUncertainty('PartonShower',binname, signalName, 1.025)
                        else:
                            c.specifyUncertainty('scale_sig',   binname, signalName, 1 + scale_cache.get({"region":r, "channel":channel.name, "PDFset":"scale"}).val)
                            c.specifyUncertainty('PDF',         binname, signalName, 1 + PDF_cache.get({"region":r, "channel":channel.name, "PDFset":PDFset}).val)
                            c.specifyUncertainty('PartonShower',binname, signalName, PS_cache.get({"region":r, "channel":channel.name, "PDFset":"PSscale"}).val) #something wrong here?
                    #c.specifyUncertainty('scale_sig',   binname, "signal", 1.05) #1.30
                    #c.specifyUncertainty('PDF',         binname, "signal", 1.04) #1.15

                uname = 'Stat_'+binname+'_%s'%signalName+postfix
                c.addUncertainty(uname, 'lnN')
                c.specifyUncertainty(uname, binname, signalName, 1 + sig.sigma/sig.val )
            else:
                uname = 'Stat_'+binname+'_%s'%signalName+postfix
                c.addUncertainty(uname, 'lnN')
                c.specifyUncertainty(uname, binname, signalName, 1 )

            
        #c.addUncertainty('Lumi'+postfix, 'lnN')
        #c.specifyFlatUncertainty('Lumi'+postfix, 1.026)
        cardFileName = c.writeToFile(cardFileName)
    else:
        logger.info("File %s found. Reusing.",cardFileName)
    
