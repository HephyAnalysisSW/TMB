import ROOT
import os
import argparse
from RootTools.core.Sample import Sample
argParser = argparse.ArgumentParser(description = "Argument parser")
argParser.add_argument('--logLevel',       action='store',      default='INFO',         nargs='?', choices=['CRITICAL', 'ERROR', 'WARNING', 'INFO', 'DEBUG', 'TRACE', 'NOTSET'],             help="Log level for logging")
argParser.add_argument("--year",           action='store',      type=int, default=2016, choices = [ 2016, 2017, 20167 ], help='Which year?')
argParser.add_argument("--lumi",           action='store',      type=float, default=35.9 , help='Which lumi?')
argParser.add_argument("--overwrite",      action='store_true', default = False,        help="Overwrite existing output files")


args = argParser.parse_args()


# Logging
import Analysis.Tools.logger as logger
logger = logger.get_logger(args.logLevel, logFile = None )
import RootTools.core.logger as logger_rt
logger_rt = logger_rt.get_logger(args.logLevel, logFile = None )

# MVA_HISTOS #
plotfilename = []
#plotfilename.append( '/mnt/hephy/cms/rosmarie.schoefbeck/www/tWZ/plots/analysisPlots/TTZ_DY_noData/Run2016/all/dilepM-onZ1-minDLmass12-njet5p-btag1p/ttz_dy_TTZ.root' )
#plotfilename.append( '/mnt/hephy/cms/rosmarie.schoefbeck/www/tWZ/plots/analysisPlots/TTZ_DY_noData/Run2016/all/dilepM-onZ1-minDLmass12-njet5p-btag1p/ttz_dy_LSTM_TTZ.root' )
plotfilename.append( '/mnt/hephy/cms/rosmarie.schoefbeck/www/tWZ/plots/analysisPlots/TTZ_2l_noData/Run2016/all/dilepM-onZ1-minDLmass12-njet5p-btag1p/ttz_2l_TTZscoarse.root' )
plotfilename.append( '/mnt/hephy/cms/rosmarie.schoefbeck/www/tWZ/plots/analysisPlots/TTZ_2l_noData/Run2016/all/dilepM-onZ1-minDLmass12-njet5p-btag1p/ttz_2l_LSTM_TTZscoarse.root' )

for filename in plotfilename: 
    f = ROOT.TFile.Open(filename)
    canvas = f.Get(f.GetListOfKeys().At(0).GetName())
    #number of bins 
    nbins = 10 
    #signal
    sig = canvas.GetListOfPrimitives().At(4)
    #backround(s)
    #bkg = [] 
    #bkg.append( canvas.GetListOfPrimitives().At(2) ) 
    #all 
    estimates = {}
    estimates['signal']= sig 
    estimates['DY']=canvas.GetListOfPrimitives().At(2)
    estimates['tt']=canvas.GetListOfPrimitives().At(3)
    print estimates 
    #estimates.append(sig)
    #estimates.append(canvas.GetListOfPrimitives().At(2))
    
    overWrite   = args.overwrite
    
    from math                               import sqrt
    from copy                               import deepcopy
    
    from TMB.Tools.user                     import combineReleaseLocation, results_directory, plot_directory
    from Analysis.Tools                     import u_float
    from Analysis.Tools.cardFileWriter      import cardFileWriter 
    #from Analysis.Tools      import CardFileWriter 
    
    year = int(args.year)
    limitDir = '/mnt/hephy/cms/rosmarie.schoefbeck/cardfiles'
        
    uncertainties = {}
    
    def wrapper(s):
        
        logger.info("Now working on %s", s)
        xSecScale = 1
        c = cardFileWriter.cardFileWriter()
    #    c = CardFileWriter.CardFileWriter()
        c.releaseLocation = combineReleaseLocation
        
        cardFileName = os.path.join(limitDir, s+'.txt')
        if not os.path.exists(cardFileName) or overWrite:
            counter=0
            c.reset()
            c.setPrecision(3)
            postfix = '_%s'%args.year
    #        c.addUncertainty('PU',                  'lnN') # correlated
            c.addUncertainty('JEC',                 'lnN') # correlated
            c.addUncertainty('JER',                 'lnN') # correlated
            c.addUncertainty('btag_heavy'+postfix,  'lnN') # uncorrelated, wait for offical recommendation
            c.addUncertainty('btag_light'+postfix,  'lnN') # uncorrelated, wait for offical recommendation
            c.addUncertainty('trigger'+postfix,     'lnN') # uncorrelated, statistics dominated
    #        c.addUncertainty('leptonSFSyst',        'lnN') # correlated
    #        c.addUncertainty('leptonTracking',      'lnN') # correlated
    #        c.addUncertainty('eleSFStat'+postfix,   'lnN') # uncorrelated
    #        c.addUncertainty('muSFStat'+postfix,    'lnN') # uncorrelated
            c.addUncertainty('scale',               'lnN') # correlated.
    #        c.addUncertainty('scale_sig',           'lnN') # correlated.
            c.addUncertainty('PDF',                 'lnN') # correlated.
    #        c.addUncertainty('PartonShower',        'lnN') # correlated.
    #        c.addUncertainty('nonprompt',           'lnN') # correlated?!
    #        c.addUncertainty('WZ_xsec',             'lnN') # correlated.
    #        c.addUncertainty('WZ_bb',               'lnN') # correlated
    #        c.addUncertainty('WZ_powheg',           'lnN') # correlated
    #        c.addUncertainty('WZ_njet',             'lnN') # correlated
    #        c.addUncertainty('ZZ_xsec',             'lnN') # correlated.
    #        c.addUncertainty('XG_xsec',             'lnN') # correlated.
    #        c.addUncertainty('rare',                'lnN') # correlated.
    #        c.addUncertainty('ttX',                 'lnN') # correlated.
            c.addUncertainty('Lumi'+postfix,        'lnN')
    
            uncList = ['PU', 'JEC', 'btag_heavy', 'btag_light', 'leptonSFSyst', 'trigger']
            for unc in uncList:
                uncertainties[unc] = []
            
            signal      = sig
            siglist     = []
    
            for b in range(nbins):
                #signalevents= sig.GetBinContent(b)
                totalBackground =  0. # u_float(0)
                sigandback = 0.
                niceName = s  
                binname = 'Bin'+str(counter)
                logger.info("Working on %s", binname)
                print binname 
                for e in estimates: print e 
                counter += 1
                c.addBin(binname, [e for e in estimates if e != "signal"], niceName)
    
                for e in estimates:  
                    expected = estimates[e].GetBinContent(b+1) 
                    if e == 'signal':
                        expected = expected - estimates['DY'].GetBinContent(b+1)

                    if args.lumi: 
                        # lumi_year = {2016: 35900.0, 2017: 41500.0, 2018: 59970.0}
                        expected = (expected/35.900)*args.lumi 

                    name = e

#                    if e == 'signal':
#                        expected = expected - (estimates['DY'].GetBinContent(b+1)/35.900)*args.lumi


                    print expected, b,estimates[e].GetBinContent(b+1)
                    logger.info("Adding expectation %s for process %s", expected, name)

                    c.specifyExpectation(binname, name, expected if expected > 0.01 else 0.01)
    
                    sigandback += expected
    
    #                # uncertainties
    #                pu          = 1 + e.PUSystematic( r, channel, setup).val            if expected.val>0.01 else 1.1
    #                jec         = 1 + e.JECSystematic( r, channel, setup).val           if expected.val>0.01 else 1.1
    #                jer         = 1 + e.JERSystematic( r, channel, setup).val           if expected.val>0.01 else 1.1
    #                btag_heavy  = 1 + e.btaggingSFbSystematic(r, channel, setup).val    if expected.val>0.01 else 1.1
    #                btag_light  = 1 + e.btaggingSFlSystematic(r, channel, setup).val    if expected.val>0.01 else 1.1
    #                trigger     = 1 + e.triggerSystematic(r, channel, setup).val        if expected.val>0.01 else 1.1
    #                leptonSFSyst= 1 + e.leptonSFSystematic(r, channel, setup).val       if expected.val>0.01 else 1.1
    #                leptonReco  = 1 + e.leptonTrackingSystematic(r, channel, setup).val if expected.val>0.01 else 1.1
    #                eleSFStat   = 1 + e.eleSFSystematic(r, channel, setup).val          if expected.val>0.01 else 1.1
    #                muSFStat    = 1 + e.muSFSystematic(r, channel, setup).val           if expected.val>0.01 else 1.1
    #
    #                c.specifyUncertainty('PU',          binname, name, 1 + e.PUSystematic( r, channel, setup).val)
    
                    c.specifyUncertainty('JEC',                 binname, name, 1.09)
                    c.specifyUncertainty('JER',                 binname, name, 1.01)
                    c.specifyUncertainty('btag_heavy'+postfix,  binname, name, 1.04)
                    c.specifyUncertainty('btag_light'+postfix,  binname, name, 1.04)
                    c.specifyUncertainty('trigger'+postfix,     binname, name, 1.01)
    #                c.specifyUncertainty('leptonSFSyst',        binname, name, leptonSFSyst)
    #                c.specifyUncertainty('leptonTracking',      binname, name, leptonReco)
    #                c.specifyUncertainty('eleSFStat'+postfix,   binname, name, eleSFStat)
    #                c.specifyUncertainty('muSFStat'+postfix,    binname, name, muSFStat)
                    c.specifyUncertainty('scale',               binname, name, 1.01) 
                    c.specifyUncertainty('PDF',                 binname, name, 1.01)
                    if year == 2016:
                        c.specifyUncertainty('Lumi'+postfix,        binname, name, 1.025 )
                    else:
                        c.specifyUncertainty('Lumi'+postfix,        binname, name, 1.023 )
    
                    #eg.
    #                if name.count('TTZ'):    c.specifyUncertainty('TTZ_xsec',     binname, name, 1.10)
    
                    #MC bkg stat (some condition to neglect the smaller ones?)
    #                uname = 'Stat_'+binname+'_'+name+postfix
    #
    #                c.addUncertainty(uname, 'lnN')
    #                if expected > 0:
    #                    c.specifyUncertainty(uname, binname, name, 1 + expected.sigma/expected.val )
    #                else:
    #                    c.specifyUncertainty(uname, binname, name, 1.01 )
                
    #            if setup.nLeptons == 3 and setupNP:
    #                nonprompt   = FakeEstimate(name="nonPromptDD_%s"%args.year, sample=setup.samples["Data"], setup=setupNP, cacheDir=setup.defaultCacheDir())
    #                np = nonprompt.cachedEstimate(r, channel, setupNP)
    #                if np.val < 0.01:
    #                    np = u_float(0.01,0.)
    #                c.specifyExpectation(binname, 'nonPromptDD', np.val ) 
    #                c.specifyUncertainty(uname,   binname, "nonPromptDD", 1 + nckgroundp.sigma/np.val )
    #                c.specifyUncertainty('nonprompt',   binname, "nonPromptDD", 1.30)
    #            else:
    #                np = u_float(0)
    #                c.specifyExpectation(binname, 'nonPromptDD', np.val)
                
                obs = sigandback # + sig.GetBinContent(b) 
                c.specifyObservation(binname, int(round(obs,0)) )
    
                signalName = 'signal'
                #sigval = sig.GetBinContent(b)
                #c.specifyExpectation(binname, 'ttZ', 0) # this is just a fake signal for combine
                #c.specifyExpectation(binname, signalName, sigval * xSecScale * xSecMod ) # this is the real signal 
    #            c.specifyExpectation(binname, signalName, sigval * xSecScale ) # this is the real signal 
                #logger.info('Adding signal %s'%(sigval * xSecScale * xSecMod))
    #            logger.info('Adding signal %s'%(sigval * xSecScale ))
                
                if False: #sigval>0:
                    if year == 2016:
                        c.specifyUncertainty('Lumi'+postfix, binname, signalName, 1.025 )
                    else:
                        c.specifyUncertainty('Lumi'+postfix, binname, signalName, 1.023 )
                       # signaluncertainties
                    #pu          = 1 + signal.PUSystematic( r, channel, setup).val
                    #jec         = 1 + signal.JECSystematic( r, channel, setup).val
                    #jer         = 1 + signal.JERSystematic( r, channel, setup).val
                    #btag_heavy  = 1 + signal.btaggingSFbSystematic(r, channel, setup).val
                    #btag_light  = 1 + signal.btaggingSFlSystematic(r, channel, setup).val
                    #trigger     = 1 + signal.triggerSystematic(r, channel, setup).val
                    #leptonSFSyst= 1 + signal.leptonSFSystematic(r, channel, setup).val
                    #leptonReco  = 1 + signal.leptonTrackingSystematic(r, channel, setup).val
                    #eleSFStat   = 1 + signal.eleSFSystematic(r, channel, setup).val
                    #muSFStat    = 1 + signal.muSFSystematic(r, channel, setup).val 
    
                    #if sig.sigma/sig.val < 0.05:
                    #    uncertainties['PU']         += [pu]
                    #    uncertainties['JEC']        += [jec]
                    #    uncertainties['btag_heavy'] += [btag_heavy]
                    #    uncertainties['btag_light'] += [btag_light]
                    #    uncertainties['trigger']    += [trigger]
                    #    uncertainties['leptonSFSyst']   += [leptonSFSyst]
    
                    #c.specifyUncertainty('PU',                  binname, signalName, pu)
                    #c.specifyUncertainty('JEC',                 binname, signalName, jec)
                    #c.specifyUncertainty('JER',                 binname, signalName, jer)
                    #c.specifyUncertainty('btag_heavy'+postfix,  binname, signalName, btag_heavy)
                    #c.specifyUncertainty('btag_light'+postfix,  binname, signalName, btag_light)
                    #c.specifyUncertainty('trigger'+postfix,     binname, signalName, trigger)
                    #c.specifyUncertainty('leptonSFSyst',        binname, signalName, leptonSFSyst)
                    #c.specifyUncertainty('leptonTracking',      binname, signalName, leptonReco)
                    #c.specifyUncertainty('eleSFStat'+postfix,   binname, signalName, eleSFStat)
    
                    uname = 'Stat_'+binname+'_%s'%signalName+postfix
                    c.addUncertainty(uname, 'lnN')
                    c.specifyUncertainty(uname, binname, signalName, 1  ) #+ sigdev/sigval )
                else:
                    uname = 'Stat_'+binname+'_%s'%signalName+postfix
                    c.addUncertainty(uname, 'lnN')
                    c.specifyUncertainty(uname, binname, signalName, 1 )
    
                
            #c.addUncertainty('Lumi'+postfix, 'lnN')
            #c.specifyFlatUncertainty('Lumi'+postfix, 1.026)
            cardFileName = c.writeToFile(cardFileName)
        else:
            logger.info("File %s found. Reusing.",cardFileName)
        
        nll = c.calcNLL()

        txtfilename = "ttz2lsc.txt"
        outfile = open(txtfilename, 'a')
        if not "LSTM" in filename : outfile.write(str(args.lumi)+' '+str(nll['nll']+nll['nll0']))
        else  : outfile.write(' '+str(nll['nll']+nll['nll0'])+'\n')
        outfile.close()
    
    sample = 'ttZ_2l_sc'
    if args.lumi : sample = sample +'_lumi' + str(args.lumi).replace('.','p')
    if not "LSTM" in filename : sample = sample
    else : sample = sample + '_lstm'

    wrapper(sample) 
