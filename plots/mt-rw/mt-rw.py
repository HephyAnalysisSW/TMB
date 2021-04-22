# Standard imports
import os
import logging
import ROOT
import array
import Analysis.Tools.syncer

#RootTools
from RootTools.core.standard import *
from TMB.Tools.user import plot_directory

from Analysis.Tools.GenSearch import GenSearch

maxN_files = 1

products = {
    'gp':      {'type':'vector<reco::GenParticle>', 'label':("genParticles")},
    'lhe':{'type':'LHEEventProduct', 'label':("externalLHEProducer")},
    }


weights = ['mt_174p0', 'mt_173p9', 'mt_173p8', 'mt_173p7', 'mt_173p6', 'mt_173p5', 'mt_173p4', 'mt_173p3', 'mt_173p2', 'mt_173p1', 'mt_173p0', 'mt_172p9', 'mt_172p8', 'mt_172p7', 'mt_172p6', 'mt_172p5', 'mt_172p4', 'mt_172p3', 'mt_172p2', 'mt_172p1', 'mt_172p0', 'mt_171p9', 'mt_171p8', 'mt_171p7', 'mt_171p6', 'mt_171p5', 'mt_171p4', 'mt_171p3', 'mt_171p2', 'mt_171p1', 'mt_171p0',] 

ROOT.gStyle.SetPalette(ROOT.kSolar)
wd = {}
for i_weight, weight in enumerate(weights):
    m = float(weight.split('_')[1].replace('p','.'))
    wd[weight] = {
        'm'   : m,
        'h_mt' : ROOT.TH1F('mt_'+weight, 'mt_'+weight, 30, 165, 180),
        'h_mlb' : ROOT.TH1F('mlb_'+weight, 'mlb_'+weight, 30, 50, 250),
        }  

    wd[weight]['h_mt'].SetLineColor(40+i_weight)
    wd[weight]['h_mlb'].SetLineColor(40+i_weight)

tt1l  = FWLiteSample.fromDAS("tt1l" , "/tt01j_1l_LO_MLM/schoef-tt01j_1l_LO_MLM-dad6b6ec8b7b0c353ead8881638d7fd8/USER", instance="phys03", maxN=maxN_files, prefix='root://eos.grid.vbc.ac.at/')
#tt1l  = FWLiteSample.fromDirectory("tt1l" , "/eos/vbc/experiments/cms/store/user/schoef/tt01j_1l_LO_MLM/tt01j_1l_LO_MLM/210421_194721/0000/", maxN = 2)
fwliteReader = tt1l.fwliteReader( products = products )
fwliteReader.start()
counter=0
while fwliteReader.run():
    if counter%100==0:
        print "At",counter
    lhe_weights = fwliteReader.products['lhe'].weights()
    #weights      = []
    #param_points = []
    rw_nominal = float('nan')
    mt_weights = {}
    for weight in lhe_weights:
        if weight.id in ['rwgt_1','dummy']: rw_nominal = weight.wgt

        # Store nominal weight (First position!)
        if weight.id.startswith('mt_'): 
            #print counter, weight.id, weight.wgt/rw_nominal
            mt_weights[weight.id] = weight.wgt

    genSearch = GenSearch(fwliteReader.event.gp)

    gp = fwliteReader.event.gp

    tops = filter( lambda p: abs(p.pdgId())==6 and genSearch.isLast(p), gp )
    for top in tops:
        top_d = genSearch.daughters(genSearch.descend(top))
        W, b = None, None
        if abs(top_d[0].pdgId())==5 and abs(top_d[1].pdgId())==24:
            b, W = top_d
        elif abs(top_d[1].pdgId())==5 and abs(top_d[0].pdgId())==24:
            W, b = top_d
        lep = None
        if W is not None:
            leps = filter(lambda p: abs(p.pdgId()) in [11,13], genSearch.daughters(genSearch.descend(W)) )
            if len(leps)==1:
                lep = leps[0]
        mlb = None
        if lep is not None and b is not None:
            mlb = (lep.p4()+b.p4()).M()
     
        for weight_id, weight in mt_weights.iteritems():
            wd[weight_id]['h_mt'].Fill( top.p4().M(), weight )
            if mlb is not None:
                wd[weight_id]['h_mlb'].Fill( mlb, weight )

    counter+=1

mlb_rw = Plot.fromHisto(name = "mlb_n"+str(len(tt1l.files)), histos = [[wd[w]['h_mlb']] for w in weights], texX = "gen mlb" , texY = "arbitrary" )


nicename = {'mt': 'gen m_{t}', 'mlb':'gen m_{lb}'}

for varname in ['mt', 'mlb']:
    plot  = Plot.fromHisto(name = varname+"_n"+str(len(tt1l.files)), histos = [[wd[w]['h_'+varname]] for w in weights], texX = nicename[varname] , texY = "arbitrary" )

    h_nom = wd['mt_172p5']['h_'+varname].Clone()
    h_nom.SetLineColor(ROOT.kBlue) 
    h_nom.SetLineWidth(2) 
    plot.histos.append( [h_nom] )
    plotting.draw(plot, 
            plot_directory = os.path.join( plot_directory, 'mt_rw'),
            yRange = 'auto',
            legend = None,
            ratio = {'yRange':(0,2), 'histos': [ (i,weights.index('mt_172p5')) for i in range(len(weights))]}, logY = True, logX = False,
    #        drawObjects = [h_nom]
        )
