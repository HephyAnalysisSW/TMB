#!/usr/bin/env python
# Standard imports

import  Analysis.Tools.syncer as syncer
from    RootTools.core.standard import *
from    TMB.Tools.user import plot_directory
import  numpy as np
import  random
import  ROOT
from    math import exp, pi
import  os
import  array

# Definition of the model dsigma/dx = const*(exp(-x/x1)+theta exp(-x/x2))^2
x1      = 1.
x2      = 2.
theta0  = 0 
xmax    = 10. 
#texX  = "x"
const   = 1
n       = 2

x12m = 1./(1./x1-1./x2)

# sample from numpy 
def get_sampled_dataset( n_events ):
    # draw events with the larger exponent
    beta_sample = max([x1, x2])
    # correct xsec if beta_sample is not x1
    xsec_corr = beta_sample/x1
    # sample np exponential
    #features = np.random.exponential(beta_sample, n_events)
    features = np.array( [np.random.exponential(beta_sample, n_events), np.random.uniform(-pi,pi,n_events)] ).reshape(-1,2)

    # we scale down the weights to the SM value exp(-x/x1) whether or not we sampled ad x1 or x2
    weights  = {tuple():   xsec_corr*const*np.exp(-features[:,0]*( 1./x1 - 1./beta_sample)),
                (0,)   : 2*xsec_corr*const*np.exp(-features[:,0]*((1./x1+1./x2)/2.-1./beta_sample))*np.cos(n*features[:,1]),
                (0,0)  : 2*xsec_corr*const*np.exp(-features[:,0]*( 1./x2 - 1./beta_sample )),
        }
    return features, weights

def make_TH1F( h ):
    # remove infs from thresholds
    vals, thrs = h
    #thrs[0]  = thrs[1] - (thrs[2]-thrs[1])
    #thrs[-1] = thrs[-2] + (thrs[-2]-thrs[-3])
    histo = ROOT.TH1F("h","h",len(thrs)-1,array.array('d', thrs))
    for i_v, v in enumerate(vals):
        histo.SetBinContent(i_v+1, v)
    return histo

dataset_size = 10**6
features, weights = get_sampled_dataset(dataset_size)

lumi = 137
xsec = .0001
n_expected  = lumi*xsec*5000
print "n_expected", n_expected

n_toys = 5000
#n_toys = 50

# precompute toy selection: n_toys toys of n_observed events that are selected from a Poissonian with mean n_expected
event_indices   = np.arange(dataset_size)
toys            = np.array( [ np.random.choice( event_indices, size=n_observed ) for n_observed in np.random.poisson(n_expected, 10*n_toys) ])
#log_event_statistic = np.log( (1 + theta*np.cos(n*features[:,1])*np.exp(features[:,0]/(2*x12m)) )**2 + (theta*np.sin(n*features[:,1])*np.exp(features[:,0]/(2*x12m)))**2 )

#def make_test_statistic( theta, mode="total"):
#    if mode=="lin":
#        return lambda weights, array, theta=theta: np.sum( weights[array]* np.log( (1 + theta*np.cos(n*features[array][:,1])*np.exp(features[array][:,0]/(2*x12m)) )**2 ))
#    elif mode=="quad":
#        return lambda weights, array, theta=theta: np.sum( weights[array]* np.log( (theta*np.sin(n*features[array][:,1])*np.exp(features[array][:,0]/(2*x12m)))**2  ))
#    else:
#        return lambda weights, array, theta=theta: np.sum( weights[array]* np.log( (1 + theta*np.cos(n*features[array][:,1])*np.exp(features[array][:,0]/(2*x12m)) )**2 + (theta*np.sin(n*features[array][:,1])*np.exp(features[array][:,0]/(2*x12m)))**2 ))

# tGraphs 
plot = False
tex = ROOT.TLatex()
tex.SetNDC()
tex.SetTextSize(0.04)
tex.SetTextAlign(11) # align right

CL   = 0.95

theta_UL_start = .05
xtol           = 0.0001

q_log_event_lin  = np.cos(n*features[:,1])*np.exp(features[:,0]/(2*x12m))
q_log_event_quad = np.sin(n*features[:,1])*np.exp(features[:,0]/(2*x12m)) 
#q_log_event_quad = np.zeros(len(features))

q_log_event              = {}
w_event_BSM              = {}
q_theta_given_theta_1mCL = {}
UL                       = {}
i_iter = {i_toy:-1 for i_toy in range(len(toys))}

import scipy.optimize
def wrapper( i_toy, plot = True ):

    import  Analysis.Tools.syncer as syncer
    #print "Toy", i_toy
    toy = toys[i_toy]
    def delta_quantiles(theta_current, plot=plot):
        i_iter[i_toy]+=1
        # compute value of q_theta test statistic under SM for each event
        if not q_log_event.has_key(theta_current):
            q_log_event[theta_current]  = np.log( (1 + theta_current*q_log_event_lin)**2 + (theta_current*q_log_event_quad)**2 )
        q_theta_SM = np.sum( weights[tuple()][toy]*q_log_event[theta_current][toy] )

        if not w_event_BSM.has_key( theta_current ):
            w_event_BSM[theta_current] = weights[tuple()] + theta_current*weights[(0,)]+.5*theta_current**2*weights[(0,0)]

        # compute 1-CL quantile of the q_theta test statistic under the theta hypothesis
        if not q_theta_given_theta_1mCL.has_key(theta_current) or plot:
            # compute distribution of test statistic q_theta under theta 
            q_theta_given_theta = np.array([np.sum(w_event_BSM[theta_current][toy_]*q_log_event[theta_current][toy_]) for toy_ in toys])
            q_theta_given_theta_1mCL[theta_current] = np.quantile( q_theta_given_theta, 1-CL)
            if plot:
                # Text on the plots
                lines = [ 
                          (0.25, 0.88, "#color[4]{%i%% qu. q_{BSM} = %3.2f}" % ( 100*(1-CL), q_theta_given_theta_1mCL[theta_current]) ),
                          (0.25, 0.83, "#color[2]{q_{SM} = %3.2f}" % ( q_theta_SM ) ),
                          (0.25, 0.78, "#theta_{current} = %5.4f" % theta_current ),
                        ]
                drawObjects = [ tex.DrawLatex(*line) for line in lines ]

                histo = make_TH1F(np.histogram(q_theta_given_theta, 100))
                quantile_line = ROOT.TLine(q_theta_given_theta_1mCL[theta_current], 0, q_theta_given_theta_1mCL[theta_current], histo.GetMaximum())
                quantile_line.SetLineColor(ROOT.kRed)
                histo.style = styles.lineStyle( ROOT.kRed )
                histo.legendText = "#color[2]{q_{BSM}}" 

                q_theta_SM_line = ROOT.TLine(q_theta_SM, 0, q_theta_SM, histo.GetMaximum())
                q_theta_SM_line.SetLineColor(ROOT.kBlue)

                q_theta_given_SM = np.array([np.sum(weights[tuple()][toy_]*q_log_event[theta_current][toy_]) for toy_ in toys])
                histo_SM         = make_TH1F(np.histogram(q_theta_given_SM, 100))
                histo_SM.legendText = "#color[4]{q_{SM}}" 
                histo_SM.style = styles.lineStyle( ROOT.kBlue )

                plot = Plot.fromHisto( "itoy_%03i_iter_%02i"%( i_toy, i_iter[i_toy] ),[[histo], [histo_SM]], texX = "q", texY = "Entries" )
                plotting.draw( plot,
                    plot_directory = os.path.join( plot_directory, "newman_%3.2f"%n_expected ),
                    #ratio          = {'yRange':(0.6,1.4)} if len(plot.stack)>=2 else None,
                    logX = False, sorting = False,
                    legend         = (0.65,0.78,0.9,0.9),
                    drawObjects    =  [quantile_line, q_theta_SM_line] + drawObjects,
                    copyIndexPHP   = True,
                    extensions     = ["png"], 
                  )            

        return q_theta_given_theta_1mCL[theta_current] - q_theta_SM
         
        ## start Newton by guessing the -1 st step. It will converge at half of theta_current
        #if delta_q_previous is None:
        #    theta_previous   = 1.5*theta_current
        #    delta_q_previous = 5*delta_q_current

        ##print q_theta_SM,  q_theta_given_theta_1mCL[theta_current]

        ## Newton step
        #theta_current, theta_previous, delta_q_previous = theta_current - delta_q_current*(theta_current-theta_previous)/(delta_q_current-delta_q_previous), theta_current, delta_q_current 

    result = scipy.optimize.minimize_scalar( delta_quantiles, method='bounded', bounds=(0, 10), options={'xtol':xtol})

        #print "theta_current", theta_current, "theta_previous", theta_previous, "delta_q_previous", delta_q_previous
    print "Toy", i_toy, "done"

    return result 

UL = [ wrapper(i_toy, plot=True) for i_toy in range(1) ]
print "median expected %i%% CL UL: %5.4f" % ( CL, np.quantile( UL, 0.5 ) ) 

#from multiprocessing import Pool
#p  = Pool(20)
#UL = p.map(wrapper, range(100))
#print "median expected %i%% CL UL: %5.4f" % ( CL, np.quantile( UL, 0.5 ) ) 

#print "UL[i_toy]",list(reversed(UL[i_toy]))
        
#theta_vals  = [.5, .75, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0]
#modes = ["lin", "quad", "total"]
#
#color       = {'lin':ROOT.kGreen, 'quad':ROOT.kRed, 'total':ROOT.kBlack }   
#
#tGraph = {}
#for mode in modes: 
#
#    tGraph[mode] = ROOT.TGraph( len(theta_vals) )
#    tGraph[mode].SetLineColor( color[mode] )
#    tGraph[mode].SetMarkerColor( color[mode] )
#    tGraph[mode].SetMarkerStyle( 0 )
#    tGraph[mode].SetLineWidth(2)
#
#    for i_theta, theta in enumerate( theta_vals ):
#        q = make_test_statistic( theta, mode=mode )
#        q_array_SM = np.array( [ q(weights[tuple()], np.random.choice( np.arange(dataset_size), size=n_observed )) for n_observed in np.random.poisson(n_expected, n_toys) ])
#    #if lin:
#    #    return lambda weights, array, theta=theta: np.sum( weights[array]* np.log( 1 + 2*theta*np.exp(features[array]/(2*x12m)) ) )
#    #else:
#    #    return lambda weights, array, theta=theta: np.sum( weights[array]* np.log( 1 + 2*theta*np.exp(features[array]/(2*x12m)) + theta**2*np.exp(features[array]/x12m)) )
#        k = np.quantile(q_array_SM, threshold)
#
#        size = np.count_nonzero(q_array_SM>k)/float(len(q_array_SM))
#
#        q_array_BSM = np.array([q(weights[tuple()] + theta*weights[(0,)]+.5*theta**2*weights[(0,0)], np.random.choice( np.arange(dataset_size), size=n_observed )) for n_observed in np.random.poisson(n_expected, n_toys) ])
#        power = np.count_nonzero(q_array_BSM>k)/float(len(q_array_BSM))
#        print "Mode?", mode, "theta", theta, "threshold", threshold, "k=%5.4f"% k, "size", size, "power", power
#        tGraph[mode].SetPoint( i_theta, theta, power ) 
#
#c1 = ROOT.TCanvas()
#ROOT.gStyle.SetOptStat(0)
#c1.SetTitle("")
#
#l = ROOT.TLegend(0.7, 0.6, 0.9, 0.8)
#l.SetFillStyle(0)
#l.SetShadowColor(ROOT.kWhite)
#l.SetBorderSize(0)
#
#first = True
#for mode in modes:
#    l.AddEntry(  tGraph[mode], mode )
#    tGraph[mode].Draw("AL" if first else "L")
#    tGraph[mode].SetTitle("")
#    tGraph[mode].GetXaxis().SetTitle("#theta")
#    tGraph[mode].GetYaxis().SetTitle("power")
#
#    first = False
#
#l.Draw()
#plot_directory = "/mnt/hephy/cms/robert.schoefbeck/www/"
#
#c1.RedrawAxis()
#c1.Print(os.path.join(plot_directory, "etc/power.png"))
#c1.Print(os.path.join(plot_directory, "etc/power.pdf"))
#c1.Print(os.path.join(plot_directory, "etc/power.root"))
#syncer.sync()
