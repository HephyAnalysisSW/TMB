#!/usr/bin/env python
# Standard imports

import Analysis.Tools.syncer as syncer
import numpy as np
import random
import ROOT
from math import exp, pi
import os

# Definition of the model dsigma/dx = const*(exp(-x/x1)+theta exp(-x/x2))^2
x1 = 1.
x2 = 2.
theta0 = 0 
xmax  = 10. 
#texX  = "x"
const = 1
n     = 2

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

dataset_size = 10**6
features, weights = get_sampled_dataset(dataset_size)

lumi = 137
xsec = .00001
n_expected  = lumi*xsec*5000
print "n_expected", n_expected

n_toys = 5000
#n_toys = 50

def make_test_statistic( theta, mode="total"):
    x12m = 1./(1./x1-1./x2)
    if mode=="lin":
        return lambda weights, array, theta=theta: np.sum( weights[array]* np.log( (1 + theta*np.cos(n*features[array][:,1])*np.exp(features[array][:,0]/(2*x12m)) )**2 ))
    elif mode=="quad":
        return lambda weights, array, theta=theta: np.sum( weights[array]* np.log( (theta*np.sin(n*features[array][:,1])*np.exp(features[array][:,0]/(2*x12m)))**2  ))
    else:
        return lambda weights, array, theta=theta: np.sum( weights[array]* np.log( (1 + theta*np.cos(n*features[array][:,1])*np.exp(features[array][:,0]/(2*x12m)) )**2 + (theta*np.sin(n*features[array][:,1])*np.exp(features[array][:,0]/(2*x12m)))**2 ))

# tGraphs 

threshold   = 0.95
theta_vals  = [.5, .75, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0]
modes = ["lin", "quad", "total"]

color       = {'lin':ROOT.kGreen, 'quad':ROOT.kRed, 'total':ROOT.kBlack }   

tGraph = {}
for mode in modes: 

    tGraph[mode] = ROOT.TGraph( len(theta_vals) )
    tGraph[mode].SetLineColor( color[mode] )
    tGraph[mode].SetMarkerColor( color[mode] )
    tGraph[mode].SetMarkerStyle( 0 )
    tGraph[mode].SetLineWidth(2)

    for i_theta, theta in enumerate( theta_vals ):
        q = make_test_statistic( theta, mode=mode )
        q_array_SM = np.array( [ q(weights[tuple()], np.random.choice( np.arange(dataset_size), size=n_observed )) for n_observed in np.random.poisson(n_expected, n_toys) ])
        k = np.quantile(q_array_SM, threshold)

        size = np.count_nonzero(q_array_SM>k)/float(len(q_array_SM))

        q_array_BSM = np.array([q(weights[tuple()] + theta*weights[(0,)]+.5*theta**2*weights[(0,0)], np.random.choice( np.arange(dataset_size), size=n_observed )) for n_observed in np.random.poisson(n_expected, n_toys) ])
        power = np.count_nonzero(q_array_BSM>k)/float(len(q_array_BSM))
        print "Mode?", mode, "theta", theta, "threshold", threshold, "k=%5.4f"% k, "size", size, "power", power
        tGraph[mode].SetPoint( i_theta, theta, power ) 

c1 = ROOT.TCanvas()
ROOT.gStyle.SetOptStat(0)
c1.SetTitle("")

l = ROOT.TLegend(0.7, 0.6, 0.9, 0.8)
l.SetFillStyle(0)
l.SetShadowColor(ROOT.kWhite)
l.SetBorderSize(0)

first = True
for mode in modes:
    l.AddEntry(  tGraph[mode], mode )
    tGraph[mode].Draw("AL" if first else "L")
    tGraph[mode].SetTitle("")
    tGraph[mode].GetXaxis().SetTitle("#theta")
    tGraph[mode].GetYaxis().SetTitle("power")

    first = False

l.Draw()
plot_directory = "/mnt/hephy/cms/robert.schoefbeck/www/"

c1.RedrawAxis()
c1.Print(os.path.join(plot_directory, "etc/power.png"))
c1.Print(os.path.join(plot_directory, "etc/power.pdf"))
c1.Print(os.path.join(plot_directory, "etc/power.root"))
syncer.sync()
