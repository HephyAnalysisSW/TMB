#!/usr/bin/env python
# Standard imports

import  Analysis.Tools.syncer as syncer
from    RootTools.core.standard import *
from    TMB.Tools.user import plot_directory

import  ROOT
from    math import exp, pi
import  os
import  array

import  numpy as np
np.random.seed(0)

# Definition of the model dsigma/dx = const*(exp(-x/x1)+theta exp(-x/x2))^2
x1      = 1
x2      = 2
theta0  = 0 
#xmax    = 10. 
#texX  = "x"
n       = 1

x12m = 1./(1./x1-1./x2)

# sample from numpy 
def get_sampled_dataset( n_events ):
    # draw events with the larger exponent
    beta_sample = max([x1, x2]) # sample slowest decaying exponential
    #beta_sample = x1 # sample SM
    # correct xsec if beta_sample is not x1
    xsec_corr = beta_sample/x1
    # sample np exponential
    #features = np.random.exponential(beta_sample, n_events)
    features = np.column_stack( (np.random.exponential(beta_sample, n_events), np.random.uniform(-pi,pi,n_events)) )

    # we scale down the weights to the SM value exp(-x/x1) whether or not we sampled ad x1 or x2
    weights  = {tuple():   xsec_corr*np.exp(-features[:,0]*( 1./x1 - 1./beta_sample)),
                (0,)   : 2*xsec_corr*np.exp(-features[:,0]*((1./x1+1./x2)/2.-1./beta_sample))*np.cos(n*features[:,1]),
                (0,0)  : 2*xsec_corr*np.exp(-features[:,0]*( 1./x2 - 1./beta_sample )),
        }
    return features, weights

###############
## Plot model #
###############
#
#colors = [ROOT.kBlack, ROOT.kBlue, ROOT.kGreen, ROOT.kMagenta, ROOT.kCyan, ROOT.kRed]
#plot_features = [ 
#    {'name':'x',   'texX':'x',    'binning':[50,0,10],   },
#    {'name':'phi', 'texX':'#phi', 'binning':[50,-pi,pi], },
#    ]
#
#h = {}
#thetas = [.05, .1, .15, .2]
#n_events_plot = 10**5
#features, weights = get_sampled_dataset(n_events_plot)
#for i_theta, theta in enumerate(thetas):
#    name = ("theta_%3.2f"%theta ).replace('.','p').replace('-','m') 
#    h[theta] = {}
#    for i_feature, feature in enumerate(plot_features):
#        h[theta][i_feature]              = ROOT.TH1F(name+'_'+feature['name'], name+'_'+feature['name'], *plot_features[i_feature]['binning'] )
#        h[theta][i_feature].style        = styles.lineStyle( colors[i_theta], width=2, dashed=False )
#        h[theta][i_feature].legendText   = "#theta = %3.2f"%theta 
#
#for i_event, event in enumerate(features):
#    if i_event%10**4==0: print "At %i/%i"%( i_event, n_events_plot )
#    for i_feature, feature in enumerate(plot_features):
#        for i_theta, theta in enumerate(thetas):
#            h[theta][i_feature].Fill(event[i_feature], weights[()][i_event] + theta*weights[(0,)][i_event] + 0.5*theta**2*weights[(0,0)][i_event])
#
#for i_feature, feature in enumerate(plot_features):
#    histos = [[h[theta][i_feature]] for theta in thetas]
#    plot   = Plot.fromHisto( feature['name'],  histos, texX=feature['texX'], texY="a.u." )
#
#    for log in [True, False]:
#
#        # Add subdirectory for lin/log plots
#        plot_directory_ = os.path.join( plot_directory, "newman", "log" if log else "lin" )
#        plotting.draw( plot,
#                       plot_directory = plot_directory_,
#                       logX = False, logY = log, sorting = False,
#                       yRange = "auto",
#                       ratio = None,
##                       drawObjects = drawObjects( lumi, offset=titleOffset ),
#                        legend=(0.2,0.7,0.9,0.9),
#                       #histModifications = histModifications,
#                       copyIndexPHP = True,
#                       )

def make_TH1F( h ):
    # remove infs from thresholds
    vals, thrs = h
    #thrs[0]  = thrs[1] - (thrs[2]-thrs[1])
    #thrs[-1] = thrs[-2] + (thrs[-2]-thrs[-3])
    histo = ROOT.TH1F("h","h",len(thrs)-1,array.array('d', thrs))
    for i_v, v in enumerate(vals):
        histo.SetBinContent(i_v+1, v)
    return histo

# we simulate a large number of events
dataset_size      = 10**6
features, weights = get_sampled_dataset(dataset_size)
event_indices     = np.arange(dataset_size)

# normalize to SM xsec
sigma_tot_sm      = 1
#weights[(0,)]     = weights[(0,)]*sigma_tot_sm/np.sum(weights[tuple()])
#weights[(0,0)]    = weights[(0,0)]*sigma_tot_sm/np.sum(weights[tuple()])
#weights[tuple()]  = weights[tuple()]*sigma_tot_sm/np.sum(weights[tuple()])

# total xsec ratio
sigma_tot_ratio_lin  = np.sum(weights[(0,)])/np.sum(weights[tuple()])
sigma_tot_ratio_quad = np.sum(weights[(0,0)])/np.sum(weights[tuple()])

# xsec ratio
def sigma_tot_ratio( theta, lin=False):
    return 1 + theta*sigma_tot_ratio_lin + ( 0.5*theta**2*sigma_tot_ratio_quad if not lin else 0 )
# total xsec
def sigma_tot( theta, lin=False):
    return sigma_tot_sm*sigma_tot_ratio( theta, lin=lin ) 


# expansion: (1 + theta*q_event_lin)**2 + (theta*q_event_quad)**2
q_event_lin  = np.cos(n*features[:,1])*np.exp(features[:,0]/(2*x12m))
q_event_quad = np.sin(n*features[:,1])*np.exp(features[:,0]/(2*x12m))

#def wrapper( i_toy, toys, plot = plot, test_statistic = "total", extended = False, verbose = True):
#for i_toy in range(5):


def make_toys( yield_per_toy, theta, n_toys):
    weights_ =  weights[tuple()]+theta*weights[(0,)]+.5*theta**2*weights[(0,0)]
    biased_sample = np.random.choice( event_indices, size=dataset_size,  p = weights_/np.sum(weights_) )

    return np.array( [ np.random.choice( biased_sample, size=n_observed ) for n_observed in np.random.poisson(yield_per_toy, n_toys) ])

#for i_toy, toy in enumerate( toys_SM[:5] ):
def wrapper( i_toy, toys, plot = True, test_statistic = "total", extended = False, verbose = True):
    q_theta_given_theta_1mCL = {}

    import  Analysis.Tools.syncer as syncer
    theta_current    = theta_UL_start 

    #print "Toy", i_toy
    toy = toys[i_toy]
    i_iter    = 0
    while True:

        # This is the second toy -> again, we guess theta
        if i_iter==1:
            delta_q_previous = q_theta_given_theta_1mCL[theta_current] - q_theta_SM 
            theta_previous   = theta_UL_start
            theta_current    = theta_UL_start/2. 

        # We are at the third iteration. Time to update.
        elif i_iter>1:

            delta_q_current = q_theta_given_theta_1mCL[theta_current] - q_theta_SM
             
            #print q_theta_SM,  q_theta_given_theta_1mCL[theta_current]
            if abs(theta_current - theta_previous)<=tolerance:
                # Converged!
                if verbose: 
                    print "Toy", i_toy, "done."
                    print "theta_current",  theta_current, "delta_q_current", delta_q_current, "i_toy", i_toy, "i_iter", i_iter
                    print
                break
                #return [ theta_current, delta_q_current, i_toy, i_iter]
            elif i_iter>=max_iter:
                if verbose: 
                    print "Toy", i_toy, "max_iter %i reached" % max_iter
                    print "theta_current", theta_current, "theta_previous", theta_previous, "delta_q_current", delta_q_current
                    print
                break

            # Newton step
            k = (theta_current-theta_previous)/(delta_q_current-delta_q_previous)
            if k>2:
                if verbose: print "k=%3.2f too large, set it to 2." % k
                k=2
            if k<-2:
                if verbose: print "k=%3.2f too small, set it to -2." % k
                k=-2
            theta_current, theta_previous, delta_q_previous = theta_current - delta_q_current*k, theta_current, delta_q_current 
            # If the predicted value is negative, cool down and half the distance (to zero)
            if theta_current<0:
                if verbose: print "Predicted negative, cooling down."
                theta_current = theta_previous/2.

        # compute value of q_theta test statistic for all events 
        #q_event = 1/theta_current * np.log( (1 + theta_current*q_event_lin)**2 )
        if test_statistic == "quad":
            q_event = np.log( theta_current**2*(q_event_lin**2+q_event_quad**2) )
        elif test_statistic == "lin":
            q_event = np.log( (1 + theta_current*q_event_lin)**2 )
        elif test_statistic == "total":
            q_event = np.log( (1 + theta_current*q_event_lin)**2 + (theta_current*q_event_quad)**2)
        else:
            raise RuntimeError( "Unknwon test statistc %s" % test_statistic )

        log_sigma_tot_ratio_subtraction = 0 #np.log(sigma_tot_ratio(theta_current)) if not extended else 0

        q_theta_SM          = np.sum( q_event[toy] - log_sigma_tot_ratio_subtraction )

        # compute 1-CL quantile of the q_theta test statistic under the theta hypothesis
        if not q_theta_given_theta_1mCL.has_key(theta_current) or plot:
            # compute distribution of test statistic q_theta under theta 
            q_theta_given_theta = np.array([np.sum( q_event[toy_] - log_sigma_tot_ratio_subtraction ) for toy_ in make_toys( lumi*sigma_tot(theta_current), theta_current, n_toys ) ])
            q_theta_given_theta_1mCL[theta_current] = np.quantile( q_theta_given_theta, 1-CL)
            if plot:
                tex = ROOT.TLatex()
                tex.SetNDC()
                tex.SetTextSize(0.04)
                tex.SetTextAlign(11) # align right

                # Text on the plots
                lines = [ 
                          (0.25, 0.88, "#color[4]{%i%% qu. q_{BSM} = %3.2f}" % ( 100*(1-CL), q_theta_given_theta_1mCL[theta_current]) ),
                          (0.25, 0.83, "#color[2]{q_{SM} = %3.2f}" % ( q_theta_SM ) ),
                          (0.25, 0.78, "#theta_{current} = %5.4f" % theta_current ),
                        ]
                drawObjects = [ tex.DrawLatex(*line) for line in lines ]

                np_histo      = np.histogram(q_theta_given_theta, 100)
                histo         = make_TH1F(np_histo)
                quantile_line = ROOT.TLine(q_theta_given_theta_1mCL[theta_current], 0, q_theta_given_theta_1mCL[theta_current], histo.GetMaximum())
                quantile_line.SetLineColor(ROOT.kRed)
                histo.style = styles.lineStyle( ROOT.kRed )
                histo.legendText = "#color[2]{q_{BSM}}" 

                q_theta_SM_line = ROOT.TLine(q_theta_SM, 0, q_theta_SM, histo.GetMaximum())
                q_theta_SM_line.SetLineColor(ROOT.kBlue)

                q_theta_given_SM = np.array([np.sum( q_event[toy_]-log_sigma_tot_ratio_subtraction ) for toy_ in toys_SM])
                histo_SM         = make_TH1F(np.histogram(q_theta_given_SM, bins=np_histo[1]))
                histo_SM.legendText = "#color[4]{q_{SM}}" 
                histo_SM.style = styles.lineStyle( ROOT.kBlue )

                plot = Plot.fromHisto( "itoy_%03i_iter_%02i"%( i_toy, i_iter ),[[histo], [histo_SM]], texX = "q", texY = "Entries" )
                plotting.draw( plot,
                    plot_directory = os.path.join( plot_directory, "newman_new2_%3.2f"% ( lumi*sigma_tot(theta_SM)) ),
                    #ratio          = {'yRange':(0.6,1.4)} if len(plot.stack)>=2 else None,
                    logX = False, sorting = False,
                    legend         = (0.65,0.78,0.9,0.9),
                    drawObjects    =  [quantile_line, q_theta_SM_line] + drawObjects,
                    copyIndexPHP   = True,
                    extensions     = ["png"], 
                  )            
        i_iter +=1
    return [ theta_current, delta_q_current, i_toy, i_iter]


UL                       = {}

#UL = [ wrapper(i_toy, toys_SM, plot=True) for i_toy in range(100) ]


from multiprocessing import Pool

theta_SM    = 0

CL          = 0.95
lumi        = 13.70
n_toys      = 50000
theta_UL_start = 1
tolerance      = 0.001
max_iter       = 50
test_statistic = "total"

plot           = True
extended       = True
verbose        = True

n_toys_UL      = 20
UL             = {}

for lumi in [ 13.70 , 2*13.70, 5*13.70, 10*13.70, 20*13.70, 50*13.70 ]: 

    toys_SM = make_toys( lumi*sigma_tot(theta_SM), theta_SM, n_toys )

    # precompute toy selection: n_toys toys of n_observed events that are selected from a Poissonian with mean n_expected
    toys            = np.array( [ np.random.choice( event_indices, size=n_observed ) for n_observed in np.random.poisson(n_expected, n_toys) ])

    UL[lumi] = {}
    for test_statistic in ["total", "lin", "quad"]:

        #wrapper_(0)

        class wrapper_(object):
            def __init__(self, *args, **kwargs ):
                self.args = args
                self.kwargs = kwargs
            def __call__(self, i_toy):
                return wrapper( i_toy, *self.args, **self.kwargs )

        p  = Pool(20)

        UL_ = np.array(p.map( wrapper_(verbose=False, toys=toys, plot=False, test_statistic=test_statistic), range(n_toys_UL)))
        UL[lumi][test_statistic] = UL_

        print "lumi %3.2f test_statistic %s median expected %i%% CL UL: %5.4f coverage %3.2f%%" % ( lumi, test_statistic, 100*CL, np.quantile( UL_[:,0], 0.5 ), 100*np.count_nonzero(UL_[:,0]>tolerance)/float(len(UL_)) ) 
