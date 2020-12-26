import ROOT
from math import *
import copy

#X -> (V->v1,v2) (V ->v3->v4) 

def VV_angles( v1, v2, v3, v4):
    ''' Assume V1=v1+v2 and V2=v3+v4 are from resonances and compute the angles phi, theta1, and theta2 defined in Fig. 1 in https://arxiv.org/pdf/2012.11631.pdf
    Permutating v1/v2 or v3/v4 reflects the phi->pi-phi ambiguity.
    '''

    v1 = copy.deepcopy(v1)
    v2 = copy.deepcopy(v2)
    v3 = copy.deepcopy(v3)
    v4 = copy.deepcopy(v4)

    # construct and boost to the cms
    cms = v1+v2+v3+v4
    boostToCMS = -cms.BoostVector()
    v1.Boost(boostToCMS)
    v2.Boost(boostToCMS)
    v3.Boost(boostToCMS)
    v4.Boost(boostToCMS)

    # V momenta in the cms (back to back)
    V1 = v1+v2
    V2 = v3+v4

    # V unit vectors (back to back)
    nV1 = V1.Vect() 
    nV1.SetMag(1)
    nV1.SetMag(nV1*v1.Vect())
    nv1_sub = v1.Vect() - nV1 

    nV2 = V2.Vect() 
    nV2.SetMag(1)
    nV2.SetMag(nV2*v3.Vect())
    nv2_sub = v3.Vect() - nV2

    cos_phi = (nv1_sub*nv2_sub)/(nv1_sub.Mag()*nv2_sub.Mag())

    # boost v1, v2 to V1 frame
    v1.Boost( -V1.BoostVector() )
    v2.Boost( -V1.BoostVector() )
    cos_theta_V1 = v1.Vect()*nV1/(v1.Vect().Mag()*nV1.Mag())

    v3.Boost( -V2.BoostVector() )
    v4.Boost( -V2.BoostVector() )
    cos_theta_V2 = v3.Vect()*nV2/(v3.Vect().Mag()*nV2.Mag())

    return cos_phi, cos_theta_V1, cos_theta_V2

if __name__=='__main__':
    v1 = ROOT.TLorentzVector(100,10,70, sqrt(2*100**2+80**2))
    v2 = ROOT.TLorentzVector(80,-100,70, sqrt(2*100**2+80**2))

    v3 = ROOT.TLorentzVector(-100,10,50, sqrt(2*100**2+80**2))
    v4 = ROOT.TLorentzVector(90,100,-70, sqrt(2*100**2+80**2))

    VV_angles( v1, v2, v3, v4)

