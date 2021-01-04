import ROOT
from math import *
import copy
from TMB.Tools.helpers import deltaPhi

#X -> (V->v1,v2) (V ->v3->v4) 

def VV_angles( v1, v2, v3, v4, debug=False):
    ''' Assume V1=v1+v2 and V2=v3+v4 are from resonances and compute the angles phi, theta1, and theta2 defined in Fig. 1 in https://arxiv.org/pdf/1708.07823.pdf
    Permutating v1,v2 or v3,v4 reflects the phi->pi-phi ambiguity.
    '''

    # make a copy of lab frame four momenta
    v1 = copy.deepcopy(v1)
    v2 = copy.deepcopy(v2)
    v3 = copy.deepcopy(v3)
    v4 = copy.deepcopy(v4)

    if debug:
        print "Input 4-fectors"
        v1.Print()
        v2.Print()
        v3.Print()
        v4.Print()
        
    # construct & boost to the cms
    cms = v1+v2+v3+v4
    if debug:
        print "c.m.s."
        cms.Print() 
    boostToCMS = -cms.BoostVector()
    v1.Boost(boostToCMS)
    v2.Boost(boostToCMS)
    v3.Boost(boostToCMS)
    v4.Boost(boostToCMS)

    r = ROOT.TLorentzVector()
    r.SetPtEtaPhiM( cms.Pt(), cms.Eta(), cms.Phi(), 0 )
    r.Boost(boostToCMS)
    cms.Boost(boostToCMS)
    if debug:
        print "After boost to c.m.s.: 4-vectors"
        v1.Print()
        v2.Print()
        v3.Print()
        v4.Print()
        print "r (light-like vector in direction of cms)"
        r.Print() 
        print "cms after boost to c.m.s."
        cms.Print()

    # V momenta in the cms (back to back)
    V1 = v1+v2
    V2 = v3+v4
    if debug:
        print "Print V1 and V2 4-vactors (back-to-back)"
        V1.Print()
        V2.Print()

    # rotate such that V1 points to +z
    z_axis   = ROOT.TVector3(0,0,1)
    nV1      = V1.Vect()
    nV1.SetMag(1)
    angle = z_axis.Angle(nV1)
    axis  = z_axis.Cross(nV1)
    for vec in [r, V1, V2, v1, v2, v3, v4 ]:
        vec.Rotate(-angle, axis)
    if debug:
        print "After rotation of V1 to +z:"
        print "V1, V2"
        V1.Print()
        V2.Print()
        print "r"
        r.Print()
        print "v1-v4"
        v1.Print()
        v2.Print()
        v3.Print()
        v4.Print()

    # rotate such that r_xy  points to +x
    x_axis = ROOT.TVector3(1,0,0)
    n_r_T  = ROOT.TVector3()
    n_r_T.SetPtEtaPhi(1,0, r.Phi() )
    if debug:
        print "transverse unit r 3-vector:"
        n_r_T.Print()
    angle = x_axis.Angle(n_r_T)
    axis  = x_axis.Cross(n_r_T)
    for vec in [r, V1, V2, v1, v2, v3, v4 ]:
        vec.Rotate(-angle, axis)
    if debug:
        print "After rotation of r into +x:"
        print "V1, V2"
        V1.Print()
        V2.Print()
        print "r"
        r.Print()
        print "v1-v4"
        v1.Print()
        v2.Print()
        v3.Print()
        v4.Print()

    # rotations complete, let's compute angles
    res = {
        'theta':r.Vect().Angle(z_axis), 
        'phi1':      v1.Phi(), 
        'phi2':     -v3.Phi(),
        }

    try:
        res['deltaPhi'] = (v1.Phi() - v3.Phi())%(2*pi) 
    except ZeroDivisionError:
        res['deltaPhi'] = float('nan')

    ## V unit vectors (back to back)

    # boost v1, v2 to V1 frame
    v1.Boost( -V1.BoostVector() )
    v2.Boost( -V1.BoostVector() )
    if debug:
        print "After boost into V1 cms: v1,v2"
        v1.Print()
        v2.Print()
         
    res['theta_V1'] = z_axis.Angle(v1.Vect())

    v3.Boost( -V2.BoostVector() )
    v4.Boost( -V2.BoostVector() )
    if debug:
        print "After boost into V2 cms: v3,v4"
        v3.Print()
        v4.Print()
    res['theta_V2'] = (-z_axis).Angle(v3.Vect()) 

    if debug:
        print
    return res

if __name__=='__main__':
    v1 = ROOT.TLorentzVector(100,10,70, sqrt(2*100**2+80**2))
    v2 = ROOT.TLorentzVector(80,-100,70, sqrt(2*100**2+80**2))

    v3 = ROOT.TLorentzVector(-100,10,50, sqrt(2*100**2+80**2))
    v4 = ROOT.TLorentzVector(90,100,-70, sqrt(2*100**2+80**2))

    res = VV_angles( v1, v2, v3, v4)

