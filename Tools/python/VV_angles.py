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
        print
        print "Input 4-fectors"
        v1.Print()
        v2.Print()
        v3.Print()
        v4.Print()
        
    # construct & boost to the cms
    cms = v1+v2+v3+v4
    if debug:
        print
        print "c.m.s."
        cms.Print() 
    boostToCMS = -cms.BoostVector()
    v1.Boost(boostToCMS)
    v2.Boost(boostToCMS)
    v3.Boost(boostToCMS)
    v4.Boost(boostToCMS)

    r = ROOT.TLorentzVector()
    #r.SetPtEtaPhiM( cms.Pt(), cms.Eta(), cms.Phi(), 0 ) #wtf
    r.SetPxPyPzE( 0, 0, 6500, 6500 )
    
    r.Boost(boostToCMS)
    cms.Boost(boostToCMS)
    if debug:
        print
        print "After boost to c.m.s.: 4-vectors"
        v1.Print()
        v2.Print()
        v3.Print()
        v4.Print()
        print "r (light-like vector in beam direction)"
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
        print
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
        print
        print "Rotation axis (in z direction)"
        axis.Print()
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
        'Theta':r.Vect().Angle(z_axis), 
        'phi1':      v1.Phi(), 
        'phi2':      -v3.Phi(),
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

# from https://github.com/HephyAnalysisSW/VH/blob/main/Tools/python/helpers.py
from ROOT import TLorentzVector, TVector3
def getTheta(lep1, lep2, H):

    beam = TLorentzVector()

    tmp_lep1, tmp_lep2, tmp_H = TLorentzVector(), TLorentzVector(), TLorentzVector()

    tmp_lep1.SetPtEtaPhiM(lep1.Pt(),lep1.Eta(),lep1.Phi(),lep1.M())
    tmp_lep2.SetPtEtaPhiM(lep2.Pt(),lep2.Eta(),lep2.Phi(),lep2.M())
    tmp_H.SetPtEtaPhiM(H.Pt(),H.Eta(),H.Phi(),H.M())

    if(lep1.Eta()<-10 or lep2.Eta()<-10 or tmp_H.Eta()<-10):
        return -100
    
    beam.SetPxPyPzE(0,0,6500,6500)
                    
    V_mom, bVH = TLorentzVector(), TVector3()
    V_mom = tmp_lep1+tmp_lep2
    bVH = (tmp_lep1+tmp_lep2+tmp_H).BoostVector()

    V_mom.Boost(-bVH)

    Theta = float('nan')

    try:
#   Theta  = acos((V_mom.Vect().Unit()).Dot(beam.Vect().Unit()))
        Theta = (V_mom.Vect().Unit()).Angle(beam.Vect().Unit())
    except Exception:
        pass
        Theta = -100

    return Theta;

def gettheta(lep1, lep2, H):

    tmp_lep1, tmp_lep2, tmp_H = TLorentzVector(), TLorentzVector(), TLorentzVector()

    tmp_lep1.SetPtEtaPhiM(lep1.Pt(),lep1.Eta(),lep1.Phi(),lep1.M())
    tmp_lep2.SetPtEtaPhiM(lep2.Pt(),lep2.Eta(),lep2.Phi(),lep2.M())
    tmp_H.SetPtEtaPhiM(H.Pt(),H.Eta(),H.Phi(),H.M())

    if(lep1.Eta()<-10 or lep2.Eta()<-10 or tmp_H.Eta()<-10):
        return -100
    
    V_mom, bVH, bV = TLorentzVector(), TVector3(), TVector3()

    bVH = (tmp_lep1 + tmp_lep2 + tmp_H).BoostVector()
    V_mom = (tmp_lep1 + tmp_lep2)

    V_mom.Boost(-bVH)
    tmp_lep1.Boost(-bVH)

    bV = V_mom.BoostVector()
    tmp_lep1.Boost(-bV)

    theta = float('nan')
    try:
        theta = (V_mom).Angle(tmp_lep1.Vect())
    except Exception:
        pass
        theta = -100

    return theta

def getphi(lep1, lep2, H):

    beam = TLorentzVector()

    tmp_lep1, tmp_lep2, tmp_H = TLorentzVector(), TLorentzVector(), TLorentzVector()

    tmp_lep1.SetPtEtaPhiM(lep1.Pt(),lep1.Eta(),lep1.Phi(),lep1.M())
    tmp_lep2.SetPtEtaPhiM(lep2.Pt(),lep2.Eta(),lep2.Phi(),lep2.M())
    tmp_H.SetPtEtaPhiM(H.Pt(),H.Eta(),H.Phi(),H.M())

    if(lep1.Eta()<-10 or lep2.Eta()<-10 or tmp_H.Eta()<-10):
        return -100

    beam.SetPxPyPzE(0,0,6500,6500)

    V_mom, bVH, n_scatter, n_decay = TLorentzVector(), TVector3(), TVector3(), TVector3()
    bVH = (tmp_lep1+tmp_lep2+tmp_H).BoostVector()
    V_mom = tmp_lep1+tmp_lep2

    tmp_lep1.Boost(-bVH)
    tmp_lep2.Boost(-bVH)
    V_mom.Boost(-bVH)

    n_scatter = ((beam.Vect().Unit()).Cross(V_mom.Vect())).Unit()
    n_decay   = (tmp_lep1.Vect().Cross(tmp_lep2.Vect())).Unit()

    sign_flip =  1 if ( ((n_scatter.Cross(n_decay))*(V_mom.Vect())) > 0 ) else -1

    try:
        phi = sign_flip*acos(n_scatter.Dot(n_decay))    
    except Exception:
        pass
        phi = -100

    return phi

def neutrino_mom(vec_lep, MET_pt, MET_phi, random):

    W_mass = 80.4

    pnu = ROOT.TLorentzVector()

    if vec_lep.E()<0:
        pnu.etPtEtaPhiM(0,-100,-100,0)
    else:
        mT2 = 2*vec_lep.Pt()*MET_pt*(1-cos(deltaPhi(vec_lep.Phi(),MET_phi)))
        Delta2 = (W_mass*W_mass - mT2)*1./(2*MET_pt*vec_lep.Pt())
        if (Delta2>=0):
            try:
                nueta = (vec_lep.Eta() + abs(acosh(1+(Delta2)))) if (random>=0.5) else (vec_lep.Eta() - abs(acosh(1+(Delta2))))
            except Exception:
                pass
                nueta = -100
            pnu.SetPtEtaPhiM(MET_pt,nueta,MET_phi,0)
        else:
            pnu.SetPtEtaPhiM(0,-100,-100,0)
    
    return pnu

#def neutrino_mom(vec_lep, MET_pt, MET_phi, random):
#
#    W_mass = 80.4
#
#    pnu = ROOT.TLorentzVector()
#
#    if vec_lep.E()<0:
#        pnu.etPtEtaPhiM(0,-100,-100,0)
#        raise RuntimeError("Negative lepton energy %3.2f. Should not happen." % vec_lep.E())
#    else:
#        mT2 = 2*vec_lep.Pt()*MET_pt*(1-cos(deltaPhi(vec_lep.Phi(),MET_phi)))
#        Delta2 = (W_mass*W_mass - mT2)*1./(2*MET_pt*vec_lep.Pt())
#        if (Delta2>=0):
#            try:
#                nueta = (vec_lep.Eta() + abs(acosh(1+(Delta2)))) if (random>=0.5) else (vec_lep.Eta() - abs(acosh(1+(Delta2))))
#            except Exception:
#                pass
#                nueta = -100
#            pnu.SetPtEtaPhiM(MET_pt,nueta,MET_phi,0)
#        else:
#            #http://cds.cern.ch/record/2757267/files/SMP-20-005-pas.pdf Sec. 6.1: If Delta<0 take eta_neu=eta_l 
#            pnu.SetPtEtaPhiM(MET_pt,vec_lep.Eta(),MET_phi,0)
#    
#    return pnu

if __name__=='__main__':

    # Let's construct 4-fectors in the V-H c.m.s

    phi_V = 0.1
    mV    = 91.2

    v1 = ROOT.TLorentzVector(mV/2.*cos(phi_V),mV/2*sin(phi_V), 0, mV/2. )
    v2 = ROOT.TLorentzVector(-mV/2.*cos(phi_V),-mV/2*sin(phi_V), 0, mV/2. )

    # Let's boost a little into the +z direction
    boost_V = ROOT.TVector3(0,0,0.1)
    v1.Boost(boost_V)
    v2.Boost(boost_V)
    
    phi_H = 0.2
    mH    = 125.

    v3 = ROOT.TLorentzVector(mH/2.*cos(phi_H),mH/2*sin(phi_H), 0, mH/2. )
    v4 = ROOT.TLorentzVector(-mH/2.*cos(phi_H),-mH/2*sin(phi_H), 0, mH/2. )

    # Let's boost a little into the -z direction
    boost_H = ROOT.TVector3(0,0,-0.1)
    v3.Boost(boost_H)
    v4.Boost(boost_H)

    # Finally, let's boost the VH system a little in the transverse plane 

    boost_VH = ROOT.TVector3(-0.1,0,0)
    for v in [ v1,v2,v3,v4 ]:
        v.Boost( boost_VH )

    #v1 = ROOT.TLorentzVector(100,10,70, sqrt(2*100**2+80**2))
    #v2 = ROOT.TLorentzVector(80,100,70, sqrt(2*100**2+80**2))
    #v3 = ROOT.TLorentzVector(-30,10,50, sqrt(2*100**2+80**2))
    #v4 = ROOT.TLorentzVector(90,100,-70, sqrt(2*100**2+80**2))

    #res = VV_angles( v1, v2, v3, v4, debug = True)

    print "suman Theta", getTheta(v1, v2, v3+v4), getTheta(v3, v4, v1+v2) #, "mine", res["Theta"]
    print "suman V1",    gettheta(v1, v2, v3+v4) #, "mine", res["theta_V1"]
    print "suman V2",    gettheta(v3, v4, v1+v2) #, "mine", res["theta_V2"]
    print "suman phi 1", getphi(v1, v2, v3+v4)   #, "mine", res["phi1"]
    print "suman phi 2", getphi(v3, v4, v1+v2)   #, "mine", res["phi2"]
