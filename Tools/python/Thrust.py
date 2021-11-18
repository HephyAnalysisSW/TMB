''' Suman's Thrust implementation from VH
    https://github.com/HephyAnalysisSW/VH/blob/main/Histogram_producer/produce_histos_v18.py#L1245-L1256
'''

from ROOT import TVector3

def SIGN_FLIP(value):
    if value>=0:
        return +1.
    else:
        return -1.

def Thrust_axis(vecs):

    axes = []
    niter = 2**(min(3,len(vecs))-1)
    for iter in range(niter):
        t_axis = TVector3()
        if iter==0:
            t_axis = vecs[0]+vecs[1]+vecs[2]
        if iter==1:
            t_axis = -vecs[0]+vecs[1]+vecs[2]
        if iter==2:
            t_axis = vecs[0]-vecs[1]+vecs[2]
        if iter==3:
            t_axis = vecs[0]+vecs[1]-vecs[2]
        axes.append(t_axis.Unit())

    final_axis = int(-1)
    maxthrust = float(-1000)
    for iax in range(len(axes)):
        thrust = float(0)
        for v in vecs:
            thrust += abs(v.Dot(axes[iax]))
        if (thrust>maxthrust):
            maxthrust = thrust
            final_axis = iax

    if final_axis>=0:
        return axes[final_axis]
    else:
        return vecs[0].Unit()

def Thrust(V,js):

    vecs = [j.Vect() for j in js]
    vecs.append(V.Vect())

    for v in vecs:
        v.SetZ(0)

    axis = TVector3()
    if len(vecs)<2: 
        return (0,0)
    else:
        if len(vecs)==2:
            axis = vecs[0].Unit()
        if len(vecs)>2:
            axis = Thrust_axis(vecs)
        beam, m_axis = TVector3(),TVector3()
        beam.SetXYZ(0,0,1)
        m_axis = axis.Cross(beam)
        thrust, thrust_min, momsum = float(0), float(0), float(0)
        for v in vecs:
            thrust += abs(v.Dot(axis))
            thrust_min += abs(v.Dot(m_axis))
            momsum += abs(v.Mag())
        return (thrust*1./max(1.e-6,momsum), thrust_min*1./max(1.e-6,momsum))
