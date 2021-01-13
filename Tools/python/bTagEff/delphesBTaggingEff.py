def getBTagSF_1a(bTagWP, unc, bJets, nonBJets):

    if unc not in ['B', 'L']: return NotImplementedError
    if unc == 'B': 
        unc_func = unc_B
    elif unc == 'L':
        unc_func = unc_L
    
    ref = 1.
    var = 1.
    for j in bJets:
        flavor = abs(j['flavor']) if abs(j['flavor']) in [4,5] else 0
        ref *= bTagEff( bTagWP, flavor, j['pt'], j['eta'] )
        var *= bTagEff( bTagWP, flavor, j['pt'], j['eta'] )*(1. + unc_func(bTagWP, flavor, j['pt']))
    for j in nonBJets:
        flavor = abs(j['flavor']) if abs(j['flavor']) in [4,5] else 0
        ref *= 1-bTagEff( bTagWP, flavor, j['pt'], j['eta'] )
        var *= 1-bTagEff( bTagWP, flavor, j['pt'], j['eta'] )*(1. + unc_func(bTagWP, flavor, j['pt']))
    if ref>0:
        return var/ref

def unc_B( bTagWP, flavor, pt ):
    if bTagWP != 'medium':
        raise NotImplementedError( "This BTag WP is not implemented: %s", bTagWP )
    if abs(flavor)==5: 
        if pt>=20 and pt<30: return 0.017
        if pt>=30 and pt<50: return 0.01
        if pt>=50 and pt<70: return 0.01
        if pt>=70 and pt<100: return 0.01
        if pt>=100 and pt<140: return 0.01
        if pt>=140 and pt<200: return 0.01
        if pt>=200 and pt<300: return 0.016
        if pt>=300 and pt<600: return 0.018
        if pt>=600 and pt<1000: return 0.023
        if pt>=1000 : return 0.046
    elif abs(flavor)==4:
        if pt>=20 and pt<30: return 0.052
        if pt>=30 and pt<50: return 0.021
        if pt>=50 and pt<70: return 0.021
        if pt>=70 and pt<100: return 0.019
        if pt>=100 and pt<140: return 0.02
        if pt>=140 and pt<200: return 0.021
        if pt>=200 and pt<300: return 0.052
        if pt>=300 and pt<600: return 0.055
        if pt>=600 and pt<1000: return 0.069
        if pt>=1000 : return 0.139
    else:
        return 0

def unc_L( bTagWP, flavor, pt ):
    if bTagWP != 'medium':
        raise NotImplementedError( "This BTag WP is not implemented: %s", bTagWP )
    if flavor == 0:
        return 0.1
    else:
        return 0

def bTagEff( bTagWP, flavor, pt, eta ):
    if bTagWP != 'medium':
        raise NotImplementedError( "This BTag WP is not implemented: %s", bTagWP )
    if abs(flavor) == 0:
        return 0.01
    elif abs(flavor) == 4: return\
        (abs(eta) >= 0.00 and abs(eta) <= 1.50) * (pt > 20.00 and pt <= 30.00) * (0.126) +\
        (abs(eta) >= 0.00 and abs(eta) <= 1.50) * (pt > 30.00 and pt <= 40.00) * (0.129) +\
        (abs(eta) >= 0.00 and abs(eta) <= 1.50) * (pt > 40.00 and pt <= 50.00) * (0.131) +\
        (abs(eta) >= 0.00 and abs(eta) <= 1.50) * (pt > 50.00 and pt <= 60.00) * (0.133) +\
        (abs(eta) >= 0.00 and abs(eta) <= 1.50) * (pt > 60.00 and pt <= 70.00) * (0.139) +\
        (abs(eta) >= 0.00 and abs(eta) <= 1.50) * (pt > 70.00 and pt <= 80.00) * (0.143) +\
        (abs(eta) >= 0.00 and abs(eta) <= 1.50) * (pt > 80.00 and pt <= 90.00) * (0.141) +\
        (abs(eta) >= 0.00 and abs(eta) <= 1.50) * (pt > 90.00 and pt <= 100.00) * (0.151) +\
        (abs(eta) >= 0.00 and abs(eta) <= 1.50) * (pt > 100.00 and pt <= 120.00) * (0.151) +\
        (abs(eta) >= 0.00 and abs(eta) <= 1.50) * (pt > 120.00 and pt <= 140.00) * (0.157) +\
        (abs(eta) >= 0.00 and abs(eta) <= 1.50) * (pt > 140.00 and pt <= 160.00) * (0.167) +\
        (abs(eta) >= 0.00 and abs(eta) <= 1.50) * (pt > 160.00 and pt <= 180.00) * (0.170) +\
        (abs(eta) >= 0.00 and abs(eta) <= 1.50) * (pt > 180.00 and pt <= 200.00) * (0.175) +\
        (abs(eta) >= 0.00 and abs(eta) <= 1.50) * (pt > 200.00 and pt <= 250.00) * (0.174) +\
        (abs(eta) >= 0.00 and abs(eta) <= 1.50) * (pt > 250.00 and pt <= 300.00) * (0.172) +\
        (abs(eta) >= 0.00 and abs(eta) <= 1.50) * (pt > 300.00 and pt <= 350.00) * (0.173) +\
        (abs(eta) >= 0.00 and abs(eta) <= 1.50) * (pt > 350.00 and pt <= 400.00) * (0.163) +\
        (abs(eta) >= 0.00 and abs(eta) <= 1.50) * (pt > 400.00 and pt <= 500.00) * (0.158) +\
        (abs(eta) >= 0.00 and abs(eta) <= 1.50) * (pt > 500.00 and pt <= 600.00) * (0.138) +\
        (abs(eta) >= 0.00 and abs(eta) <= 1.50) * (pt > 600.00 and pt <= 700.00) * (0.127) +\
        (abs(eta) >= 0.00 and abs(eta) <= 1.50) * (pt > 700.00 and pt <= 800.00) * (0.127) +\
        (abs(eta) >= 0.00 and abs(eta) <= 1.50) * (pt > 800.00 and pt <= 1000.00) * (0.112) +\
        (abs(eta) >= 0.00 and abs(eta) <= 1.50) * (pt > 1000.00 and pt <= 1400.00) * (0.101) +\
        (abs(eta) >= 0.00 and abs(eta) <= 1.50) * (pt > 1400.00 and pt <= 2000.00) * (0.093) +\
        (abs(eta) >= 0.00 and abs(eta) <= 1.50) * (pt > 2000.00 and pt <= 3000.00) * (0.078) +\
        (abs(eta) >= 0.00 and abs(eta) <= 1.50) * (pt > 3000.00) * (0.078) +\
        (abs(eta) > 1.50 and abs(eta) <= 2.50) * (pt > 20.00 and pt <= 30.00) * (0.093) +\
        (abs(eta) > 1.50 and abs(eta) <= 2.50) * (pt > 30.00 and pt <= 40.00) * (0.100) +\
        (abs(eta) > 1.50 and abs(eta) <= 2.50) * (pt > 40.00 and pt <= 50.00) * (0.108) +\
        (abs(eta) > 1.50 and abs(eta) <= 2.50) * (pt > 50.00 and pt <= 60.00) * (0.116) +\
        (abs(eta) > 1.50 and abs(eta) <= 2.50) * (pt > 60.00 and pt <= 70.00) * (0.119) +\
        (abs(eta) > 1.50 and abs(eta) <= 2.50) * (pt > 70.00 and pt <= 80.00) * (0.122) +\
        (abs(eta) > 1.50 and abs(eta) <= 2.50) * (pt > 80.00 and pt <= 90.00) * (0.124) +\
        (abs(eta) > 1.50 and abs(eta) <= 2.50) * (pt > 90.00 and pt <= 100.00) * (0.132) +\
        (abs(eta) > 1.50 and abs(eta) <= 2.50) * (pt > 100.00 and pt <= 120.00) * (0.131) +\
        (abs(eta) > 1.50 and abs(eta) <= 2.50) * (pt > 120.00 and pt <= 140.00) * (0.132) +\
        (abs(eta) > 1.50 and abs(eta) <= 2.50) * (pt > 140.00 and pt <= 160.00) * (0.138) +\
        (abs(eta) > 1.50 and abs(eta) <= 2.50) * (pt > 160.00 and pt <= 180.00) * (0.130) +\
        (abs(eta) > 1.50 and abs(eta) <= 2.50) * (pt > 180.00 and pt <= 200.00) * (0.138) +\
        (abs(eta) > 1.50 and abs(eta) <= 2.50) * (pt > 200.00 and pt <= 250.00) * (0.132) +\
        (abs(eta) > 1.50 and abs(eta) <= 2.50) * (pt > 250.00 and pt <= 300.00) * (0.125) +\
        (abs(eta) > 1.50 and abs(eta) <= 2.50) * (pt > 300.00 and pt <= 350.00) * (0.116) +\
        (abs(eta) > 1.50 and abs(eta) <= 2.50) * (pt > 350.00 and pt <= 400.00) * (0.123) +\
        (abs(eta) > 1.50 and abs(eta) <= 2.50) * (pt > 400.00 and pt <= 500.00) * (0.113) +\
        (abs(eta) > 1.50 and abs(eta) <= 2.50) * (pt > 500.00 and pt <= 600.00) * (0.102) +\
        (abs(eta) > 1.50 and abs(eta) <= 2.50) * (pt > 600.00 and pt <= 700.00) * (0.091) +\
        (abs(eta) > 1.50 and abs(eta) <= 2.50) * (pt > 700.00 and pt <= 800.00) * (0.085) +\
        (abs(eta) > 1.50 and abs(eta) <= 2.50) * (pt > 800.00 and pt <= 1000.00) * (0.072) +\
        (abs(eta) > 1.50 and abs(eta) <= 2.50) * (pt > 1000.00 and pt <= 1400.00) * (0.072) +\
        (abs(eta) > 1.50 and abs(eta) <= 2.50) * (pt > 1400.00 and pt <= 2000.00) * (0.072) +\
        (abs(eta) > 1.50 and abs(eta) <= 2.50) * (pt > 2000.00 and pt <= 3000.00) * (0.072) +\
        (abs(eta) > 1.50 and abs(eta) <= 2.50) * (pt > 3000.00) * (0.072) +\
        (abs(eta) > 2.50 and abs(eta) <= 3.50) * (pt > 20.00 and pt <= 30.00) * (0.079) +\
        (abs(eta) > 2.50 and abs(eta) <= 3.50) * (pt > 30.00 and pt <= 40.00) * (0.086) +\
        (abs(eta) > 2.50 and abs(eta) <= 3.50) * (pt > 40.00 and pt <= 50.00) * (0.085) +\
        (abs(eta) > 2.50 and abs(eta) <= 3.50) * (pt > 50.00 and pt <= 60.00) * (0.090) +\
        (abs(eta) > 2.50 and abs(eta) <= 3.50) * (pt > 60.00 and pt <= 70.00) * (0.088) +\
        (abs(eta) > 2.50 and abs(eta) <= 3.50) * (pt > 70.00 and pt <= 80.00) * (0.088) +\
        (abs(eta) > 2.50 and abs(eta) <= 3.50) * (pt > 80.00 and pt <= 90.00) * (0.094) +\
        (abs(eta) > 2.50 and abs(eta) <= 3.50) * (pt > 90.00 and pt <= 100.00) * (0.091) +\
        (abs(eta) > 2.50 and abs(eta) <= 3.50) * (pt > 100.00 and pt <= 120.00) * (0.089) +\
        (abs(eta) > 2.50 and abs(eta) <= 3.50) * (pt > 120.00 and pt <= 140.00) * (0.089) +\
        (abs(eta) > 2.50 and abs(eta) <= 3.50) * (pt > 140.00 and pt <= 160.00) * (0.088) +\
        (abs(eta) > 2.50 and abs(eta) <= 3.50) * (pt > 160.00 and pt <= 180.00) * (0.085) +\
        (abs(eta) > 2.50 and abs(eta) <= 3.50) * (pt > 180.00 and pt <= 200.00) * (0.087) +\
        (abs(eta) > 2.50 and abs(eta) <= 3.50) * (pt > 200.00 and pt <= 250.00) * (0.098) +\
        (abs(eta) > 2.50 and abs(eta) <= 3.50) * (pt > 250.00 and pt <= 300.00) * (0.072) +\
        (abs(eta) > 2.50 and abs(eta) <= 3.50) * (pt > 300.00 and pt <= 350.00) * (0.098) +\
        (abs(eta) > 2.50 and abs(eta) <= 3.50) * (pt > 350.00 and pt <= 400.00) * (0.053) +\
        (abs(eta) > 2.50 and abs(eta) <= 3.50) * (pt > 400.00 and pt <= 500.00) * (0.081) +\
        (abs(eta) > 2.50 and abs(eta) <= 3.50) * (pt > 500.00 and pt <= 600.00) * (0.081) +\
        (abs(eta) > 2.50 and abs(eta) <= 3.50) * (pt > 600.00 and pt <= 700.00) * (0.081) +\
        (abs(eta) > 2.50 and abs(eta) <= 3.50) * (pt > 700.00 and pt <= 800.00) * (0.081) +\
        (abs(eta) > 2.50 and abs(eta) <= 3.50) * (pt > 800.00 and pt <= 1000.00) * (0.081) +\
        (abs(eta) > 2.50 and abs(eta) <= 3.50) * (pt > 1000.00 and pt <= 1400.00) * (0.081) +\
        (abs(eta) > 2.50 and abs(eta) <= 3.50) * (pt > 1400.00 and pt <= 2000.00) * (0.081) +\
        (abs(eta) > 2.50 and abs(eta) <= 3.50) * (pt > 2000.00 and pt <= 3000.00) * (0.081) +\
        (abs(eta) > 2.50 and abs(eta) <= 3.50) * (pt > 3000.00) * (0.081)
    elif abs(flavor) == 5: return\
        (abs(eta) >= 0.00 and abs(eta) <= 1.50) * (pt > 20.00 and pt <= 30.00) * (0.663) +\
        (abs(eta) >= 0.00 and abs(eta) <= 1.50) * (pt > 30.00 and pt <= 40.00) * (0.715) +\
        (abs(eta) >= 0.00 and abs(eta) <= 1.50) * (pt > 40.00 and pt <= 50.00) * (0.740) +\
        (abs(eta) >= 0.00 and abs(eta) <= 1.50) * (pt > 50.00 and pt <= 60.00) * (0.752) +\
        (abs(eta) >= 0.00 and abs(eta) <= 1.50) * (pt > 60.00 and pt <= 70.00) * (0.761) +\
        (abs(eta) >= 0.00 and abs(eta) <= 1.50) * (pt > 70.00 and pt <= 80.00) * (0.765) +\
        (abs(eta) >= 0.00 and abs(eta) <= 1.50) * (pt > 80.00 and pt <= 90.00) * (0.765) +\
        (abs(eta) >= 0.00 and abs(eta) <= 1.50) * (pt > 90.00 and pt <= 100.00) * (0.770) +\
        (abs(eta) >= 0.00 and abs(eta) <= 1.50) * (pt > 100.00 and pt <= 120.00) * (0.764) +\
        (abs(eta) >= 0.00 and abs(eta) <= 1.50) * (pt > 120.00 and pt <= 140.00) * (0.761) +\
        (abs(eta) >= 0.00 and abs(eta) <= 1.50) * (pt > 140.00 and pt <= 160.00) * (0.754) +\
        (abs(eta) >= 0.00 and abs(eta) <= 1.50) * (pt > 160.00 and pt <= 180.00) * (0.749) +\
        (abs(eta) >= 0.00 and abs(eta) <= 1.50) * (pt > 180.00 and pt <= 200.00) * (0.737) +\
        (abs(eta) >= 0.00 and abs(eta) <= 1.50) * (pt > 200.00 and pt <= 250.00) * (0.720) +\
        (abs(eta) >= 0.00 and abs(eta) <= 1.50) * (pt > 250.00 and pt <= 300.00) * (0.693) +\
        (abs(eta) >= 0.00 and abs(eta) <= 1.50) * (pt > 300.00 and pt <= 350.00) * (0.674) +\
        (abs(eta) >= 0.00 and abs(eta) <= 1.50) * (pt > 350.00 and pt <= 400.00) * (0.643) +\
        (abs(eta) >= 0.00 and abs(eta) <= 1.50) * (pt > 400.00 and pt <= 500.00) * (0.615) +\
        (abs(eta) >= 0.00 and abs(eta) <= 1.50) * (pt > 500.00 and pt <= 600.00) * (0.581) +\
        (abs(eta) >= 0.00 and abs(eta) <= 1.50) * (pt > 600.00 and pt <= 700.00) * (0.539) +\
        (abs(eta) >= 0.00 and abs(eta) <= 1.50) * (pt > 700.00 and pt <= 800.00) * (0.521) +\
        (abs(eta) >= 0.00 and abs(eta) <= 1.50) * (pt > 800.00 and pt <= 1000.00) * (0.468) +\
        (abs(eta) >= 0.00 and abs(eta) <= 1.50) * (pt > 1000.00 and pt <= 1400.00) * (0.433) +\
        (abs(eta) >= 0.00 and abs(eta) <= 1.50) * (pt > 1400.00 and pt <= 2000.00) * (0.387) +\
        (abs(eta) >= 0.00 and abs(eta) <= 1.50) * (pt > 2000.00 and pt <= 3000.00) * (0.345) +\
        (abs(eta) >= 0.00 and abs(eta) <= 1.50) * (pt > 3000.00) * (0.345) +\
        (abs(eta) > 1.50 and abs(eta) <= 2.50) * (pt > 20.00 and pt <= 30.00) * (0.520) +\
        (abs(eta) > 1.50 and abs(eta) <= 2.50) * (pt > 30.00 and pt <= 40.00) * (0.598) +\
        (abs(eta) > 1.50 and abs(eta) <= 2.50) * (pt > 40.00 and pt <= 50.00) * (0.642) +\
        (abs(eta) > 1.50 and abs(eta) <= 2.50) * (pt > 50.00 and pt <= 60.00) * (0.663) +\
        (abs(eta) > 1.50 and abs(eta) <= 2.50) * (pt > 60.00 and pt <= 70.00) * (0.672) +\
        (abs(eta) > 1.50 and abs(eta) <= 2.50) * (pt > 70.00 and pt <= 80.00) * (0.676) +\
        (abs(eta) > 1.50 and abs(eta) <= 2.50) * (pt > 80.00 and pt <= 90.00) * (0.679) +\
        (abs(eta) > 1.50 and abs(eta) <= 2.50) * (pt > 90.00 and pt <= 100.00) * (0.680) +\
        (abs(eta) > 1.50 and abs(eta) <= 2.50) * (pt > 100.00 and pt <= 120.00) * (0.674) +\
        (abs(eta) > 1.50 and abs(eta) <= 2.50) * (pt > 120.00 and pt <= 140.00) * (0.661) +\
        (abs(eta) > 1.50 and abs(eta) <= 2.50) * (pt > 140.00 and pt <= 160.00) * (0.654) +\
        (abs(eta) > 1.50 and abs(eta) <= 2.50) * (pt > 160.00 and pt <= 180.00) * (0.629) +\
        (abs(eta) > 1.50 and abs(eta) <= 2.50) * (pt > 180.00 and pt <= 200.00) * (0.614) +\
        (abs(eta) > 1.50 and abs(eta) <= 2.50) * (pt > 200.00 and pt <= 250.00) * (0.593) +\
        (abs(eta) > 1.50 and abs(eta) <= 2.50) * (pt > 250.00 and pt <= 300.00) * (0.558) +\
        (abs(eta) > 1.50 and abs(eta) <= 2.50) * (pt > 300.00 and pt <= 350.00) * (0.537) +\
        (abs(eta) > 1.50 and abs(eta) <= 2.50) * (pt > 350.00 and pt <= 400.00) * (0.509) +\
        (abs(eta) > 1.50 and abs(eta) <= 2.50) * (pt > 400.00 and pt <= 500.00) * (0.498) +\
        (abs(eta) > 1.50 and abs(eta) <= 2.50) * (pt > 500.00 and pt <= 600.00) * (0.467) +\
        (abs(eta) > 1.50 and abs(eta) <= 2.50) * (pt > 600.00 and pt <= 700.00) * (0.443) +\
        (abs(eta) > 1.50 and abs(eta) <= 2.50) * (pt > 700.00 and pt <= 800.00) * (0.392) +\
        (abs(eta) > 1.50 and abs(eta) <= 2.50) * (pt > 800.00 and pt <= 1000.00) * (0.361) +\
        (abs(eta) > 1.50 and abs(eta) <= 2.50) * (pt > 1000.00 and pt <= 1400.00) * (0.361) +\
        (abs(eta) > 1.50 and abs(eta) <= 2.50) * (pt > 1400.00 and pt <= 2000.00) * (0.361) +\
        (abs(eta) > 1.50 and abs(eta) <= 2.50) * (pt > 2000.00 and pt <= 3000.00) * (0.361) +\
        (abs(eta) > 1.50 and abs(eta) <= 2.50) * (pt > 3000.00) * (0.361) +\
        (abs(eta) > 2.50 and abs(eta) <= 3.50) * (pt > 20.00 and pt <= 30.00) * (0.357) +\
        (abs(eta) > 2.50 and abs(eta) <= 3.50) * (pt > 30.00 and pt <= 40.00) * (0.442) +\
        (abs(eta) > 2.50 and abs(eta) <= 3.50) * (pt > 40.00 and pt <= 50.00) * (0.481) +\
        (abs(eta) > 2.50 and abs(eta) <= 3.50) * (pt > 50.00 and pt <= 60.00) * (0.496) +\
        (abs(eta) > 2.50 and abs(eta) <= 3.50) * (pt > 60.00 and pt <= 70.00) * (0.505) +\
        (abs(eta) > 2.50 and abs(eta) <= 3.50) * (pt > 70.00 and pt <= 80.00) * (0.504) +\
        (abs(eta) > 2.50 and abs(eta) <= 3.50) * (pt > 80.00 and pt <= 90.00) * (0.506) +\
        (abs(eta) > 2.50 and abs(eta) <= 3.50) * (pt > 90.00 and pt <= 100.00) * (0.502) +\
        (abs(eta) > 2.50 and abs(eta) <= 3.50) * (pt > 100.00 and pt <= 120.00) * (0.492) +\
        (abs(eta) > 2.50 and abs(eta) <= 3.50) * (pt > 120.00 and pt <= 140.00) * (0.487) +\
        (abs(eta) > 2.50 and abs(eta) <= 3.50) * (pt > 140.00 and pt <= 160.00) * (0.463) +\
        (abs(eta) > 2.50 and abs(eta) <= 3.50) * (pt > 160.00 and pt <= 180.00) * (0.470) +\
        (abs(eta) > 2.50 and abs(eta) <= 3.50) * (pt > 180.00 and pt <= 200.00) * (0.446) +\
        (abs(eta) > 2.50 and abs(eta) <= 3.50) * (pt > 200.00 and pt <= 250.00) * (0.447) +\
        (abs(eta) > 2.50 and abs(eta) <= 3.50) * (pt > 250.00 and pt <= 300.00) * (0.393) +\
        (abs(eta) > 2.50 and abs(eta) <= 3.50) * (pt > 300.00 and pt <= 350.00) * (0.390) +\
        (abs(eta) > 2.50 and abs(eta) <= 3.50) * (pt > 350.00 and pt <= 400.00) * (0.350) +\
        (abs(eta) > 2.50 and abs(eta) <= 3.50) * (pt > 400.00 and pt <= 500.00) * (0.398) +\
        (abs(eta) > 2.50 and abs(eta) <= 3.50) * (pt > 500.00 and pt <= 600.00) * (0.398) +\
        (abs(eta) > 2.50 and abs(eta) <= 3.50) * (pt > 600.00 and pt <= 700.00) * (0.398) +\
        (abs(eta) > 2.50 and abs(eta) <= 3.50) * (pt > 700.00 and pt <= 800.00) * (0.398) +\
        (abs(eta) > 2.50 and abs(eta) <= 3.50) * (pt > 800.00 and pt <= 1000.00) * (0.398) +\
        (abs(eta) > 2.50 and abs(eta) <= 3.50) * (pt > 1000.00 and pt <= 1400.00) * (0.398) +\
        (abs(eta) > 2.50 and abs(eta) <= 3.50) * (pt > 1400.00 and pt <= 2000.00) * (0.398) +\
        (abs(eta) > 2.50 and abs(eta) <= 3.50) * (pt > 2000.00 and pt <= 3000.00) * (0.398) +\
        (abs(eta) > 2.50 and abs(eta) <= 3.50) * (pt > 3000.00) * (0.398)
