from occurrencerateclass import *


def estimate_transmission_metric(Ps, rps, Jmags, Mss, Rss, Teffs):
    '''from Kempton+2018 Eq 1'''
    # define S/N scale factors
    scale_factors = {'terrestrial':.19, 'super-Earth':1.26, 'sub-Neptune':1.28,
                     'giant':1.15}

    # make arrays
    Ps, rps, Jmags, Mss, Rss, Teffs = np.array([Ps]).flatten(), \
                                      np.array([rps]).flatten(), \
                                      np.array([Jmags]).flatten(), \
                                      np.array([Mss]).flatten(), \
                                      np.array([Rss]).flatten(), \
                                      np.array([Teffs]).flatten()
    
    # compute masses
    mps = rad2mass(rps)
    Ks = rvs.RV_K(Ps, Mss, mps)
    
    # compute Teq assuming zero albedo
    smas = rvs.semimajoraxis(Ps, Mss, mps)
    Teqs = Teffs * np.sqrt(rvs.Rsun2m(Rss) / (2*rvs.AU2m(smas)))

    # get scale factors to scale analytic S/N to the simulations from Louie+2018
    SFs = np.zeros(Ps.size)
    for i in range(Ps.size):
        if rps[i] < 1.5:
            SFs[i] = scale_factors['terrestrial']
        elif 1.5 < rps[i] < 2.75:
            SFs[i] = scale_factors['super-Earth']
        elif 2.75 < rps[i] < 4:
            SFs[i] = scale_factors['sub-Neptune']
        elif 4 < rps[i] < 11:
            SFs[i] = scale_factors['giant']
        else:
            raise ValueError('%.3f Rearth is not a valid rp'%rps[i])
    
    # compute the transmission spectroscopy metric which is propto S/N
    TMSs = SFs * rps**3 * Teqs * 10**(-.2*Jmags) / (mps * Rss**2)
    return mps, Ks, Teqs, SFs, TMSs



def estimate_eclipse_metric(Ps, rps, Kmags, Mss, Rss, Teffs):
    '''from Kempton+2018 Eq 4'''
    # make arrays
    Ps, rps, Kmags, Mss, Rss, Teffs = np.array([Ps]).flatten(), \
                                      np.array([rps]).flatten(), \
                                      np.array([Kmags]).flatten(), \
                                      np.array([Mss]).flatten(), \
                                      np.array([Rss]).flatten(), \
                                      np.array([Teffs]).flatten()

    # compute mass
    mps = rad2mass(rps)
    Ks = rvs.RV_K(Ps, Mss, mps)
    
    # compute Teq assuming zero albedo
    smas = rvs.semimajoraxis(Ps, Mss, mps)
    Teqs = Teffs * np.sqrt(rvs.Rsun2m(Rss) / (2*rvs.AU2m(smas)))
    
    # evaluate the Planck function for the star and planet-day-side (ie 1.1*Teq)
    B7d5p = Planck(1.1*Teqs, 7.5e-6)
    B7d5s = Planck(Teffs, 7.5e-6)
    Bratio = (B7d5p/B7d5s)
    
    # compute the eclipse spectroscopy metric which is propto S/N
    EMSs = 4.29e6 * Bratio * (rvs.Rearth2m(rps)/rvs.Rsun2m(Rss))**2 * \
           10**(-.2*Kmags)
    return mps, Ks, Teqs, Bratio, EMSs



def rad2mass(rps, maxE=1.23, maxN=14.26):
    '''Convert the input radii in Earth radii to masses using the
    analytical model from Chen & Kipping 2017'''
    mps = np.zeros(rps.size)
    for i in range(rps.size):
        if rps[i] < maxE:
            mps[i] = .9781*rps[i]**(3.58)
        elif rps[i] < maxN:
            mps[i] = 1.436 * rps[i]**(1.7)
    return mps


def Planck(T_K, wl_meters):
    h = 6.626070e-34
    c = 2.997e8
    kb = 1.380648e-23
    return 2.*h*c / (wl_meters**5) * 1./(np.exp(h*c/(wl_meters*kb*T_K)) - 1)
