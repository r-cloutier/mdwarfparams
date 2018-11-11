# use the results from the Kepler-GAIA cross-match here: https://gaia-kepler.fun
from imports import *
from query_gaia_2MASS import *


def get_initial_KepID_data(fout):
    '''Get the available data for the Kepler-GAIA cross-matched stars.'''
    # get cross-matach results with different search radii
    fs = np.sort(glob.glob('input_data/Keplertargets/kepler_dr2_*fits'))

    # get data from the cross-match
    KepIDs = np.zeros(0)
    ras = np.zeros(0)
    decs = np.zeros(0)
    ls = np.zeros(0)
    bs = np.zeros(0)
    GBPmags = np.zeros(0)
    GRPmags = np.zeros(0) 
    e_GBPmags = np.zeros(0)
    e_GRPmags = np.zeros(0) 
    Kepmags = np.zeros(0)
    pars = np.zeros(0)
    e_pars = np.zeros(0)
    Teff_tmp = np.zeros(0)
    Jmags = np.zeros(0)
    Hmags = np.zeros(0)
    Kmags = np.zeros(0)
    for i in range(fs.size):
        
        f = fits.open(fs[i])[1].data
        KepIDs = np.append(KepIDs, f['kepid'])
        ras = np.append(ras, f['ra'])
        decs = np.append(decs, f['dec'])
        ls = np.append(ls, f['l'])
        bs = np.append(bs, f['b'])

        GBPmags = np.append(GBPmags, f['phot_bp_mean_mag'])
        FBP = f['phot_bp_mean_flux']
        eFBP = f['phot_bp_mean_flux_error']
        e_GBPmags = np.append(e_GBPmags, -2.5*np.log10(FBP / (FBP+eFBP)))

        GRPmags = np.append(GRPmags, f['phot_rp_mean_mag'])
        FRP = f['phot_rp_mean_flux']
        eFRP = f['phot_rp_mean_flux_error']
        e_GRPmags = np.append(e_GRPmags, -2.5*np.log10(FRP / (FRP+eFRP)))

        Kepmags = np.append(Kepmags, f['kepmag'])
        pars = np.append(pars, f['parallax']) + .029  # systemic correction
        e_pars = np.append(e_pars, f['parallax_error'])
        Teff_tmp = np.append(Teff_tmp, f['teff'])
        Jmags = np.append(Jmags, f['jmag'])
        Hmags = np.append(Hmags, f['hmag'])
        Kmags = np.append(Kmags, f['kmag'])

    # make intitial list of M dwarfs
    g = Teff_tmp <= 4200
    KepIDs, ras, decs, ls, bs = KepIDs[g], ras[g], decs[g], ls[g], bs[g]
    GBPmags, e_GBPmags, GRPmags, e_GRPmags = GBPmags[g], e_GBPmags[g], GRPmags[g], e_GRPmags[g]
    Kepmags, pars, e_pars, Teff_tmp = Kepmags[g], pars[g], e_pars[g], Teff_tmp[g]
    Jmags, Hmags, Kmags = Jmags[g], Hmags[g], Kmags[g]

    # get 2MASS uncertainties which are not included in the cross-match data
    e_Jmags, e_Hmags, e_Kmags = get_2MASS(ras, decs, Jmags, Hmags, Kmags)
    hdr = 'KepID,ra_deg,dec_deg,Kepmag,parallax_mas,e_parallax,Jmag,'+ \
          'e_Jmag,Hmag,e_Hmag,Kmag,e_Kmag'
    outarr_tmp = np.array([KepIDs, ras, decs, Kepmags, pars, e_pars, Jmags,
                           e_Jmags, Hmags, e_Hmags, Kmags, e_Kmags]).T
    np.savetxt(fout, outarr_tmp, delimiter=',', header=hdr, fmt='%.8e')
    
    
    # compute quantities
    dists, mus = compute_distance_modulus(unp.uarray(pars,e_pars))
    AKs = compute_AK_mwdust(ls, bs, unp.nominal_values(dists),
                            unp.std_devs(dists))
    MKs = compute_MK(unp.uarray(Kmags, e_Kmags), mus, AKs)
    Rss = MK2Rs(MKs)
    Teffs = gaia2Teff(unp.uarray(GBPmags, e_GBPmags),
                      unp.uarray(GRPmags, e_GRPmags),
                      unp.uarray(Jmags, e_Jmags),
                      unp.uarray(Hmags, e_Hmags))
    Mss = MK2Ms(MK)
    loggs = compute_logg(Mss, Rss)

    # save results
    hdr = 'KepID,ra_deg,dec_deg,Kepmag,parallax_mas,e_parallax,Jmag,'+ \
          'e_Jmag,Hmag,e_Hmag,Kmag,e_Kmag,dist_pc,e_dist,mu,e_mu,AK,e_AK,'+ \
          'MK,e_MK,Rs_RSun,e_Rs,Teff_K,e_Teff,Ms_MSun,e_Ms,logg_dex,e_logg'
    outarr = np.array([KepIDs, ras, decs, Kepmags, pars, e_pars, Jmags,
                       e_Jmags, Hmags, e_Hmags, Kmags, e_Kmags,
                       unp.nominal_values(dists), unp.std_devs(dists),
                       unp.nominal_values(mus), unp.std_devs(mus),
                       unp.nominal_values(AKs), unp.std_devs(AKs),
                       unp.nominal_values(MKs), unp.std_devs(MKs),
                       unp.nominal_values(Rss), unp.std_devs(Rss),
                       unp.nominal_values(Teffs), unp.std_devs(Teffs),
                       unp.nominal_values(Mss), unp.std_devs(Mss),
                       unp.nominal_values(loggs), unp.std_devs(loggs)]).T
    np.savetxt(fout, outarr, delimiter=',', header=hdr, fmt='%.8e')
    
        

def get_2MASS(ras_deg, decs_deg, Jmags, Hmags, Kmags,
              radius_deg=.017, phot_rtol=.02):
    '''Match Kepler stars with GAIA data to the 2MASS point-source catlog to 
    retrieve photometric uncertainties.'''
    # get 2MASS data
    d = np.load('input_data/Keplertargets/fp_2mass.fp_psc12298.npy')
    inds = np.array([0,1,3,5,6,8,9,11])
    ras2M, decs2M, J2M, eJ2M, H2M, eH2M, K2M, eK2M = d[:,inds].T
    
    # match each star individually
    Nstars = ras_deg.size
    e_Jmags, e_Hmags, e_Kmags = np.zeros(Nstars), np.zeros(Nstars), \
                                np.zeros(Nstars)
    for i in range(Nstars):

        if i % 1e3 == 0:
            print float(i) / Nstars
        
        # get matching photometry between Kepler-GAIA and 2MASS
        g = (ras2M >= ras_deg[i] - radius_deg) & \
            (ras2M <= ras_deg[i] + radius_deg) & \
            (decs2M >= decs_deg[i] - radius_deg) & \
            (decs2M <= decs_deg[i] + radius_deg) & \
            np.isclose(J2M, Jmags[i], rtol=phot_rtol) & \
            np.isclose(H2M, Hmags[i], rtol=phot_rtol) & \
            np.isclose(K2M, Kmags[i], rtol=phot_rtol)

        if g.sum() > 0:
            g2 = abs(K2M[g]-Kmags[i]) == np.min(abs(K2M[g]-Kmags[i]))
            e_Jmags[i] = eJ2M[g][g2]
            e_Hmags[i] = eH2M[g][g2]
            e_Kmags[i] = eK2M[g][g2]

        else:
            e_Jmags[i], e_Hmags[i], e_Kmags[i] = np.repeat(np.nan, 3)

    return e_Jmags, e_Hmags, e_Kmags 



def gaia2Teff(GBPmag, GRPmag, Jmag, Hmag):
    '''Use the relation from Mann+2015 (table 2)
    http://adsabs.harvard.edu/abs/2015ApJ...804...64M
    '''
    a, b, c, d, e, f, g = 3.172, -2.475, 1.082, -.2231, .01738, .08776, -.04355
    pG = np.poly1d((e,d,c,b,a))
    p2 = np.poly1d((g,f,0))
    Teff = 35e2 * (pG(GBPmag-GRPmag) + p2(Jmag-Hmag))
    eTeff = np.sqrt(unp.std_devs(Teff)**2 + 49**2)
    return unp.uarray(unp.nominal_values(Teff), eTeff)


if __name__ == '__main__':
    fout = 'input_data/Keplertargets/KepMdwarfsv5.csv'
    get_initial_KepID_data(fout)
