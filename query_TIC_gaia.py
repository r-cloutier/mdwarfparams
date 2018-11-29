# get GAIA DR2 parallaxes and photometry for stars in the input list of
# TIC M dwarfs
from imports import *
import astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.gaia import Gaia
from astroquery.vizier import Vizier
import rvs
from send_email import *
from uncertainties import unumpy as unp
from planetsearch import *
import mwdust


def get_stellar_data_TIC(tics, fout, radius_arcsec=10,
                         overwrite=False):

    # get stellar data from the full TIC
    ticfull_fname = 'input_data/TESStargets/TICv7_Mdwarfsv1.csv'
    cols = (0,13,14,42,43,44,45,46,47,58,60,84,87)
    d = np.loadtxt(ticfull_fname, delimiter=',', skiprows=5, usecols=cols)
    ticsF, rasF, decsF, JmagsF, e_JmagsF, HmagsF, e_HmagsF, KmagsF, e_KmagsF, GAIAmagsF, TESSmagsF, contratioF, priorityF = d.T

    Nstars = tics.size
    radius_arcsec_orig = radius_arcsec+0
    ras, decs = np.zeros(Nstars)+np.nan, np.zeros(Nstars)+np.nan
    GAIAmags, TESSmags = np.zeros(Nstars), np.zeros(Nstars)
    Jmags, e_Jmags = np.zeros(Nstars), np.zeros(Nstars)
    Hmags, e_Hmags = np.zeros(Nstars), np.zeros(Nstars)
    Kmags, e_Kmags = np.zeros(Nstars), np.zeros(Nstars)
    GBPmags, GRPmags = np.zeros((Nstars,2)), np.zeros((Nstars,2))
    pars = np.zeros((Nstars,2))
    dists, mus = np.zeros((Nstars,2)), np.zeros((Nstars,2))
    AKs, MKs = np.zeros((Nstars,2)), np.zeros((Nstars,2))
    Teffs, Mss = np.zeros((Nstars,2)), np.zeros((Nstars,2))
    Rss, loggs = np.zeros((Nstars, 2)), np.zeros((Nstars, 2))
    for i in range(Nstars):

        print float(i) / Nstars        
        # search gaia and 2MASS until we find a likely match
        # based on photometry
        if np.in1d(tics[i],ticsF):

            g = np.where(np.in1d(ticsF,tics[i]))[0]
            ras[i], decs[i] = rasF[g], decsF[g]
            TESSmags[i], GAIAmags[i] = TESSmagsF[g], GAIAmagsF[g]
            Jmags[i], e_Jmags[i] = JmagsF[g], e_JmagsF[g]
            Hmags[i], e_Hmags[i] = HmagsF[g], e_HmagsF[g]
            Kmags[i], e_Kmags[i] = KmagsF[g], e_KmagsF[g]
            
            pars[i] = np.nan, np.nan
            radius_arcsec = radius_arcsec_orig
            while (np.isnan(pars[i,0])) and (radius_arcsec <= 60):
                theta = ras[i], decs[i], GAIAmags[i], Jmags[i], e_Jmags[i], \
                        Hmags[i], e_Hmags[i], Kmags[i], e_Kmags[i]
                p = query_one_TIC(theta, radius_arcsec=radius_arcsec)
                pars[i],GBPmags[i],GRPmags[i],dists[i],mus[i],AKs[i],MKs[i],Rss[i],Teffs[i],Mss[i]=p
                logg = compute_logg(unp.uarray(Mss[i,0],Mss[i,1]),
                                    unp.uarray(Rss[i,0],Rss[i,1]))
                loggs[i] = unp.nominal_values(logg), unp.std_devs(logg)
                radius_arcsec += 5

        else:
            p = np.repeat(np.nan,20).reshape(10,2)
            pars[i],GBPmags[i],GRPmags[i],dists[i],mus[i],AKs[i],MKs[i],Rss[i],Teffs[i],Mss[i]=p
            logg = compute_logg(unp.uarray(Mss[i,0],Mss[i,1]),
                                unp.uarray(Rss[i,0],Rss[i,1]))
            loggs[i] = unp.nominal_values(logg), unp.std_devs(logg)
            
    # save results
    hdr = 'TIC,ra_deg,dec_deg,GBPmag,e_GBPmag,GRPmag,e_GRPmag,TESSmag,Jmag,'+ \
          'e_Jmag,Hmag,e_Hmag,Kmag,e_Kmag,parallax_mas,e_parallax,dist_pc,'+ \
          'ehi_dist,elo_dist,mu,ehi_mu,elo_mu,AK,e_AK,MK,ehi_MK,elo_MK,'+ \
          'Rs_RSun,ehi_Rs,elo_Rs,Teff_K,ehi_Teff,elo_Teff,Ms_MSun,ehi_Ms,'+ \
          'elo_Ms,logg_dex,ehi_logg,elo_logg'
    g = np.isfinite(ras)
    outarr = np.array([tics[g], ras[g], decs[g], GBPmags[g,0], GBPmags[g,1],
                       GRPmags[g,0], GBPmags[g,1], TESSmags[g], Jmags[g],
                       e_Jmags[g], Hmags[g], e_Hmags[g], Kmags[g], e_Kmags[g],
                       pars[g,0], pars[g,1], dists[g,0], dists[g,1], dists[g,1],
                       mus[g,0], mus[g,1], mus[g,1], AKs[g,0], AKs[g,1],
                       MKs[g,0], MKs[g,1], MKs[g,1], Rss[g,0], Rss[g,1],
                       Rss[g,1], Teffs[g,0], Teffs[g,1], Teffs[g,1], Mss[g,0],
                       Mss[g,1], Mss[g,1], loggs[g,0], loggs[g,1],loggs[g,1]]).T

    if os.path.exists(fout) and not overwrite:
        inarr = np.loadtxt(fout, delimiter=',')
        assert inarr.shape[1] == outarr.shape[1]
        np.savetxt(fout, np.append(inarr, outarr, axis=0),
                   delimiter=',', header=hdr, fmt='%.8e')
    else:
        np.savetxt(fout, outarr, delimiter=',', header=hdr, fmt='%.8e')




def query_one_TIC(theta, radius_arcsec=10):
    
    # get gaia data
    assert len(theta) == 9
    ra_deg,dec_deg,GAIAmag,Jmag,e_Jmag,Hmag,e_Hmag,Kmag,e_Kmag = theta
    coord = SkyCoord(ra=ra_deg, dec=dec_deg,
                     unit=(u.degree, u.degree), frame='icrs')
    rad = u.Quantity(float(radius_arcsec), u.arcsec)
    rG = Gaia.query_object_async(coordinate=coord, radius=rad)
    # return NaN if no GAIA star is found
    if len(rG) == 0:
        return np.repeat(np.nan, 20).reshape(10,2)

    # step through possible matches
    for i in range(len(rG)):

        # get gaia parameters
        par, epar = rG['parallax'][i]+.029, rG['parallax_error'][i]
        GAIAmag_dr2 = rG['phot_g_mean_mag'][i]
        GBPmag = rG['phot_bp_mean_mag'][i]
        FBP = rG['phot_bp_mean_flux'][i]
        eFBP = rG['phot_bp_mean_flux_error'][i]
        eGBPmag = -2.5*np.log10(FBP / (FBP+eFBP))
        GRPmag = rG['phot_rp_mean_mag'][i]
        FRP = rG['phot_rp_mean_flux'][i]
        eFRP = rG['phot_rp_mean_flux_error'][i]
        eGRPmag = -2.5*np.log10(FRP / (FRP+eFRP))

        # change bad values to NaN
        if np.ma.is_masked(par) or par<=0: par, epar = np.nan, np.nan 
        if np.ma.is_masked(epar) or epar<=0: par, epar = np.nan, np.nan 
        if np.ma.is_masked(GBPmag): GBPmag, eGBPmag = np.nan, np.nan
        if np.ma.is_masked(eGBPmag): GBPmag, eGBPmag = np.nan, np.nan
        if np.ma.is_masked(GRPmag): GRPmag, eGRPmag = np.nan, np.nan
        if np.ma.is_masked(eGRPmag): GRPmag, eGRPmag = np.nan, np.nan
            
        # check that GAIA and 2MASS photometry aproximately match
        if np.isnan(par) or np.isnan(GBPmag) or np.isnan(GRPmag):
            match = False
        elif GAIAmag == GAIAmag_dr2:
            match = True
        else:
            matches = np.zeros(4).astype(bool)
            matches[0],_,_ = does_G_K_match(GAIAmag_dr2, Hmag, Kmag)
            matches[1],_,_ = does_GBP_K_match(GBPmag, Hmag, Kmag)
            matches[2],_,_ = does_GRP_K_match(GRPmag, Hmag, Kmag)
            matches[3],_,_ = does_GBP_GRP_match(GBPmag, GRPmag, Hmag, Kmag)
            match = np.all(matches)

        if match:
            dist, mu = compute_distance_modulus(unp.uarray(par,epar))
            l, b = coord.galactic.l.deg, coord.galactic.b.deg
	    AK = compute_AK_mwdust(l, b, unp.nominal_values(dist),
                                   unp.std_devs(dist))
            MK = compute_MK(unp.uarray(Kmag, e_Kmag), mu, AK)
            Rs = MK2Rs(MK)
            Teff = gaia2Teff(unp.uarray(GBPmag, eGBPmag),
                             unp.uarray(GRPmag, eGRPmag),
                             unp.uarray(Jmag, e_Jmag),
                             unp.uarray(Hmag, e_Hmag))
            Ms = MK2Ms(MK)
            return [par,epar], [GRPmag,eGRPmag], [GBPmag,eGBPmag], \
                [unp.nominal_values(dist), unp.std_devs(dist)], \
                [unp.nominal_values(mu), unp.std_devs(mu)], \
                [unp.nominal_values(AK), unp.std_devs(AK)], \
                [unp.nominal_values(MK), unp.std_devs(MK)], \
                [unp.nominal_values(Rs), unp.std_devs(Rs)], \
                [unp.nominal_values(Teff), unp.std_devs(Teff)], \
                [unp.nominal_values(Ms), unp.std_devs(Ms)], \
                  
    # return NaN if not found a match by now
    return np.repeat(np.nan,20).reshape(10,2)



def does_G_K_match(Gmag, Hmag, Kmag):
    '''Evans+2018 Tables A1 and A2
    (http://adsabs.harvard.edu/abs/2018A%26A...616A...4E)'''
    H_K = Hmag-Kmag
    if -.1 < H_K < .7:
        p = np.poly1d((-1.359, 12.0733, 0.6613))
        G_K = p(H_K)
        return np.isclose(G_K, Gmag-Kmag, atol=3*.3692), Gmag-Kmag, G_K
    else:
        return False, Gmag-Kmag, None


def does_GBP_K_match(GBPmag, Hmag, Kmag):
    '''Evans+2018 Tables A1 and A2
    (http://adsabs.harvard.edu/abs/2018A%26A...616A...4E)'''
    H_K = Hmag-Kmag
    if -.1 < H_K < .6:
        p = np.poly1d((2.711, 14.66, 0.8174))
        GBP_K = p(H_K)
        return np.isclose(GBP_K, GBPmag-Kmag, atol=3*.4839), GBPmag-Kmag, GBP_K
    else:
        return False, GBPmag-Kmag, None


def does_GRP_K_match(GRPmag, Hmag, Kmag):
    '''Evans+2018 Tables A1 and A2
    (http://adsabs.harvard.edu/abs/2018A%26A...616A...4E)'''
    H_K = Hmag-Kmag
    if -.1 < H_K < .7:
        p = np.poly1d((-1.721, 9.287, 0.3560))
        GRP_K = p(H_K)
        return np.isclose(GRP_K, GRPmag-Kmag, atol=3*.2744), GRPmag-Kmag, GRP_K
    else:
        return False, GRPmag-Kmag, None


def does_GBP_GRP_match(GBPmag, GRPmag, Hmag, Kmag):
    '''Evans+2018 Tables A1 and A2
    (http://adsabs.harvard.edu/abs/2018A%26A...616A...4E)'''
    H_K = Hmag-Kmag
    if -.1 < H_K < .55:
        p = np.poly1d((1.991, 6.098, 0.4238))
        GBP_GRP = p(H_K)
        return np.isclose(GBP_GRP, GBPmag-GRPmag, atol=3*.2144), \
            GBPmag-GRPmag, GBP_GRP
    else:
        return False, GBPmag-GRPmag, None

    

def compute_logg(Ms, Rs):
    G = 6.67e-11
    return unp.log10(G*rvs.Msun2kg(Ms)*1e2 / rvs.Rsun2m(Rs)**2)
    
    
def compute_distance_modulus(par_mas):
    if par_mas > 0:
        dist_pc = 1. / (par_mas*1e-3)
        mu = 5.*unp.log10(dist_pc) - 5
        return dist_pc, mu
    else:
        return unp.uarray(np.nan, np.nan), unp.uarray(np.nan, np.nan) 


def compute_AK_mwdust(ls, bs, dist, edist, eAK_frac=.3):
    '''Using the EB-V map from 2014MNRAS.443.2907S and the extinction vector
    RK = 0.31 from Schlafly and Finkbeiner 2011 (ApJ 737, 103)'''
    dustmap = mwdust.Combined15(filter='2MASS Ks')
    dist_kpc, edist_kpc = np.ascontiguousarray(dist)*1e-3, \
                          np.ascontiguousarray(edist)*1e-3
    ls, bs = np.ascontiguousarray(ls), np.ascontiguousarray(bs)
    AK, eAK = np.zeros(ls.size), np.zeros(ls.size)
    for i in range(ls.size):
        v = dustmap(ls[i], bs[i],
                    np.array([dist_kpc[i], dist_kpc[i]+edist_kpc[i]]))
        AK[i], eAK[i] = v[0], np.sqrt(abs(np.diff(v))**2 + (eAK_frac*v[0])**2)
    return unp.uarray(AK, eAK)


def compute_MK(Kmags, mus, AKs):
    return Kmags - mus - AKs


def MK2Rs(MK):
    '''Use relation from Mann+2015 (table 1)
    http://adsabs.harvard.edu/abs/2015ApJ...804...64M
    '''
    if 4.6 < float(unp.nominal_values(MK)) < 9.8:
        a, b, c, Rs_sigma_frac = 1.9515, -.3520, .01680, .0289
        p = np.poly1d((c,b,a))
        Rss = p(MK)
        eRss = np.sqrt(unp.std_devs(Rss)**2 + \
                       (unp.nominal_values(Rss)*Rs_sigma_frac)**2)
        return unp.uarray(unp.nominal_values(Rss), eRss)
    else:
        return unp.uarray(np.nan, np.nan)


def gaia2Teff(GBPmag, GRPmag, Jmag, Hmag):
    '''Use the relation from Mann+2015 (table 2)
    http://adsabs.harvard.edu/abs/2015ApJ...804...64M
    '''
    a, b, c, d, e, f, g = 3.172, -2.475, 1.082, -.2231, .01738, .08776, -.04355
    p1 = np.poly1d((e,d,c,b,a))
    p2 = np.poly1d((g,f,0))
    Teff = 35e2 * (p1(GBPmag-GRPmag) + p2(Jmag-Hmag))
    eTeff = np.sqrt(unp.std_devs(Teff)**2 + 49**2)
    return unp.uarray(unp.nominal_values(Teff), eTeff)


def MK2Ms(MK):
    '''Use relation from Benedict+2016
    http://arxiv.org/abs/1608.04775'''
    if 5 <= unp.nominal_values(MK) <= 10:
        c0, c1, c2, c3, c4, x0 = unp.uarray(.2311,4e-4), \
                                 unp.uarray(-.1352,7e-4), \
                                 unp.uarray(.04,5e-4), \
                                 unp.uarray(.0038,2e-4), \
                                 unp.uarray(-.0032,1e-4), 7.5
        dMK = MK - x0
        Ms = c0 + c1*dMK + c2*dMK**2 + c3*dMK**3 + c4*dMK**4
        return Ms
    else:
        return unp.uarray(np.nan, np.nan)
    
    

if __name__ == '__main__':
    t0 = time.time()
    fout = 'input_data/TESStargets/TESSMdwarfs_sector1.csv'
    tics = np.loadtxt('input_data/TESStargets/all_targets_S001_v1.csv',
                      delimiter=',', skiprows=6)[:,0]
    tics = tics[0:200000]
    get_stellar_data_TIC(tics, fout, overwrite=True)
    print 'Took %.3f min'%((time.time()-t0)/60.)
    send_email()
