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


def get_stellar_data_K2(epicnums, fout, radius_arcsec=10, overwrite=False):
    
    Nstars = epicnums.size
    radius_arcsec_orig = radius_arcsec+0
    ras, decs = np.zeros(Nstars), np.zeros(Nstars)
    Kepmags, K2campaigns = np.zeros(Nstars), np.zeros(Nstars)
    pars, Kmags = np.zeros((Nstars,2)), np.zeros((Nstars,2))
    dists, mus = np.zeros((Nstars,2)), np.zeros((Nstars,2))
    AKs, MKs = np.zeros((Nstars,2)), np.zeros((Nstars,2))
    Teffs, Mss = np.zeros((Nstars,2)), np.zeros((Nstars,2))
    Rss, loggs = np.zeros((Nstars, 2)), np.zeros((Nstars, 2))
    for i in range(epicnums.size):

        print float(i) / epicnums.size
        
        # download fits header and get coordinates
        ras[i],decs[i],Kepmags[i],K2campaigns[i] = \
                                            get_data_from_fits_K2(epicnums[i])
        
        # search gaia and 2MASS until we find a likely match
        # based on photometry
        if np.isfinite(ras[i]) and np.isfinite(decs[i]):

            pars[i] = np.nan, np.nan
            radius_arcsec = radius_arcsec_orig
            while (np.isnan(pars[i,0])) and (radius_arcsec <= 60):
                p = query_one_star(ras[i], decs[i], radius_arcsec=radius_arcsec)
                pars[i],Kmags[i],dists[i],mus[i],AKs[i],MKs[i],Rss[i],Teffs[i],Mss[i]=p
                logg = compute_logg(unp.uarray(Mss[i,0],Mss[i,1]),
                                    unp.uarray(Rss[i,0],Rss[i,1]))
                loggs[i] = unp.nominal_values(logg), unp.std_devs(logg)
                radius_arcsec += 5

        else:
            p = np.repeat(np.nan,18).reshape(9,2)
            pars[i],Kmags[i],dists[i],mus[i],AKs[i],MKs[i],Rss[i],Teffs[i],Mss[i] = p
            logg = compute_logg(unp.uarray(Mss[i,0],Mss[i,1]),
                                unp.uarray(Rss[i,0],Rss[i,1]))
            loggs[i] = unp.nominal_values(logg), unp.std_devs(logg)
            
    # save results
    hdr = 'EPIC,ra_deg,dec_deg,campaign,Kepmag,parallax_mas,e_parallax,Kmag,'+ \
          'e_Kmag,dist_pc,e_dist,mu,e_mu,AK,e_AK,MK,e_MK,Rs_RSun,e_Rs,'+ \
          'Teff_K,e_Teff,Ms_MSun,e_Ms,logg_dex,e_logg'
    outarr = np.array([epicnums, ras, decs, K2campaigns, Kepmags, pars[:,0],
                       pars[:,1], Kmags[:,0], Kmags[:,1], dists[:,0],
                       dists[:,1], mus[:,0], mus[:,1], AKs[:,0], AKs[:,1],
                       MKs[:,0], MKs[:,1], Rss[:,0], Rss[:,1], Teffs[:,0],
                       Teffs[:,1], Mss[:,0],Mss[:,1], loggs[:,0], loggs[:,1]]).T

    if os.path.exists(fout) and not overwrite:
        inarr = np.loadtxt(fout, delimiter=',')
        assert inarr.shape[1] == outarr.shape[1]
        np.savetxt(fout, np.append(inarr, outarr, axis=0),
                   delimiter=',', header=hdr, fmt='%.8e')
    else:
        np.savetxt(fout, outarr, delimiter=',', header=hdr, fmt='%.8e')


        
def get_stellar_data_Kep(KICids, fout, radius_arcsec=10, overwrite=False):
    
    Nstars = KICids.size
    radius_arcsec_orig = radius_arcsec+0
    ras, decs = np.zeros(Nstars), np.zeros(Nstars)
    Kepmags = np.zeros(Nstars)
    pars, Kmags = np.zeros((Nstars,2)), np.zeros((Nstars,2))
    dists, mus = np.zeros((Nstars,2)), np.zeros((Nstars,2))
    AKs, MKs = np.zeros((Nstars,2)), np.zeros((Nstars,2))
    Teffs, Mss = np.zeros((Nstars,2)), np.zeros((Nstars,2))
    Rss, loggs = np.zeros((Nstars, 2)), np.zeros((Nstars, 2))
    for i in range(KICids.size):

        print float(i) / KICids.size
        
        # download fits header and get coordinates
        ras[i],decs[i],Kepmags[i] = get_data_from_fits_Kep(KICids[i])
        
        # search gaia and 2MASS until we find a likely match
        # based on photometry
        if np.isfinite(ras[i]) and np.isfinite(decs[i]):

            pars[i] = np.nan, np.nan
            radius_arcsec = radius_arcsec_orig
            while (np.isnan(pars[i,0])) and (radius_arcsec <= 60):
                p = query_one_star(ras[i], decs[i], radius_arcsec=radius_arcsec)
                pars[i],Kmags[i],dists[i],mus[i],AKs[i],MKs[i],Rss[i],Teffs[i],Mss[i]=p
                logg = compute_logg(unp.uarray(Mss[i,0],Mss[i,1]),
                                    unp.uarray(Rss[i,0],Rss[i,1]))
                loggs[i] = unp.nominal_values(logg), unp.std_devs(logg)
                radius_arcsec += 5

        else:
            p = np.repeat(np.nan,18).reshape(9,2)
            pars[i],Kmags[i],dists[i],mus[i],AKs[i],MKs[i],Rss[i],Teffs[i],Mss[i] = p
            logg = compute_logg(unp.uarray(Mss[i,0],Mss[i,1]),
                                unp.uarray(Rss[i,0],Rss[i,1]))
            loggs[i] = unp.nominal_values(logg), unp.std_devs(logg)
            
    # save results
    hdr = 'KICid,ra_deg,dec_deg,Kepmag,parallax_mas,e_parallax,'+ \
          'Kmag,e_Kmag,dist_pc,e_dist,mu,e_mu,AK,e_AK,MK,e_MK,Rs_RSun,e_Rs,'+ \
          'Teff_K,e_Teff,Ms_MSun,e_Ms,logg_dex,e_logg'
    outarr = np.array([KICids, ras, decs, Kepmags, pars[:,0],
                       pars[:,1], Kmags[:,0], Kmags[:,1], dists[:,0],
                       dists[:,1], mus[:,0], mus[:,1], AKs[:,0], AKs[:,1],
                       MKs[:,0], MKs[:,1], Rss[:,0], Rss[:,1], Teffs[:,0],
                       Teffs[:,1], Mss[:,0],Mss[:,1], loggs[:,0], loggs[:,1]]).T

    if os.path.exists(fout) and not overwrite:
        inarr = np.loadtxt(fout, delimiter=',')
        assert inarr.shape[1] == outarr.shape[1]
        np.savetxt(fout, np.append(inarr, outarr, axis=0),
                   delimiter=',', header=hdr, fmt='%.8e')
    else:
        np.savetxt(fout, outarr, delimiter=',', header=hdr, fmt='%.8e')



def get_data_from_fits_K2(epicnum):
    # download tar file
    # from https://archive.stsci.edu/hlsps/k2sff/
    campaigns = [0,1,2,3,4,5,6,7,8,91,92,102,111,112,12,13,14,15,16,17][::-1]
    for j in range(len(campaigns)):
        folder = 'c%.2d/%.4d00000/%.5d'%(campaigns[j],
                                         int(str(epicnum)[:4]),
                                         int(str(epicnum)[4:9]))
        fname = 'hlsp_k2sff_k2_lightcurve_%.9d-c%.2d_kepler_v1_llc.fits'%(int(epicnum), campaigns[j])
        url = 'https://archive.stsci.edu/hlsps/k2sff/%s/%s'%(folder, fname)
        os.system('wget %s'%url)
        if os.path.exists(fname):
            hdr = fits.open(fname)[0].header
            os.system('rm %s'%fname)
	    ra, dec, campaign = hdr['RA_OBJ'], hdr['DEC_OBJ'], hdr['CAMPAIGN']
	    Kepmag = np.nan if isinstance(hdr['KEPMAG'], fits.card.Undefined) \
                     else hdr['KEPMAG']
            return ra, dec, Kepmag, campaign

    return np.repeat(np.nan, 4)


def get_data_from_fits_Kep(KICid):
    # download tar file
    dir2 = '%.9d'%KICid
    dir1 = dir2[:4]
    url = 'https://archive.stsci.edu/pub/kepler/lightcurves/'
    url += '%s/%s'%(dir1,dir2)

    # get fits files in the directory and download the LCs
    fs = listFD(url, ext='.fits')
    if fs.size > 0:
        os.system('wget %s'%fs[0])
        fname = str(fs[0]).split('/')[-1]
        if os.path.exists(fname):
            hdr = fits.open(fname)[0].header
            os.system('rm %s'%fname)
            ra, dec, Kepmag = hdr['RA_OBJ'], hdr['DEC_OBJ'], hdr['KEPMAG']
            return ra, dec, Kepmag
        else:
            return np.repeat(np.nan, 3)
    else:
        return np.repeat(np.nan, 3)



def query_one_star(ra_deg, dec_deg, radius_arcsec=10):
    
    # get gaia data
    coord = SkyCoord(ra=ra_deg, dec=dec_deg,
                     unit=(u.degree, u.degree), frame='icrs')
    rad = u.Quantity(float(radius_arcsec), u.arcsec)
    rG = Gaia.query_object_async(coordinate=coord, radius=rad)
    # return NaN if no GAIA star is found
    if len(rG) == 0:
        return np.repeat(np.nan, 18).reshape(9,2)

    # get 2MASS data
    r2 = Vizier.query_region(coord, radius=rad, catalog='II/246')
    # return NaN if no 2MASS star is found
    if (r2 == None) or (len(r2) == 0):
        return np.repeat(np.nan, 18).reshape(9,2)

    # step through possible matches
    for i in range(len(rG)):
        for j in range(len(r2)):
            
            # get giaa parameters
            par, epar = rG['parallax'][i]+.03, rG['parallax_error'][i]
            Gmag = rG['phot_g_mean_mag'][i]
            GBPmag = rG['phot_bp_mean_mag'][i]
            FBP = rG['phot_bp_mean_flux'][i]
            eFBP = rG['phot_bp_mean_flux_error'][i]
            eGBPmag = -2.5*np.log10(FBP / (FBP+eFBP))
            GRPmag = rG['phot_rp_mean_mag'][i]
            FRP = rG['phot_rp_mean_flux'][i]
            eFRP = rG['phot_rp_mean_flux_error'][i]
            eGRPmag = -2.5*np.log10(FRP / (FRP+eFRP))
            
            # get 2MASS photometry
            Hmag, Kmag, eKmag = r2[j]['Hmag'][0], r2[j]['Kmag'][0], \
                                r2[j]['e_Kmag'][0]

            # change bad values to NaN
            if np.ma.is_masked(par) or par<=0: par, epar = np.nan, np.nan 
            if np.ma.is_masked(epar) or epar<=0: par, epar = np.nan, np.nan 
            if np.ma.is_masked(Gmag): Gmag = np.nan
            if np.ma.is_masked(GBPmag): GBPmag, eGBPmag = np.nan, np.nan
            if np.ma.is_masked(eGBPmag): GBPmag, eGBPmag = np.nan, np.nan
            if np.ma.is_masked(GRPmag): GRPmag, eGRPmag = np.nan, np.nan
            if np.ma.is_masked(eGRPmag): GRPmag, eGRPmag = np.nan, np.nan
            if np.ma.is_masked(Hmag): Hmag = np.nan
            if np.ma.is_masked(Kmag): Kmag, eKmag = np.nan, np.nan
            if np.ma.is_masked(eKmag): Kmag, eKmag = np.nan, np.nan
            
            # check that GAIA and 2MASS photometry aproximately match
            if np.ma.is_masked(par) or np.ma.is_masked(Gmag):
                match = False
            else:
                matches = np.zeros(4).astype(bool)
                matches[0],_,_ = does_G_K_match(Gmag, Hmag, Kmag)
                matches[1],_,_ = does_GBP_K_match(GBPmag, Hmag, Kmag)
                matches[2],_,_ = does_GRP_K_match(GRPmag, Hmag, Kmag)
                matches[3],_,_ = does_GBP_GRP_match(GBPmag, GRPmag, Hmag, Kmag)
                match = np.all(matches)

            if match:
                dist, mu = compute_distance_modulus(unp.uarray(par,epar))
                print unp.nominal_values(dist)
                l, b = coord.galactic.l.deg, coord.galactic.b.deg
		AK = compute_AK_mwdust(l, b, unp.nominal_values(dist),
                                       unp.std_devs(dist))
                MK = compute_MK(unp.uarray(Kmag, eKmag), mu, AK)
                Rs = MK2Rs(MK)
                Teff = gaia2Teff(unp.uarray(GBPmag, eGBPmag),
                                 unp.uarray(GRPmag, eGRPmag))
                Ms = MK2Ms(MK)
                return [par,epar], [Kmag,eKmag], \
                    [unp.nominal_values(dist), unp.std_devs(dist)], \
                    [unp.nominal_values(mu), unp.std_devs(mu)], \
                    [unp.nominal_values(AK), unp.std_devs(AK)], \
                    [unp.nominal_values(MK), unp.std_devs(MK)], \
                    [unp.nominal_values(Rs), unp.std_devs(Rs)], \
                    [unp.nominal_values(Teff), unp.std_devs(Teff)], \
                    [unp.nominal_values(Ms), unp.std_devs(Ms)], \
                  
    # return NaN if not found a match by now
    return np.repeat(np.nan,18).reshape(9,2)


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
        print ls[i], bs[i], dist_kpc[i], edist_kpc[i]
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


def gaia2Teff(GBPmag, GRPmag):
    '''Use the relation from Mann+2015 (table 2)
    http://adsabs.harvard.edu/abs/2015ApJ...804...64M
    '''
    a, b, c, d, e  = 3.245, -2.4309, 1.043, -.2127, .01649
    p = np.poly1d((e,d,c,b,a))
    Teff = 35e2 * p(GBPmag-GRPmag)
    eTeff = np.sqrt(unp.std_devs(Teff)**2 + 55**2)
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
    #fname = 'input_data/K2targets/K2Mdwarfsv1_leq3300.csv' 
    #epicnums = np.loadtxt(fname, delimiter=',')[:,0]
    #epicnums = epicnums[1000:]
    #t0 = time.time()
    #fout = 'input_data/K2targets/K2Mdwarfsv6_leq3300.csv'
    #get_stellar_data_K2(epicnums, fout, overwrite=False)
    #print 'Took %.3f min'%((time.time()-t0)/60.)
    #send_email()

    fname = 'input_data/Keplertargets/kic10_search.csv'
    fout = 'input_data/Keplertargets/KepMdwarfsv1.csv'
    KICids = np.loadtxt(fname, delimiter=',', usecols=range(1))

    KICids = KICids[1001:10000]
    t0 = time.time()
    get_stellar_data_Kep(KICids, fout, overwrite=False)
    print 'Took %.3f min'%((time.time()-t0)/60.)
    send_email()
