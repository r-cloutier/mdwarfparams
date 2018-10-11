from imports import *
import astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.gaia import Gaia
from astroquery.vizier import Vizier
import rvs
from uncertainties import unumpy as unp
from planetsearch import get_star


def get_stellar_data(epicnums, radius_arcsec=10):
    
    Nstars = epicnums.size
    ras, decs, Kepmags = np.zeros(Nstars), np.zeros(Nstars), np.zeros(0)
    pars, Kmags = np.zeros((Nstars,2)), np.zeros((Nstars,2))
    dists, mus = np.zeros((Nstars,2)), np.zeros((Nstars,2))
    MKs, Rss = np.zeros((Nstars,2)), np.zeros((Nstars,2))
    Teffs, Mss = np.zeros((Nstars,2)), np.zeros((Nstars,2))
    loggs = np.zeros((Nstars, 2))
    for i in range(epicnums.size):

        # download fits header and get coordinates
        ras[i], decs[i], Kepmags[i] = get_data_from_fits(epicnums[i])

        # search gaia and 2MASS until we find a likely match
        # based on photometry
        if np.isfinite(ras[i]) and np.isfinite(decs[i]):
            
            pars[i] = np.nan, np.nan 
            while (np.isnan(pars[i,0])) and (radius_arcsec < 3600):
                p = query_one_star(ras[i], decs[i], radius_arcsec=radius_arcsec)
                pars[i],Kmags[i],dists[i],mus[i],MKs[i],Rss[i],Teffs[i],Mss[i]=p
                loggs[i] = compute_logg(unp.uarray(Mss[i]), unp.uarray(Rss[i]))
                radius_arcsec += 2

        else:
            p = np.repeat(np.nan,16).reshape(8,2)
            pars[i],Kmags[i],dists[i],mus[i],MKs[i],Rss[i],Teffs[i],Mss[i] = p
            loggs[i] = compute_logg(unp.uarray(Mss[i]), unp.uarray(Rss[i]))

    # save results
    hdr = 'EPIC,ra_deg,dec_deg,Kepmag,parallax_mas,e_parallax,Kmag,e_Kmag,'+ \
          'dist_pc,e_dist,mu,e_mu,MK,e_MK,Rs_RSun,e_Rs,Teff_K,e_Teff,'+ \
          'Ms_MSun,e_Ms'
    outarr = np.array([epicnums, ras, decs, Kepmags, pars[:,0], pars[:,1], 
                       Kmags[:,0], Kmags[:,1], dists[:,0], dists[:,1],
                       mus[:,0], mus[:,1], MKs[:,0], MKs[:,1], Rss[:,0],
                       Rss[:,1], Teffs[:,0], Teffs[:,1], Mss[:,0], Mss[:,1]]).T

    fout = 'input_data/K2targets/K2Mdwarf_radii.csv'
    if os.path.exists(fout):
        inarr = np.loadtxt(fout, delimiter=',')
        assert inarr.shape[1] == outarr.shape[1]
        np.savetxt(fout, np.append(inarr, outarr, axis=0),
                   delimiter=',', header=hdr, fmt='%.8e')
    else:
        np.savetxt(fout, outarr, delimiter=',', header=hdr, fmt='%.8e')



def get_data_from_fits(epicnum):
    # download tar file
    # from https://archive.stsci.edu/hlsps/k2sff/
    campaigns = [0,1,2,3,4,5,6,7,8,102,111,112,12,13,14,15,16,17,91,92]
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
            return hdr['RA_OBJ'], hdr['DEC_OBJ'], hdr['KEPMAG']

    return np.repeat(np.nan, 3)


def query_one_star(ra_deg, dec_deg, radius_arcsec=10):
    
    # get gaia data
    coord = SkyCoord(ra=ra_deg, dec=dec_deg,
                     unit=(u.degree, u.degree), frame='icrs')
    rad = u.Quantity(float(radius_arcsec), u.arcsec)
    rG = Gaia.query_object_async(coordinate=coord, radius=rad)
    # return NaN if no GAIA star is found
    if len(rG) == 0:
        return np.repeat(np.nan, 16).reshape(8,2)

    # get 2MASS data
    r2 = Vizier.query_region(coord, radius=rad, catalog='II/246')
    # return NaN if no 2MASS star is found
    if len(r2) == 0:
        return np.repeat(np.nan, 16).reshape(8,2)

    # step through possible matches
    for i in range(len(rG)):
        for j in range(len(r2)):
            Gmag, GBPmag, GRPmag = rG['phot_g_mean_mag'][i], \
                                   rG['phot_bp_mean_mag'][i], \
                                   rG['phot_rp_mean_mag'][i]
            par, epar = rG['parallax'][i], rG['parallax_error'][i]
            Hmag, Kmag, eKmag = r2[j]['Hmag'][0], r2[j]['Kmag'][0], \
                                r2[j]['e_Kmag'][0]
            
            # check that GAIA and 2MASS photometry aproximately match
            matches = np.zeros(4).astype(bool)
            matches[0],_,_ = does_G_K_match(Gmag, Hmag, Kmag)
            matches[1],_,_ = does_GBP_K_match(GBPmag, Hmag, Kmag)
            matches[2],_,_ = does_GRP_K_match(GRPmag, Hmag, Kmag)
            matches[3],_,_ = does_GBP_GRP_match(GBPmag, GRPmag, Hmag, Kmag)
            match = np.all(matches)

            if match:
                dist, mu = compute_distance_modulus(unp.uarray(par,epar))
                MK = compute_MK(unp.uarray(Kmag, eKmag), mu)
                Rs = MK2Rs(MK)
                Teff = MK2Teff(MK)
                Ms = MK2Ms(MK)
                return [par,epar], [Kmag,eKmag], \
                    [unp.nominal_values(dist), unp.std_devs(dist)], \
                    [unp.nominal_values(mu), unp.std_devs(mu)], \
                    [unp.nominal_values(MK), unp.std_devs(MK)], \
                    [unp.nominal_values(Rs), unp.std_devs(Rs)], \
                    [unp.nominal_values(Teff), unp.std_devs(Teff)], \
                    [unp.nominal_values(Ms), unp.std_devs(Ms)], \
                  
    # return NaN if not found a match by now
    return np.repeat(np.nan,16).reshape(8,2)


def does_G_K_match(Gmag, Hmag, Kmag):
    '''Evans+2018 Tables A1 and A2 
    (http://adsabs.harvard.edu/abs/2018A%26A...616A...4E)'''
    H_K = Hmag-Kmag
    if -.1 < H_K < .7:
        p = np.poly1d((-1.359, 12.0733, 0.6613))
        G_K = p(H_K)
        return np.isclose(G_K, Gmag-Kmag, atol=2*.3692), Gmag-Kmag, G_K
    else:
        return False, Gmag-Kmag, G_K


def does_GBP_K_match(GBPmag, Hmag, Kmag):
    '''Evans+2018 Tables A1 and A2
    (http://adsabs.harvard.edu/abs/2018A%26A...616A...4E)'''
    H_K = Hmag-Kmag
    if -.1 < H_K < .6:
        p = np.poly1d((2.711, 14.66, 0.8174))
        GBP_K = p(H_K)
        return np.isclose(GBP_K, GBPmag-Kmag, atol=2*.4839), GBPmag-Kmag, GBP_K
    else:
        return False, GBPmag-Kmag, GBP_K


def does_GRP_K_match(GRPmag, Hmag, Kmag):
    '''Evans+2018 Tables A1 and A2
    (http://adsabs.harvard.edu/abs/2018A%26A...616A...4E)'''
    H_K = Hmag-Kmag
    if -.1 < H_K < .7:
        p = np.poly1d((-1.721, 9.287, 0.3560))
        GRP_K = p(H_K)
        return np.isclose(GRP_K, GRPmag-Kmag, atol=2*.2744), GRPmag-Kmag, GRP_K
    else:
        return False, GRPmag-Kmag, GRP_K

    
def does_GBP_GRP_match(GBPmag, GRPmag, Hmag, Kmag):
    '''Evans+2018 Tables A1 and A2
    (http://adsabs.harvard.edu/abs/2018A%26A...616A...4E)'''
    H_K = Hmag-Kmag
    if -.1 < H_K < .55:
        p = np.poly1d((1.991, 6.098, 0.4238))
        GBP_GRP = p(H_K)
        return np.isclose(GBP_GRP, GBPmag-GRPmag, atol=2*.2144), \
            GBPmag-GRPmag, GBP_GRP
    else:
        return False, GBPmag-GRPmag, GBP_GRP


def compute_logg(Ms, Rs):
    G = 6.67e-11
    return unp.log10(G*rvs.Msun2kg(Ms)*1e2 / rvs.Rsun2m(Rs)**2)
    
    
def compute_distance_modulus(pars_mas):
    dists_pc = 1. / (pars_mas*1e-3)
    mus = 5.*unp.log10(dists_pc) - 5
    return dists_pc, mus


def compute_MK(Kmags, mus):
    return Kmags - mus


def MK2Rs(MKs):
    '''Use relation from Mann+2015
    http://adsabs.harvard.edu/abs/2015ApJ...804...64M
    '''
    g = 4.6 < MKs < 9.8
    a, b, c, Rs_sigma_frac = 1.9515, -.3520, .01680, .0289
    p = np.poly1d((c,b,a))
    Rss = p(MKs)
    eRss = np.sqrt(unp.std_devs(Rss)**2 + \
                   (unp.nominal_values(Rss)*Rs_sigma_frac)**2)
    return unp.uarray(unp.nominal_values(Rss), eRss)


def MK2Ms(MKs):
    '''Use relation from Benedict+2016
    http://arxiv.org/abs/1608.04775'''
    c0, c1, c2, c3, c4, x0 = unp.uarray(.2311,4e-4), unp.uarray(-.1352,7e-4), \
                             unp.uarray(.04,5e-4), unp.uarray(.0038,2e-4), \
                             unp.uarray(-.0032,1e-4), 7.5
    p = np.poly1d((c4,c3,c2,c1,c0))
    Ms = p(MKs-x0)
    return Ms
    
    

if __name__ == '__main__':
    fs = np.array(glob.glob('MAST/K2/EPIC*'))
    epicnums = np.zeros(fs.size)
    for i in range(fs.size):
        epicnums[i] = int(fs[i].split('EPIC')[-1])
    
    #get_stellar_data(epicnums)
