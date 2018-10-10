from imports import *
import astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.gaia import Gaia
from astroquery.vizier import Vizier
import rvs
from uncertainties import unumpy as unp


def get_gaia_data(epicnums, radius_arcsec=10):

    Nstars = epicnums.size
    ras, decs = np.zeros(Nstars), np.zeros(Nstars)
    pars, Kmags = range(Nstars), range(Nstars)
    for i in range(epicnums.size):

        # get coordinates from fits header
        f = np.array(glob.glob('MAST/K2/EPIC%i/*.fits'%epicnums[i]))
        try:
            hdr = fits.open(f[0])[0].header
        except IOError:
            pass

        ras[i], decs[i] = hdr['RA_OBJ'], hdr['DEC_OBJ']

        # query gaia for the parallax
        radius_arcsec_orig, p = radius_arcsec+0, None
        while p == None and radius_arcsec < 3600:
            p, ep = one_gaia_query(ras[i], decs[i], radius_arcsec=radius_arcsec)
            radius_arcsec += 1
        pars[i] = unp.uarray(p, ep) if p != None else unp.uarray(np.nan, np.nan)

        # query 2MASS for the photometry
        radius_arcsec, Kmag = radius_arcsec_orig, None
        while Kmag == None and radius_arcsec < 3600:
            Kmag, eKmag = one_2MASS_query(ras[i], decs[i],
                                          radius_arcsec=radius_arcsec)
            radius_arcsec += 1
        Kmags[i] = unp.uarray(Kmag, eKmag) if Kmag != None else unp.uarray(np.nan, np.nan)
            
    # compute quantities of interest
    pars = np.ascontiguousarray(pars)
    Kmags = np.ascontiguousarray(Kmags)
    dists, mus = compute_distance_modulus(pars)
    MKs = compute_MK(Kmags, mus)
    Rss = MK2Rs(MKs)

    # save results
    hdr = 'EPIC,ra_deg,dec_deg,parallax_mas,e_parallax,Kmag,e_Kmag,dist_pc,e_dist,mu,e_mu,MK,e_MK,Rs_RSun,e_Rs'
    outarr = np.array([epicnums, ras, dec, unp.nominal_values(pars), unp.std_devs(pars),
		       unp.nominal_values(Kmags), unp.std_devs(Kmags), unp.nominal_values(dists), 
		       unp.std_devs(dists), unp.nominal_values(mus), unp.std_devs(mus),
		       unp.nominal_values(MKs), unp.std_devs(MKs), unp.nominal_values(Rss), 
		       unp.std_devs(Rss)]).T
    np.savetxt('input_data/K2targets/K2Mdwarf_radii.csv', outarr, delimiter=',', hdr=hdr, fmt='%.8e')


    
def one_gaia_query(ra_deg, dec_deg, radius_arcsec=10):
    coord = SkyCoord(ra=ra_deg, dec=dec_deg,
                     unit=(u.degree, u.degree), frame='icrs')
    width = u.Quantity(float(radius_arcsec), u.arcsec)
    height = u.Quantity(float(radius_arcsec), u.arcsec)
    r = Gaia.query_object_async(coordinate=coord, width=width, height=height)
    pars, epars = r['parallax'], r['parallax_error']
    if pars.size == 0:
        return None, None
    else:
        i = np.where(pars > 0)[0][0]
        return pars[i], epars[i]


    
def one_2MASS_query(ra_deg, dec_deg, radius_arcsec=10):
    '''Query 2MASS catalog from Cutri+2003 
    http://adsabs.harvard.edu/abs/2003tmc..book.....C'''
    coord = SkyCoord(ra=ra_deg, dec=dec_deg,
                     unit=(u.degree, u.degree), frame='icrs')
    width = u.Quantity(float(radius_arcsec), u.arcsec)
    height = u.Quantity(float(radius_arcsec), u.arcsec)
    r = Vizier.query_region(coord, width=width, height=height, catalog='II/246')
    i = 0
    Kmag, eKmag =  r[i]['Kmag'][0], r[i]['e_Kmag'][0]
    return Kmag, eKmag



def compute_distance_modulus(pars):
    dists_pc = 1. / (pars_mas*1e-3)
    mus = 5.*unp.log10(dists_pc) - 5
    return dists_pc, mus


def compute_MK(Kmags, mus):
    return Kmags - mus


def MK2Rs(MKs):
    '''Use relation from Mann+2015
    http://adsabs.harvard.edu/abs/2015ApJ...804...64M
    '''
    a, b, c, Rs_sigma_frac = 1.9515, -.3520, .01680, .0289
    p = np.poly1d((c,b,a))
    Rss = p(MKs)
    return unp.uarray(Rss, Rss*Rs_sigma_frac)
