# use the results from the Kepler-GAIA cross-match here: https://gaia-kepler.fun
from imports import *

def get_initial_KICs():
    '''Get the data for the Kepler-GAIA cross-matched stars.'''
    fs = np.sort(glob.glob('input_data/Keplertargets/kepler_dr2_*fits'))

    # define columns of interest
    cols = ['kepid','ra','dec','kepmag','parallax','parallax_error',
            'kmag',nan,'r_est','r_lo','r_hi']

    for i in range(
    outdata = 
    
    # get results from the 1 arcsec search radius
    f0 = fits.open(fs[0])[1].data
    out = np.array([f0.data['kepid'], f0.data['ra'], f0.data['dec'],
                    f0.data['kepmag'], f0.data['parallax']])
