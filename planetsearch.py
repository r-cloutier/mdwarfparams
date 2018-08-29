from KepK2LCclass import *
#from TESS_search import *
#import linear_lnlike as llnl
#import batman, ellc
#from massradius import radF2mass
#from truncate_cmap import *
#from joint_LCmodel import *
from imports import *


def read_K2_data(fits_path):
    '''Get one light curve and from the "best" aperture.'''
    hdu = fits.open(fits_path)
    for i in range(1,len(hdu)):
	if hdu[i].header['EXTNAME'] == 'BESTAPER':
    	    bjd, f = hdu[i].data['T']+hdu[i].header['BJDREFI'], hdu[i].data['FCOR']
	    name = hdu[i].header['OBJECT']
	    break
    ef = np.zeros(bjd.size)
    s = np.argsort(bjd)
    return name, bjd[s], f[s], ef[s]


def read_Kepler_data(fits_dir, maxdays=2e2):
    '''Combine the LCs from each Kepler quarter into a single LC.'''
    fs = np.array(glob.glob('%s/kplr*llc.fits'%fits_dir))
    assert fs.size > 1
    name = fits.open(fs[0])[0].header['OBJECT']
    bjd, f, ef = np.zeros(0), np.zeros(0), np.zeros(0)
    for i in range(fs.size):
  	hdu = fits.open(fs[i])
   	assert len(hdu) == 3
 	bjd = np.append(bjd, hdu[1].data['TIME']+hdu[1].header['BJDREFI'])
	ftmp = hdu[1].data['PDCSAP_FLUX']
	f = np.append(f, ftmp-np.nanmedian(ftmp))
	ef = np.append(ef, hdu[1].data['PDCSAP_FLUX_ERR'])
    # remove bad values
    g = np.isfinite(f) & np.isfinite(ef)
    bjd, f, ef = bjd[g], f[g], ef[g]
    # sort chronologically
    s = np.argsort(bjd)
    bjd, f, ef = bjd[s], f[s], ef[s]
    # resistrict to a certain time period
    g = bjd-bjd.min() <= maxdays
    return name, bjd[g], f[g], ef[g]


def planet_search(fits_path):
    '''Run a planet search on an input Kepler or K2 light curve using the pipeline defined in 
    compute_sensitivity to search for planets.'''
    
    # get data and only run the star if it is of interest
    name, bjd, f, ef = read_K2_data(fits_path) if 'fits' in fits_path else read_Kepler_data(fits_path)
    self = KepK2LC(name)
    self.bjd, self.f, self.ef = bjd, f, ef
    self.DONE = False
    self._pickleobject()

    # fit initial GP
    thetaGPall, resultsGPall, thetaGPin, thetaGPout = do_optimize_0(bjd, f, ef)
    self.thetaGPall, self.resultsGPall = thetaGPall, resultsGPall
    self.thetaGPin, self.thetaGPout = thetaGPin, thetaGPout
    self._pickleobject()
 
    # search for transits in the corrected LC and get the transit parameters
    # guesses
    print 'Searching for transit-like events...\n'
    params, EBparams, maybeEBparams = find_transits(self, bjd, f, ef,
                                                    thetaGPout, hdr,
                                                    fname_short)
    self.params_guess = params
    self.params_guess_labels = np.array(['Ps', 'T0s', 'depths [Z]', \
                                         'durations [D]'])
    self.EBparams_guess, self.maybeEBparams_guess = EBparams, maybeEBparams
    self._pickleobject()

    # are planets detected?
    '''self.is_detected = np.array([int(np.any(np.isclose(params[:,0], Ps[i],
                                                       rtol=.02)))
                                 for i in range(Ps.size)])
    assert tesslc.is_detected.size == tesslc.params_true.shape[0]
    tesslc.is_FP = np.array([int(np.invert(np.any(np.isclose(tesslc.params_guess[i,0],Ps,rtol=.02))))
                             for i in range(tesslc.params_guess.shape[0])])
    assert tesslc.is_FP.size == tesslc.params_guess.shape[0]
    tesslc.paramsFP_guess = tesslc.params_guess[tesslc.is_FP.astype(bool)]
    tesslc.pickleobject()

    # do joint GP+transit model
    if np.any(tesslc.is_detected.astype(bool)):
        params, transit_params, resultsGPfin, mufin, sigfin, LDcoeffs = \
                                                            joint_LC_fit(tesslc)
        tesslc.params_guessfin = params
        tesslc.transit_params = transit_params
        tesslc.u1, tesslc.u2 = LDcoeffs
        tesslc.resultsGPfin = resultsGPfin
        tesslc.mufin, tesslc.sigfin = mufin, sigfin
    '''
    self.DONE = True
    self.pickleobject()
    

if __name__ == '__main__':
    fits_path = sys.argv[1]
    planet_search(fits_path)
