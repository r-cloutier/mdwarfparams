from K2LCclass import *
from TESS_search import *
import linear_lnlike as llnl
import batman
import sys
from mcmc1 import run_emcee
from massradius import radF2mass
from truncate_cmap import *

global K2Mdwarffile
K2Mdwarffile = 'input_data/K2targets/K2Mdwarfsv5.csv'


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
    return name, Rs, bjd[g], f[g], ef[g]


def read_K2_data(epicnum):
    # make directories
    try:
	os.mkdir('MAST')
    except OSError:
	pass
    try:
	os.mkdir('MAST/K2')
    except OSError:
	pass

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
            folder2 = 'MAST/K2/EPIC%i'%epicnum
            try:
                os.mkdir(folder2)
            except OSError:
                pass
            os.system('mv %s %s'%(fname, folder2))
            break

    # read fits file
    hdu = fits.open('%s/%s'%(folder2, fname))
    assert len(hdu) > 1
     
    for i in range(1,len(hdu)):
        if hdu[i].header['EXTNAME'] == 'BESTAPER':
            bjd = hdu[i].data['T']+hdu[i].header['BJDREFI']
            f = hdu[i].data['FCOR']
            name = hdu[i].header['OBJECT'].replace(' ','_')
            star_dict = get_star(epicnum)
            # estimate flux error on 96 hour timescales
	    ef_6hrs = hdu[i].header['QCDPP6']*1e-6  # ppm to normalized flux
            ef_96hrs = ef_6hrs * np.sqrt(96./6)
            break

    ef = np.repeat(ef_96hrs, bjd.size)
    s = np.argsort(bjd)
    return name, star_dict, bjd[s], f[s], ef[s]


def get_star(epicnum):
    d = np.loadtxt(K2Mdwarffile, delimiter=',')
    epicnums = d[:,0]
    g = epicnums == epicnum
    assert g.sum() == 1
    star_info = d[g].reshape(23)
    star_dict = {'epicnum': int(star_info[0]), 'ra': star_info[1],
                 'dec': star_info[2], 'K2campaign': star_info[3],
                 'Kepmag': star_info[4], 'par': star_info[5],
                 'e_par': star_info[6], 'Kmag': star_info[7],
                 'e_Kmag': star_info[8], 'dist': star_info[9],
                 'e_dist': star_info[10], 'mu': star_info[11],
                 'e_mu': star_info[12], 'MK': star_info[13],
                 'e_MK': star_info[14], 'Rs': star_info[15],
                 'e_Rs': star_info[16], 'Teff': star_info[17],
                 'e_Teff': star_info[18], 'Ms': star_info[19],
                 'e_Ms': star_info[20], 'logg': star_info[21],
                 'e_logg': star_info[22]}
    return star_dict


def is_star_of_interest(epicnum):
    '''Return True is star obeys the desired conditions'''
    star_dict = get_star(epicnum)
    # 3083 K2 M dwarfs w/ Kepmag<15.2
    return (star_dict['Kepmag'] < 15.2) & (star_dict['Ms'] <= .75) & \
        (star_dict['Rs'] <= .75) & (star_dict['logg'] > 3) & \
        (star_dict['Teff'] >= 2700) & (star_dict['Teff'] <= 4000)


def planet_search(epicnum):
    '''Run a planet search on an input Kepler or K2 light curve using the 
    pipeline defined in compute_sensitivity to search for planets.'''
    
    # get data and only run the star if it is of interest
    if not is_star_of_interest(epicnum):
        return None
    name, star_dict, bjd, f, ef = read_K2_data(epicnum)

    # save stellar data and time-series
    self = K2LC(name, -99)  # -99 is unique to planet_search
    self.bjd, self.f, self.ef = bjd, f, ef
    for attr in star_dict.keys():
        setattr(self, attr, star_dict[attr])
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
                                                    thetaGPout)
    self.params_guess = params
    self.params_guess_labels = np.array(['Ps', 'T0s', 'depths [Z]', \
                                         'durations [D]'])
    self.EBparams_guess, self.maybeEBparams_guess = EBparams, maybeEBparams
    self._pickleobject()

    # fit joint model to light curve
    self.Ndet = self.params_guess.shape[0]
    self.params_optimized = np.zeros((self.Ndet, 5))
    self.params_optimized_labels = np.array(['P','T0','a/Rs','rp/Rs','inc'])
    self.params_results = np.zeros((self.Ndet, 3, 5))
    nwalkers, burnin, nsteps = 100, 200, 400
    self.params_samples = np.zeros((self.Ndet, nwalkers*nsteps, 5))
    for i in range(self.Ndet):
        _,_,_,_,_,theta = llnl.fit_params(self.params_guess[i], self.bjd,
                                          self.fcorr, self.ef, self.Ms,
                                          self.Rs, self.Teff)
        self.params_optimized[i] = theta[:5]
        self.u1, self.u2 = self.params_optimized[i,-2:]
        initialize = [1e-3, 1e-3, 1e-1, 1e-2, 1e-2]
        _,samples,results = run_emcee(self.params_optimized[i,:5], self.bjd,
                                      self.fcorr, self.ef, initialize, self.u1,
                                      self.u2, self.Ms, self.Rs, nwalkers=nwalkers,
                                      burnin=burnin, nsteps=nsteps, a=2)
        self.params_samples[i] = samples
        self.params_results[i] = results

    self.DONE = True
    self._pickleobject()
    


def do_i_run_this_star(epicnum):
    # first check that the star is available
    epics= np.loadtxt(K2Mdwarffile, delimiter=',')[:,0]
    g = epics == epicnum
    if g.sum() != 1:
	return False
    # check if star is already done
    fname = 'PipelineResults/EPIC_%i/K2LC'%epicnum
    if os.path.exists(fname):
        return not loadpickle(fname).DONE
    else:
        return True


if __name__ == '__main__':
    startind = int(sys.argv[1])
    endind = int(sys.argv[2])
    epics= np.loadtxt(K2Mdwarffile, delimiter=',')[:,0]
    #epics = np.loadtxt('input_data/K2targets/K2knownMdwarfplanets.csv', delimiter=',')
    for i in range(startind, endind):
	print epics[i]
	if do_i_run_this_star(epics[i]):
            planet_search(epics[i])
