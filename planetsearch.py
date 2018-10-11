from K2LCclass import *
from TESS_search import *
import linear_lnlike as llnl
import batman
import sys
from massradius import radF2mass
from truncate_cmap import *
#from joint_LCmodel import *

global K2Mdwarffile
K2Mdwarffile = 'input_data/K2targets/K2Mdwarfsv4.csv'

def read_K2_data_OLD(fits_path):
    '''Get one light curve and from the "best" aperture.'''
    hdu = fits.open(fits_path)
    for i in range(1,len(hdu)):
	if hdu[i].header['EXTNAME'] == 'BESTAPER':
    	    bjd, f = hdu[i].data['T']+hdu[i].header['BJDREFI'], hdu[i].data['FCOR']
	    name = hdu[i].header['OBJECT'].replace(' ','_')
            print name
            names,_,Rss,Mss,Teffs = np.genfromtxt('MAST/K2_epic', delimiter=',', dtype='|S50').T
            Rs = float(Rss[names==name].astype(float))
            Ms = float(Mss[names==name].astype(float))
            Teff = float(Teffs[names==name].astype(float))
	    # read or get flux error
	    effile = 'MAST/K2_eflux'
	    names, efs = np.genfromtxt(effile, delimiter=',',dtype='|S50').T
	    if name in names:	
		ef = np.repeat(efs[names==name].astype(float), bjd.size)
	    else:
		bad = True
		while bad:
		    plt.plot(np.arange(bjd.size), (f/np.nanmedian(f)-1)*1e6, '.'), plt.show()
		    start = int(raw_input('First index to estimate ef = '))
		    end = int(raw_input('Last index to estimate ef = '))
		    ef = float((f/np.nanmedian(f))[start:end].std())
		    bad = not bool(int(raw_input('Is %.4e ppm a reasonable ef value (1=yes, 0=no)? '%(ef*1e6))))
		h = open(effile, 'r')
		g = h.read()
		h.close()
		g += '%s,%6e\n'%(name, ef)
		h = open(effile, 'w')
		h.write(g)
		h.close()
		ef = np.repeat(ef, bjd.size)
    	    break
    s = np.argsort(bjd)
    return name, Rs, Ms, Teff, bjd[s], f[s], ef[s]


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
    campaigns = [0,1,2,3,4,5,6,7,8,102,111,112,12,13,14,15,16,17,91,92]  # from https://archive.stsci.edu/hlsps/k2sff/
    for j in range(len(campaigns)):
        folder = 'c%.2d/%.4d00000/%.5d'%(campaigns[j], int(str(epicnum)[:4]), int(str(epicnum)[4:9]))
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
            bjd, f = hdu[i].data['T']+hdu[i].header['BJDREFI'], hdu[i].data['FCOR']
            name = hdu[i].header['OBJECT'].replace(' ','_')
            Kepmag, logg, Ms, Rs, Teff = get_star(epicnum)
            # estimate flux error on 96 hour timescales
	    ef_6hrs = hdu[i].header['QCDPP6']  * 1e-6    # ppm to normalized flux
            ef_96hrs = ef_6hrs * np.sqrt(96./6)
            break

    ef = np.repeat(ef_96hrs, bjd.size)
    s = np.argsort(bjd)
    return name, Kepmag, logg, Ms, Rs, Teff, bjd[s], f[s], ef[s]


def get_star(epicnum):
    epics, Kepmags, Teffs, loggs,  Rss, Mss = np.loadtxt(K2Mdwarffile, delimiter=',').T
    g = epics == epicnum
    assert g.sum() == 1
    return float(Kepmags[g]), float(loggs[g]), float(Mss[g]), float(Rss[g]), float(Teffs[g])


def is_star_of_interest(epicnum):
    '''Return True is star obeys the desired conditions'''
    Kepmag, logg, Ms, Rs, Teff = get_star(epicnum)
    return (Kepmag<=14) & (Ms<=.75) & (Rs<=.75) & (logg>3) & (Teff>=2700) & (Teff<=4000)  # 3528 K2 M dwarfs w/ Kepmag<14
    #return (Ms<=.75) & (Rs<=.75) & (logg>3) & (Teff<=4000)


def planet_search(epicnum):
    '''Run a planet search on an input Kepler or K2 light curve using the pipeline defined in 
    compute_sensitivity to search for planets.'''
    
    # get data and only run the star if it is of interest
    if not is_star_of_interest(epicnum):
        return None
    name, Kepmag, logg, Ms, Rs, Teff, bjd, f, ef = read_K2_data(epicnum)
    self = K2LC(name)
    self.bjd, self.f, self.ef = bjd, f, ef
    self.Kepmag, self.logg, self.Ms, self.Rs, self.Teff = Kepmag, logg, Ms, Rs, Teff
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
