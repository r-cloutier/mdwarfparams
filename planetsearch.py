from LCclass import *
from TESS_search import *
import linear_lnlike as llnl
import batman, sys
from mcmc1 import run_emcee
from priors import get_results
from massradius import radF2mass
from truncate_cmap import *
from perrakisFML import *
import requests
from bs4 import BeautifulSoup


global K2Mdwarffile, KepMdwarffile, TESSMdwarffile, threshBayesfactor
K2Mdwarffile = 'input_data/K2targets/K2lowmassstars_sens.csv'
KepMdwarffile = 'input_data/Keplertargets/KepMdwarfsv11.csv'
#TESSMdwarffile = 'input_data/TESStargets/TESSMdwarfs_sector1_v2.csv'
TESSMdwarffile = 'input_data/TESStargets/AUMic.csv'
#threshBayesfactor = 1e2


def read_Kepler_data(Kepid):
    # make directories
    try:
	os.mkdir('MAST')
    except OSError:
	pass
    try:
	os.mkdir('MAST/Kepler')
    except OSError:
	pass
    folder2 = 'MAST/Kepler/KepID%.9d'%Kepid
    try:
       	os.mkdir(folder2)
    except OSError:
       	pass

    # download LC files if not already
    fnames = np.array(glob.glob('%s/kplr*_llc.fits'%folder2))
    if fnames.size == 0:
        dir2 = '%.9d'%Kepid
        dir1 = dir2[:4]
        url = 'https://archive.stsci.edu/pub/kepler/lightcurves/'
        url += '%s/%s'%(dir1,dir2)

        # get fits files in the directory and download the 
 	# long cadence LCs
        fs = listFD(url, ext='_llc.fits')
        for i in range(fs.size):
            os.system('wget %s'%fs[i])
            fname = str(fs[i]).split('/')[-1]
            if os.path.exists(fname):
                os.system('mv %s %s'%(fname, folder2))

    # get data from fits files
    fnames = np.sort(np.array(glob.glob('%s/kplr*_llc.fits'%folder2)))
    bjd, f, ef, quarters = np.zeros(0), np.zeros(0), np.zeros(0), np.zeros(0)
    for i in range(fnames.size):
        hdu = fits.open(fnames[i])[1]
        bjd = np.append(bjd, hdu.data['TIME']+2454833)
        ftmp = hdu.data['PDCSAP_FLUX']
        eftmp = hdu.data['PDCSAP_FLUX_ERR']
        eftmp /= np.nanmedian(ftmp)
        ftmp /= np.nanmedian(ftmp)
        f = np.append(f, ftmp)
        ef = np.append(ef, eftmp)
        quarters = np.append(quarters, np.zeros(ftmp.size)+i)

    # limit Kepler baseline to save on computational expense
    print 'KepID ', Kepid
    bjd, f, ef, quarters = _reduce_Kepler_baseline(bjd, f, ef, quarters)
        
    star_dict = get_star(Kepid, Kep=True)
    return 'KepID_%i'%Kepid, star_dict, bjd, f, ef, quarters


def _reduce_Kepler_baseline(bjd, f, ef, quarters, tmax=270):
    '''Q0-Q17 baselines are ~1440 days (~3.9 years) so searching for individual 
    transits in the linear search phase is cripplingly expensive (compared to 
    K2 with 80 day baselines). For my M-dwarf occurrence rate calculations I am 
    only concered with P<=100 days. In a light curve with 100 ppm measurement 
    uncertainties, an Earth-sized planet around a 0.5 RSun star has a
    SNR = 3.35*sqrt(N_transits) so only 2-3 transits are needed for a detection 
    so the maximum baseline is set to 2-3 times Pmax or ~300 days. Let's make 
    it 270 since many Kepler quarters are ~90 days in length.'''  
    assert bjd.size == f.size
    assert bjd.size == ef.size
    assert bjd.size == quarters.size

    # clean
    g = (np.isfinite(bjd)) & (np.isfinite(f)) & (np.isfinite(ef))
    s = np.argsort(bjd[g])
    bjd, f, ef, quarters = bjd[g][s], f[g][s], ef[g][s], quarters[g][s]

    # get each quarters duration and median uncertainty
    if bjd.max()-bjd.min() < tmax:
	return np.repeat(np.nan, 4)

    NQ = np.unique(quarters).size
    med_ef_quarter = np.zeros(NQ)
    bjdmin_quarter = np.zeros(NQ)
    bjdmax_quarter = np.zeros(NQ)
    for i in range(NQ):
        g = quarters == i
        med_ef_quarter[i] = np.median(ef[g])
        bjdmin_quarter[i] = bjd[g].min()
        bjdmax_quarter[i] = bjd[g].max()

    # take the least noisy consecutive quarters that encapsulate tmax days
    tspan, med_ef_span = np.zeros((NQ,NQ)), np.zeros((NQ,NQ))
    for i in range(NQ):
        for j in range(NQ):
            tspan[i,j] = bjdmax_quarter[i] - bjdmin_quarter[j]
            med_ef_span[i,j] = med_ef_quarter[i]

    g_rows, g_cols = np.where(np.isclose(tspan, tmax, rtol=.1))
    if g_rows.size == 0:
        return np.repeat(np.nan, 4)
    ef_avg = np.zeros(g_rows.size)
    for i in range(g_rows.size):
        ef_avg[i] = med_ef_span[g_cols[i]:g_rows[i]+1,g_cols[i]].sum()

    g = np.where(ef_avg == np.min(ef_avg))[0][0]
    g1 = np.in1d(quarters, range(g_cols[g], g_rows[g]+1))
    bjd, f, ef, quarters = bjd[g1], f[g1], ef[g1], quarters[g1]
    quarters -= quarters.min()   # start indexing from 0 

    return bjd, f, ef, quarters

        

def listFD(url, ext=''):
    '''List the contents of an https directory.'''
    page = requests.get(url).text
    soup = BeautifulSoup(page, 'html.parser')
    return np.array([url + '/' + node.get('href') \
                     for node in soup.find_all('a') \
                     if node.get('href').endswith(ext)])
                    

def read_K2_data_k2sff(epicnum):
    # make directories
    try:
	os.mkdir('MAST')
    except OSError:
	pass
    try:
	os.mkdir('MAST/K2')
    except OSError:
	pass
    folder2 = 'MAST/K2/EPIC%i'%epicnum
    try:
       	os.mkdir(folder2)
    except OSError:
       	pass

    # download tar file if not already
    # from https://archive.stsci.edu/hlsps/k2sff/
    fnames = np.array(glob.glob('%s/hlsp*.fits'%folder2))
    if fnames.size == 0:
        campaigns=[0,1,2,3,4,5,6,7,8,91,92,102,111,112,12,13,14,15,16,17][::-1]
        for j in range(len(campaigns)):
            folder = 'c%.2d/%.4d00000/%.5d'%(campaigns[j],
                                             int(str(epicnum)[:4]),
                                             int(str(epicnum)[4:9]))
            fname = 'hlsp_k2sff_k2_lightcurve_%.9d-c%.2d_kepler_v1_llc.fits'%(int(epicnum), campaigns[j])
            url = 'https://archive.stsci.edu/hlsps/k2sff/%s/%s'%(folder, fname)
            os.system('wget %s'%url)
            if os.path.exists(fname):
                os.system('mv %s %s'%(fname, folder2))
                break

        # read fits file
        try:
	    hdu = fits.open('%s/%s'%(folder2, fname))
	except IOError:
	    return '', {}, np.zeros(0), np.zeros(0), np.zeros(0), np.zeros(0)
        assert len(hdu) > 1

    # file already downloaded
    else:
        hdu = fits.open(fnames[0])
        assert len(hdu) > 1

    # get data from fits file
    for i in range(1,len(hdu)):
        if hdu[i].header['EXTNAME'] == 'BESTAPER':
            bjd = hdu[i].data['T']+hdu[i].header['BJDREFI']
            f = hdu[i].data['FCOR']
            name = hdu[i].header['OBJECT'].replace(' ','_')
            star_dict = get_star(epicnum, K2=True)
            # estimate flux error on 96 hour timescales
	    ef_6hrs = hdu[i].header['QCDPP6']*1e-6  # ppm to normalized flux
            ef_96hrs = ef_6hrs * np.sqrt(96./6)
            break

    ef = np.repeat(ef_96hrs, bjd.size)
    s = np.argsort(bjd)
    quarters = np.zeros(bjd.size)
    return name, star_dict, bjd[s], f[s], ef[s], quarters[s]


def read_K2_data_everest(epicnum):
    # make directories
    try:
	os.mkdir('MAST')
    except OSError:
	pass
    try:
	os.mkdir('MAST/K2_EVEREST')
    except OSError:
	pass
    folder2 = 'MAST/K2_EVEREST/EPIC%i'%epicnum
    try:
       	os.mkdir(folder2)
    except OSError:
       	pass

    # download tar file if not already
    fnames = np.array(glob.glob('%s/hlsp*.fits'%folder2))
    if fnames.size == 0:
        campaigns=[0,1,2,3,4,5,6,7,8,102,111,112,12,13,14,15,16,17,18][::-1]
        for j in range(len(campaigns)):
            folder = 'c%.2d/%.4d00000/%.5d'%(campaigns[j],
                                             int(str(epicnum)[:4]),
                                             int(str(epicnum)[4:9]))
            fname = 'hlsp_everest_k2_llc_%.9d-c%.2d_kepler_v2.0_lc.fits'%(int(epicnum), campaigns[j])
            url = 'https://archive.stsci.edu/hlsps/everest/v2/%s/%s'%(folder, fname)
            os.system('wget %s'%url)
            if os.path.exists(fname):
                os.system('mv %s %s'%(fname, folder2))
                break

        # read fits file
        try:
	    hdu = fits.open('%s/%s'%(folder2, fname))
	except IOError:
	    return '', {}, np.zeros(0), np.zeros(0), np.zeros(0), np.zeros(0)
        assert len(hdu) > 1

    # file already downloaded
    else:
        hdu = fits.open(fnames[0])
        assert len(hdu) > 1

    # get data from fits file
    g = hdu[1].data['QUALITY'] == 0
    bjd = hdu[1].data['TIME'][g] + hdu[1].header['BJDREFI']
    f = hdu[1].data['FCOR'][g]
    ef = hdu[1].data['FRAW_ERR'][g]
    name = hdu[1].header['OBJECT'].replace(' ','_')
    star_dict = get_star(epicnum, K2=True)
    
    s = np.argsort(bjd)
    quarters = np.zeros(bjd.size)
    return name, star_dict, bjd[s], f[s], ef[s], quarters[s]



def read_TESS_data(tic):
    # make directories
    try:
	os.mkdir('MAST')
    except OSError:
	pass
    try:
	os.mkdir('MAST/TESS')
    except OSError:
	pass
    folder2 = 'MAST/TESS/TIC%i'%tic
    try:
       	os.mkdir(folder2)
    except OSError:
       	pass

    # setup directories to fetch the data
    tic_str = '%.16d'%int(tic)
    tid1 = tic_str[:4]
    tid2 = tic_str[4:8]
    tid3 = tic_str[8:12]
    tid4 = tic_str[12:]

    # download tar file if not already
    fnames = np.array(glob.glob('%s/tess*.fits'%folder2))
    if fnames.size == 0:
        sectors = [1,2]
        for j in range(len(sectors)):
            sctr = 's%.4d'%sectors[j]
            url = 'https://archive.stsci.edu/missions/tess/tid/'
            folder = '%s/%s/%s/%s/%s/'%(sctr,tid1,tid2,tid3,tid4)
    	    print url+folder
            fs = listFD(url+folder, ext='_lc.fits')
            for i in range(fs.size):
                os.system('wget %s'%fs[i])
                fname = str(fs[i]).split('/')[-1]
                if os.path.exists(fname):
                    os.system('mv %s %s'%(fname, folder2))

    # get data from fits files
    fnames = np.sort(np.array(glob.glob('%s/tess*_lc.fits'%folder2)))
    bjd, f, ef, quarters = np.zeros(0), np.zeros(0), np.zeros(0), np.zeros(0)
    for i in range(fnames.size):
        hdus = fits.open(fnames[i])
        bjd = np.append(bjd, hdus[1].data['TIME'] + 2457000)
        ftmp = hdus[1].data['PDCSAP_FLUX']
        eftmp = hdus[1].data['PDCSAP_FLUX_ERR']
        eftmp /= np.nanmedian(ftmp)
        ftmp /= np.nanmedian(ftmp)
        f = np.append(f, ftmp)
        ef = np.append(ef, eftmp)
        quarters = np.append(quarters, np.zeros(ftmp.size)+i)

    # get stellar data
    name = hdus[0].header['OBJECT'].replace(' ','_')
    star_dict = get_star(tic, TESS=True)
    s = np.argsort(bjd)
    bjd, f, ef, quarters = bjd[s], f[s], ef[s], quarters[s]
    g = np.isfinite(bjd) & np.isfinite(f) & np.isfinite(ef) 
    return name, star_dict, bjd[g], f[g], ef[g], quarters[g]



def get_star(IDnum, Kep=False, K2=False, TESS=False):
    if Kep:
        infile = KepMdwarffile
        IDtype, magname = 'kepid', 'Kepmag'
    elif K2:
        infile = K2Mdwarffile
        IDtype, magname = 'epicnum', 'Kepmag'
    elif TESS:
        infile = TESSMdwarffile
        IDtype, magname = 'tic', 'TESSmag'
    else:
        raise ValueError('Must select one of Kep, K2, or TESS')

    d = np.loadtxt(infile, delimiter=',')
    IDnums = d[:,0]
    g = IDnums == IDnum
    assert g.sum() == 1
    star_info = d[g].reshape(39)
    star_dict = {IDtype: int(star_info[0]), 'ra': star_info[1],
                 'dec': star_info[2], 'GBPmag': star_info[3],
                 'e_GBPmag': star_info[4], 'GRPmag': star_info[5],
                 'e_GRPmag': star_info[6], magname: star_info[7],
                 'Jmag': star_info[8], 'e_Jmag': star_info[9],
                 'Hmag': star_info[10], 'e_Hmag': star_info[11],
                 'Kmag': star_info[12], 'e_Kmag': star_info[13],
                 'par': star_info[14], 'e_par': star_info[15],
                 'dist': star_info[16], 'ehi_dist': star_info[17],
                 'elo_dist': star_info[18], 'mu': star_info[19],
                 'ehi_mu': star_info[20], 'elo_mu': star_info[21],
                 'AK': star_info[22], 'e_AK': star_info[23],
                 'MK': star_info[24], 'ehi_MK': star_info[25],
                 'elo_MK': star_info[26], 'Rs': star_info[27],
                 'ehi_Rs': star_info[28], 'elo_Rs': star_info[29],
                 'Teff': star_info[30], 'ehi_Teff': star_info[31],
                 'elo_Teff': star_info[32], 'Ms': star_info[33],
                 'ehi_Ms': star_info[34], 'elo_Ms': star_info[35],
                 'logg': star_info[36], 'ehi_logg': star_info[37],
                 'elo_logg': star_info[38]}
    
    return star_dict


def is_star_of_interest(IDnum, Kep=False, K2=False, TESS=False):
    '''Return True if star obeys the desired conditions'''
    star_dict = get_star(IDnum, Kep=Kep, K2=K2, TESS=TESS)
    return (star_dict['Ms']-star_dict['elo_Ms'] <= .75) & \
        (star_dict['Rs']-star_dict['elo_Rs'] <= .75) & \
	(star_dict['logg']+star_dict['ehi_logg'] > 3.5) & \
        (star_dict['Teff']+star_dict['ehi_Teff'] >= 2700) & \
	(star_dict['Teff']-star_dict['elo_Teff'] <= 4000)
    #print 'TEMP: Only by-passed for AU Mic'
    return True


def run_mcmc(self, nwalkers=100, burnin=200, nsteps=400, Kep=False, TESS=False):
    '''Run an MCMC on the detrended light curve for 1 planet models to obtain
    the posterior PDFs and model parameter uncertainties.'''
    self.Ndet = self.params_guess.shape[0]
    self.params_optimized = np.zeros((self.Ndet, 5))
    self.params_optimized_labels = np.array(['P','T0','a/Rs','rp/Rs','inc'])
    self.fmodels = np.zeros((self.Ndet, self.bjd.size))
    nwalkers, burnin, nsteps = int(nwalkers), int(burnin), int(nsteps)
    self.params_lnprobs = np.zeros((self.Ndet, nwalkers*nsteps))
    #self.params_samples = np.zeros((self.Ndet, nwalkers*nsteps, 5))
    self.params_results = np.zeros((self.Ndet, 3, 5))

    for i in range(self.Ndet):

        _,_,_,_,fmodel,theta = llnl.fit_params(self.params_guess[i], self.bjd,
                                               self.fcorr, self.ef, self.Ms,
                                               self.Rs, self.Teff, Kep=Kep,
					       TESS=TESS)
        self.params_optimized[i] = theta[:5]
	self.fmodels[i] = fmodel
        u1, u2 = theta[-2:]

	# run MCMC on transit LC
        initialize = [1e-3, 1e-3, 1e-1, 1e-2, 1e-2]
        print 'Running MCMC on 1-planet model with P=%.3f days'%(theta[0])
        sampler,samples = run_emcee(self.params_optimized[i], self.bjd,
                                    self.fcorr, self.ef, initialize,
                                    u1, u2, self.Ms, self.Rs,
                                    a=2, nwalkers=nwalkers, burnin=burnin,
                                    nsteps=nsteps, zeroplanetmodel=False)
        #self.params_samples[i] = samples
        self.params_lnprobs[i] = sampler.lnprobability.flatten()
        results = get_results(samples)
        self.params_results[i] = results


def compute_model_evidences(samples, lnprobs):
    '''Compute the model evidences for 0 and 1-planet models for each detected
    planet.'''
    Nplanets = samples.shape[0]
    assert samples.shape[:2] == (Nplanets, 2)
    assert lnprobs.shape[:2] == (Nplanets, 2)

    # compute lnZ for 0 and 1-planet models for each detected planet 
    lnevidences = np.zeros((Nplanets, 2))
    for i in range(Nplanets):
        for j in range(2):
            print j
            g = np.invert(np.all(samples[i,j]==0, axis=0))
            lnevidences[i,j] = compute_perrakis_FML(samples[i,j,:,g].T,
                                                    lnprobs[i,j])
        
    return lnevidences


def compute_evidence_ratios(lnZs):
    '''Compute the model evidence ratio for 0 to 1-planets models for each 
    detected planet (crudely) assuming uniform model priors.'''
    prior_ratio_1_0 = 1.
    Nplanets = lnZs.shape[0]
    evidence_ratios = np.zeros(Nplanets)
    for i in range(Nplanets):
        evidence_ratios[i] = np.exp(lnZs[i,1] - lnZs[i,0]) * prior_ratio_1_0

    return evidence_ratios


def get_planet_detections(evidence_ratios, params_guess):
    '''Flag planets as either detected or otherwise based on their model 
    evidences and the threshold evidence (i.e. a free parameter).'''
    Nplanets = evidence_ratios.size
    assert params_guess.shape[0] == Nplanets

    detected = np.zeros(Nplanets).astype(bool)
    detected[evidence_ratios >= threshBayesfactor] = True
    params_guess_out = params_guess[order][detected]

    s = np.argsort(params_guess_out[:,0])
    return params_guess_out[s]



def compute_depths(transit_params, u1, u2):
    assert len(transit_params.shape) == 2
    Np = transit_params.shape[0]
    func = llnl.transit_model_func_curve_fit(u1,u2)
    depths = np.zeros(Np)
    for i in range(Np):
        P,T0 = transit_params[i,:2]
        tarr = np.linspace(T0-P/5, T0+P/5, 1000)
        fmodel = func(tarr, *transit_params[i])
        depths[i] = 1. - np.nanmin(fmodel)
    return depths



def planet_search(folder, IDnum, Kep=False, K2=False, TESS=False):
    '''Run a planet search on an input Kepler or K2 light curve using the 
    pipeline defined in compute_sensitivity to search for planets.'''
    
    # get data and only run the star if it is of interest
    if not is_star_of_interest(IDnum, Kep=Kep, K2=K2, TESS=TESS):
	print 'star not of interest: ', IDnum
        return None

    # get parameters depending on the input light curves
    if Kep:
        Kepnum, Nopt = IDnum, 5
        name, star_dict, bjd, f, ef, quarters = read_Kepler_data(Kepnum)
        K2, TESS = False, False
    elif K2:
	epicnum, Nopt = IDnum, 10
        name, star_dict, bjd, f, ef, quarters = read_K2_data_k2sff(epicnum)
        Kep, TESS = False, False
    elif TESS:
        TICnum, Nopt = IDnum, 10
        name, star_dict, bjd, f, ef, quarters = read_TESS_data(TICnum)
        Kep, K2 = False, False
    else:
        raise ValueError('Must select one of Kep, K2, or TESS')

    # stop if bad time-series (e.g. sometimes baseline is too short)
    if np.any(np.isnan(bjd)):
        print 'bad bjds'
        return None

    # save stellar data and time-series
    self = LCclass(folder, name, -99)  # -99 is unique to planet_search
    self.bjd, self.f, self.ef, self.quarters = bjd, f, ef, quarters
    for attr in star_dict.keys():
        setattr(self, attr, star_dict[attr])
    self.DONE = False
    self._pickleobject()

    # fit initial GP hyperparams for each quarter
    thetaGPin, thetaGPout = do_optimize_0(bjd, f, ef, quarters, N=Nopt)
    self.thetaGPin, self.thetaGPout = thetaGPin, thetaGPout
    self._pickleobject()

    # search for transits in the corrected LC and get the transit parameters
    # guesses
    print '\nSearching for transit-like events...\n'
    if K2 or Kep:
        Kep, TESS = True, False
	self.u1, self.u2 = llnl.get_LDcoeffs_Kepler(self.Ms, self.Rs, self.Teff)
    else:
        Kep, TESS = False, True
        self.u1, self.u2 = llnl.get_LDcoeffs_TESS(self.Ms, self.Rs, self.Teff)
    params, EBparams, maybeEBparams = find_transits(self, self.bjd, self.f,
                                                    self.ef, self.quarters,
                                                    self.thetaGPout, Kep=Kep,
						    TESS=TESS)
    self.params_guess = params
    self.params_guess_labels = np.array(['Ps', 'T0s', 'depths [Z]', \
                                         'durations [D]'])
    self.EBparams_guess, self.maybeEBparams_guess = EBparams, maybeEBparams
    self._pickleobject()

    # fit transit models to light curve and compute diagnostics
    run_mcmc(self, Kep=Kep, TESS=TESS)

    # compute transit S/N of planet candidates
    Np = self.params_optimized.shape[0]
    self.CDPPs = np.array([llnl.estimate_CDPP(self.bjd,self.fcorr,self.ef,
                                              self.params_guess[i,3])
                           for i in range(Np)])
    self.depths = compute_depths(self.params_optimized, self.u1, self.u2)
    self.Ntransits=np.array([llnl.compute_Ntransits(self.bjd,
                                                    *self.params_optimized[i,:2]) 
			       for i in range(Np)])
    self.SNRtransits = self.depths / self.CDPPs * np.sqrt(self.Ntransits)

    # compute model evidences for 0-Ndet planet models
    #self.lnevidences = compute_model_evidences(self.params_samples,
    #                                           self.params_lnprobs)
    #self.evidence_ratios = compute_evidence_ratios(self.lnevidences)
    #self.params_guess_final = get_planet_detections(self.evidence_ratios,
    #                                                self.params_guess)
    #self.Ndet = self.params_guess_final.shape[0]
    self.DONE = True
    self._pickleobject()
    


def do_i_run_this_star(folder, ID, K2=False, Kep=False, TESS=False):
    # first check that the star is available
    if K2:
        epics = np.loadtxt(K2Mdwarffile, delimiter=',')[:,0]
        g = epics == ID
        prefix = 'EPIC'
        Kep = False
    elif Kep:
        kepids = np.loadtxt(KepMdwarffile, delimiter=',')[:,0]
        prefix = 'KepID'
        g = kepids == ID
    elif TESS:
	tics = np.loadtxt(TESSMdwarffile, delimiter=',')[:,0]
	prefix = 'TIC'
	g = tics == ID
    else:
        return None

    if g.sum() != 1:
	return False
    # check if star is already done
    fname = glob.glob('%s/%s_%i/LC_-00099'%(folder, prefix, ID))
    if len(fname) == 1:
        return not loadpickle(fname[0]).DONE
    else:
        return True


if __name__ == '__main__':
    startind = int(sys.argv[1])
    endind = int(sys.argv[2])
    folder = sys.argv[3]
    if 'KepID' in folder:
	Kep, K2, TESS = True, False, False
	fname = KepMdwarffile
    if 'EPIC' in folder:
	Kep, K2, TESS = False, True, False
	fname = K2Mdwarffile
    if 'TIC' in folder:
	Kep, K2, TESS = False, False, True
	fname = TESSMdwarffile
    
    ids = np.loadtxt(fname, delimiter=',')[:,0]
    for i in range(startind, endind):
	print ids[i]
	if do_i_run_this_star(folder, ids[i], Kep=Kep, K2=K2, TESS=TESS):
            planet_search(folder, ids[i], Kep=Kep, K2=K2, TESS=TESS)
