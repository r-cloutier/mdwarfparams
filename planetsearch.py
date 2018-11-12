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


global K2Mdwarffile, threshBayesfactor
K2Mdwarffile = 'input_data/K2targets/K2Mdwarfsv7.csv'
KepMdwarffile = 'input_data/Keplertargets/KepMdwarfsv9.csv'
threshBayesfactor = 1e2


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
    fnames = np.array(glob.glob('%s/kplr*.fits'%folder2))
    if fnames.size == 0:
        dir2 = '%.9d'%Kepid
        dir1 = dir2[:4]
        url = 'https://archive.stsci.edu/pub/kepler/lightcurves/'
        url += '%s/%s'%(dir1,dir2)

        # get fits files in the directory and download the LCs
        fs = listFD(url, ext='.fits')
        for i in range(fs.size):
            os.system('wget %s'%fs[i])
            fname = str(fs[i]).split('/')[-1]
            if os.path.exists(fname):
                os.system('mv %s %s'%(fname, folder2))

    # get data from fits files
    fnames = np.array(glob.glob('%s/kplr*.fits'%folder2))
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

    g = np.where(ef_avg == np.min(ef_avg))[0]
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
    folder2 = 'MAST/K2/EPIC%i'%epicnum
    try:
       	os.mkdir(folder2)
    except OSError:
       	pass

    # download tar file if not already
    # from https://archive.stsci.edu/hlsps/k2sff/
    fnames = np.array(glob.glob('%s/hlsp*.fits'))
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
        hdu = fits.open('%s/%s'%(folder2, fname))
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


def get_star(IDnum, Kep=False, K2=False, TESS=False):
    if Kep:
        infile = KepMdwarffile
        IDtype, sector = 'kepid', 'N/A'
    elif K2:
        infile = K2Mdwarffile
        IDtype, sector = 'epicnum', 'K2campaign'
    elif TESS:
        infile = TESSMdwarffile
        IDtype, sector = 'ticid', 'TESSsector'
    else:
        raise ValueError('Must select one of Kep, K2, or TESS')

    d = np.loadtxt(infile, delimiter=',')
    # if Kepler, add a fake column to replace K2campaign/TESSsector
    if Kep:
	IDnums = d[:,0]
        g = IDnums == IDnum
        assert g.sum() == 1
        star_info = d[g].reshape(32)
	star_dict = {IDtype: int(star_info[0]), 'ra': star_info[1],
                     'dec': star_info[2], 'Kepmag': star_info[3],
                     'par': star_info[4], 'e_par': star_info[5],
                     'Jmag': star_info[6], 'e_Jmag': star_info[7],
                     'Hmag': star_info[8], 'e_Hmag': star_info[9],
                     'Kmag': star_info[10], 'e_Kmag': star_info[11],
                     'GBPmag': star_info[12], 'e_GBPmag': star_info[13],
                     'GRPmag': star_info[14], 'e_GRPmag': star_info[15],
                     'dist': star_info[16], 'e_dist': star_info[17],
                     'mu': star_info[18], 'e_mu': star_info[19],
                     'AK': star_info[20], 'e_AK': star_info[21],
                     'MK': star_info[22], 'e_MK': star_info[23],
                     'Rs': star_info[24], 'e_Rs': star_info[25],
                     'Teff': star_info[26], 'e_Teff': star_info[27],
                     'Ms': star_info[28], 'e_Ms': star_info[29],
                     'logg': star_info[30], 'e_logg': star_info[31]}
    elif K2:
	IDnums = d[:,0]
        g = IDnums == IDnum
        assert g.sum() == 1
        star_info = d[g].reshape(25)
    	star_dict = {IDtype: int(star_info[0]), 'ra': star_info[1],
                     'dec': star_info[2], sector: star_info[3],
                     'Kepmag': star_info[4], 'par': star_info[5],
                     'e_par': star_info[6], 'Kmag': star_info[7],
                     'e_Kmag': star_info[8], 'dist': star_info[9],
                     'e_dist': star_info[10], 'mu': star_info[11],
                     'e_mu': star_info[12], 'AK': star_info[13],
                     'e_AK': star_info[14], 'MK': star_info[15],
                     'e_MK': star_info[16], 'Rs': star_info[17],
                     'e_Rs': star_info[18], 'Teff': star_info[19],
                     'e_Teff': star_info[20], 'Ms': star_info[21],
                     'e_Ms': star_info[22], 'logg': star_info[23],
                     'e_logg': star_info[24]}
        
    return star_dict


def is_star_of_interest(IDnum, Kep=False, K2=False, TESS=False):
    '''Return True if star obeys the desired conditions'''
    star_dict = get_star(IDnum, Kep=Kep, K2=K2, TESS=TESS)
    return (star_dict['Ms'] <= .75) & \
        (star_dict['Rs'] <= .75) & (star_dict['logg'] > 3) & \
        (star_dict['Teff'] >= 2700) & (star_dict['Teff'] <= 4000)


def run_mcmc(self, nwalkers=100, burnin=200, nsteps=400):
    '''Run an MCMC on the detrended light curve for 1 planet models to obtain
    the posterior PDFs and model parameter uncertainties.'''
    self.Ndet = self.params_guess.shape[0]
    self.params_optimized = np.zeros((self.Ndet, 5))
    self.params_optimized_labels = np.array(['P','T0','a/Rs','rp/Rs','inc'])
    self.fmodels = np.zeros((self.Ndet, self.bjd.size))
    nwalkers, burnin, nsteps = int(nwalkers), int(burnin), int(nsteps)
    self.params_lnprobs = np.zeros((self.Ndet, 2, nwalkers*nsteps))
    #self.params_samples = np.zeros((self.Ndet, 2, nwalkers*nsteps, 5))
    self.params_results = np.zeros((self.Ndet, 3, 5))

    for i in range(self.Ndet):

        _,_,_,_,fmodel,theta = llnl.fit_params(self.params_guess[i], self.bjd,
                                               self.fcorr, self.ef, self.Ms,
                                               self.Rs, self.Teff)
        self.params_optimized[i] = theta[:5]
	self.fmodels[i] = fmodel
        self.u1, self.u2 = theta[-2:]

        # do 0 and 1-planet models
        for j in range(2):

            zeroplanetmodel = True if j == 0 else False        
            initialize = [1e-3, 1e-3, 1e-1, 1e-2, 1e-2]
            print 'Running MCMC on %i-planet model with P=%.3f days'%(j,
                                                                      theta[0])
            sampler,samples = run_emcee(self.params_optimized[i], self.bjd,
                                        self.fcorr, self.ef, initialize,
                                        self.u1, self.u2, self.Ms, self.Rs,
                                        a=2, nwalkers=nwalkers, burnin=burnin,
                                        nsteps=nsteps,
                                        zeroplanetmodel=zeroplanetmodel)
            #self.params_samples[i,j] = samples
            self.params_lnprobs[i,j] = sampler.lnprobability.flatten()
            if j == 1:
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
        

def planet_search(IDnum, Kep=False, K2=False, TESS=False):
    '''Run a planet search on an input Kepler or K2 light curve using the 
    pipeline defined in compute_sensitivity to search for planets.'''
    
    # get data and only run the star if it is of interest
    if not is_star_of_interest(IDnum, Kep=Kep):
        return None

    # get parameters depending on the inut light curves
    if Kep:
        Kepnum, Nopt = IDnum, 5
        name, star_dict, bjd, f, ef, quarters = read_Kepler_data(Kepnum)
        K2, TESS = False, False
    elif K2:
	epicnum, Nopt = IDnum, 10
        name, star_dict, bjd, f, ef, quarters = read_K2_data(epicnum)
        Kep, TESS = False, False
    elif TESS:
        TICnum, Nopt = IDnum, 10
        name, star_dict, bjd, f, ef, quarters = read_TESS_data(TICnum)
        Kep, K2 = False, False
    else:
        raise ValueError('Must select one of Kep, K2, or TESS')

    # stop if bad time-series (e.g. sometimes baseline is too short)
    if np.any(np.isnan(bjd)):
	return None

    # save stellar data and time-series
    self = LCclass(name, -99)  # -99 is unique to planet_search
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
    params, EBparams, maybeEBparams = find_transits(self, self.bjd, self.f,
                                                    self.ef, self.quarters,
                                                    self.thetaGPout)
    self.params_guess = params
    self.params_guess_labels = np.array(['Ps', 'T0s', 'depths [Z]', \
                                         'durations [D]'])
    self.EBparams_guess, self.maybeEBparams_guess = EBparams, maybeEBparams
    self._pickleobject()

    # fit transit models to light curve
    run_mcmc(self)

    # compute model evidences for 0-Ndet planet models
    #self.lnevidences = compute_model_evidences(self.params_samples,
    #                                           self.params_lnprobs)
    #self.evidence_ratios = compute_evidence_ratios(self.lnevidences)
    #self.params_guess_final = get_planet_detections(self.evidence_ratios,
    #                                                self.params_guess)
    #self.Ndet = self.params_guess_final.shape[0]
    self.DONE = True
    self._pickleobject()
    


def do_i_run_this_star(ID, K2=False, Kep=False):
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
    else:
        return None

    if g.sum() != 1:
	return False
    # check if star is already done
    fname = glob.glob('PipelineResults/%s_%i/LC_-00099'%(prefix, ID))
    if len(fname) == 1:
        return not loadpickle(fname[0]).DONE
    else:
        return True


if __name__ == '__main__':
    startind = int(sys.argv[1])
    endind = int(sys.argv[2])
    #epics= np.loadtxt(K2Mdwarffile, delimiter=',')[:,0]
    kepids = np.loadtxt(KepMdwarffile, delimiter=',')[:,0]
    for i in range(startind, endind):
	print kepids[i] #epics[i]
	if do_i_run_this_star(kepids[i], Kep=True):
            planet_search(kepids[i], Kep=True)
