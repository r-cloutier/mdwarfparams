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
KepMdwarffile = 'input_data/Keplertargets/KepMdwarfsv2.csv'
threshBayesfactor = 1e2


def read_Kepler_data(KICid):
    # make directories
    try:
	os.mkdir('MAST')
    except OSError:
	pass
    try:
	os.mkdir('MAST/Kepler')
    except OSError:
	pass
    folder2 = 'MAST/Kepler/KIC%.9d'%KICid
    try:
       	os.mkdir(folder2)
    except OSError:
       	pass

    # download LC files if not already
    fnames = np.array(glob.glob('%s/kplr*.fits'%folder2))
    if fnames.size == 0:
        dir2 = '%.9d'%KICid
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

    star_dict = get_star(KICid, Kep=True)
    g = (np.isfinite(bjd)) & (np.isfinite(f)) & (np.isfinite(ef))
    s = np.argsort(bjd[g])
    return 'KIC_%i'%KICid, star_dict, bjd[g][s], f[g][s], ef[g][s], quarters[g][s]


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


def get_star(IDnum, K2=False, Kep=False):
    infile = K2Mdwarffile if K2 else KepMdwarffile
    print infile
    d = np.loadtxt(infile, delimiter=',')
    IDnums = d[:,0]
    g = IDnums == IDnum
    assert g.sum() == 1
    # add fake column if Kepler (to replace K2campaign)
    if d.shape[1] == 24:
	d = np.insert(d, 3, np.zeros(d.shape[0]), axis=1)
    star_info = d[g].reshape(25)
    IDtype = 'epicnum' if K2 else 'kicid'
    star_dict = {IDtype: int(star_info[0]), 'ra': star_info[1],
                 'dec': star_info[2], 'K2campaign': star_info[3],
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


def is_star_of_interest(IDnum, Kep=False):
    '''Return True is star obeys the desired conditions'''
    K2 = not Kep
    star_dict = get_star(IDnum, K2=K2)
    # 3083 K2 M dwarfs w/ Kepmag<15.2
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
        

def planet_search(IDnum, K2=False, Kep=False):
    '''Run a planet search on an input Kepler or K2 light curve using the 
    pipeline defined in compute_sensitivity to search for planets.'''
    
    # get data and only run the star if it is of interest
    if not is_star_of_interest(IDnum, Kep=Kep):
        return None
    if K2:
	epicnum = IDnum
        name, star_dict, bjd, f, ef, quarters = read_K2_data(epicnum)
        Kepler = False
    elif Kep:
	KICnum = IDnum
        name, star_dict, bjd, f, ef, quarters = read_Kepler_data(KICnum)
    else:
        return None
        
    # save stellar data and time-series
    self = LCclass(name, -99)  # -99 is unique to planet_search
    self.bjd, self.f, self.ef, self.quarters = bjd, f, ef, quarters
    for attr in star_dict.keys():
        setattr(self, attr, star_dict[attr])
    self.DONE = False
    self._pickleobject()

    # fit initial GP hyperparams for each quarter
    thetaGPin, thetaGPout = do_optimize_0(bjd, f, ef, quarters)
    self.thetaGPin, self.thetaGPout = thetaGPin, thetaGPout
    self._pickleobject()
 
    # search for transits in the corrected LC and get the transit parameters
    # guesses
    print 'Searching for transit-like events...\n'
    params, EBparams, maybeEBparams = find_transits(self, bjd, f, ef,
                                                    quarters, thetaGPout)
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
        kicids = np.loadtxt(KepMdwarffile, delimiter=',')[:,0]
        prefix = 'KICid'
        g = kicids == ID
    else:
        return None

    if g.sum() != 1:
	return False
    # check if star is already done
    fname = glob.glob('PipelineResults/%s_%i/*LC_-00099'%(prefix, ID))
    if len(fname) == 1:
        return not loadpickle(fname[0]).DONE
    else:
        return True


if __name__ == '__main__':
    startind = int(sys.argv[1])
    endind = int(sys.argv[2])
    #epics= np.loadtxt(K2Mdwarffile, delimiter=',')[:,0]
    kicids = np.loadtxt(KepMdwarffile, delimiter=',')[:,0]
    for i in range(startind, endind):
	print kicids[i] #epics[i]
	if do_i_run_this_star(kicids[i], Kep=True):
            planet_search(kicids[i], Kep=True)
