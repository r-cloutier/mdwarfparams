from imports import *
from periodogram import compute_LSperiodogram
import GPmcmc0 as mcmc0
import GPmcmcN as mcmcN
import linear_lnlike as llnl
import rvs
from LCclass import *
from scipy.stats import normaltest
from scipy.signal import medfilt

global SNRthresh
SNRthresh = 5

def download_one_fits(fname):
    '''Get the full fits file for a TOI.'''
    urlprefix = 'https://archive.stsci.edu/missions/tess/ete-6/tid/00/000/000/057/'
    hdu = fits.open(urlprefix+fname)
    try:
    	hdu.writeto('MAST_Data/%s'%fname)
    except IOError:
	pass
    return hdu


def get_lc(hdu):
    '''Get the light curve from a fits file.'''
    hdr = hdu[0].header
    t, f, ef = hdu[1].data['TIME'], \
	       hdu[1].data['PDCSAP_FLUX'], \
               hdu[1].data['PDCSAP_FLUX_ERR']
    g = np.isfinite(t) & np.isfinite(f) & np.isfinite(ef)
    norm = np.median(f[g])
    t, f, ef = t[g], f[g]/norm, ef[g]/norm
    return hdr, t, f, ef


def initialize_GP_hyperparameters(bjd, f, ef, Pindex=0,
                                  Npntsmin=5e2, Npntsmax=1e3):
    '''Guess initial values of the GP hyperparameters.'''
    # P and l from periodogram
    per,_,pwrn = compute_LSperiodogram(bjd, f, ef,
                                       plims=(1.1,bjd.max()-bjd.min()))
    ##plt.plot(per, pwrn, 'k-'), plt.xscale('log'), plt.show()
    Pgp, inds, i = 1, np.argsort(pwrn)[::-1], int(Pindex)
    # ensure Pgp is not within 5% of 1 or 2 days
    while np.isclose(Pgp,1,rtol=.05) or np.isclose(Pgp,2,rtol=.05):
	Pgp = per[inds[i]]
	i += 1
    lnl, lnG, lnP = np.log(Pgp*3), 0., np.log(Pgp)
    # s is just a fraction of the photometric uncertainty
    s = np.median(ef)*.1
    # a from binned LC    
    Npnts_per_timescale = 8.
    timescale_to_resolve = Pgp / Npnts_per_timescale
    # bin the light curve
    Ttot = bjd.max() - bjd.min()
    if Ttot/timescale_to_resolve < Npntsmin:
        dt = Ttot / Npntsmin
    elif Ttot/timescale_to_resolve > Npntsmax:
        dt = Ttot / Npntsmax
    else:
        dt = timescale_to_resolve
    tbin, fbin, efbin = llnl.boxcar(bjd, f, ef, dt=dt)
    lna = np.log(np.max(abs(fbin-fbin.mean())) * .75)
    return lna, lnl, lnG, lnP#, s


def do_optimize_0(bjd, f, ef, quarters, N=10,
                  Npntsmin=5e2, Npntsmax=1e3, Nsig=3, medkernel=9):
    '''First fit the PDC LC with GP and a no planet model using an optimization 
    routine and tested with N different initializations.'''
    # fit GP separately to each quarter (or K2 campaign)
    assert bjd.size == f.size
    assert bjd.size == ef.size
    assert bjd.size == quarters.size
    NGP = np.unique(quarters).size
    
    # test various hyperparameter initializations and keep the one resulting
    # in the most gaussian-like residuals
    N = int(N)
    thetaGPs_in_tmp, thetaGPs_out_tmp = np.zeros((NGP,N,4)), \
                                        np.zeros((NGP,N,4))
    thetaGPs_in, thetaGPs_out = np.zeros((NGP,4)), np.zeros((NGP,4))
    for i in range(NGP):
        g1 = quarters == i
       
        pvalues = np.zeros(N) 
        for j in range(N):
            # get initial GP parameters (P, and l from periodogram)
            thetaGPs_in_tmp[i,j] = initialize_GP_hyperparameters(bjd[g1], f[g1],
                                                                 ef[g1],
                                                                 Pindex=j)
            Npnts_per_timescale = 8.
            inds = np.array([1,3])
            timescale_to_resolve = np.exp(thetaGPs_in_tmp[i,j,inds]).min() / \
                                   Npnts_per_timescale
            # bin the light curve
            Ttot = bjd[g1].max() - bjd[g1].min()
            if Ttot/timescale_to_resolve < Npntsmin:
                dt = Ttot / Npntsmin
            elif Ttot/timescale_to_resolve > Npntsmax:
                dt = Ttot / Npntsmax
            else:
                dt = timescale_to_resolve

            # trim outliers and median filter to avoid fitting deep transits
            g = abs(f[g1]-np.median(f[g1])) <= Nsig*np.std(f[g1])
            tbin, fbin, efbin = llnl.boxcar(bjd[g1][g], medfilt(f[g1][g],medkernel),
                                            ef[g1][g], dt=dt)
            gp,mu,sig,thetaGPs_out_tmp[i,j] = fit_GP_0(thetaGPs_in_tmp[i,j],
                                                       tbin, fbin, efbin)
            # compute residuals and the normality test p-value
            _,pvalues[j] = normaltest(fbin-mu)
            
        # select the most gaussian-like residuals
        g = np.argsort(pvalues)[::-1][0]
        thetaGPs_in[i]  = thetaGPs_in_tmp[i,g]
        thetaGPs_out[i] = thetaGPs_out_tmp[i,g]

    return thetaGPs_in, thetaGPs_out


def fit_GP_0(thetaGP, tbin, fbin, efbin):
    '''optimize the hyperparameters of this quasi-periodic GP'''
    assert len(thetaGP) == 4
    a, l, G, Pgp = np.exp(thetaGP)
    k1 = george.kernels.ExpSquaredKernel(l)
    k2 = george.kernels.ExpSine2Kernel(G,Pgp)
    gp = george.GP(a*(k1+k2))
    bnds = ((-20,0),(-3,10),(-5,5),(-3,10))
    results = gp.optimize(tbin, fbin, efbin, **{'bounds':bnds})
    try:
        gp.compute(tbin, efbin)
    except (ValueError, np.linalg.LinAlgError):
        return np.repeat(None,4)
    mu, cov = gp.predict(fbin, tbin)
    sig = np.sqrt(np.diag(cov))
    return gp, mu, sig, results[0]
    


def do_mcmc_N(thetaGP, params, bjd, f, ef, Nmcmc_pnts=3e2,
              nwalkers=100, burnin=200, nsteps=400, a=2):
    '''Fit the LC with a GP and a full transit model.'''
    Ntransits = params.size / 4
    theta = params.reshape(Ntransits*4)
    assert thetaGP.size == 5
    initialize = np.array([.1,.1,.01,.1,thetaGP[4]*.1])
    for i in range(Ntransits):
	initialize = np.append(initialize, [.1,.1,params[i,2]*.1,params[i,3]*.1])
    assert 0 not in initialize
    Prot = np.exp(thetaGP[3])
    if (Prot/4. > (bjd.max()-bjd.min())/1e2) or (np.arange(bjd.min(), bjd.max(), Prot/4.).size > Nmcmc_pnts):
        dt = (bjd.max()-bjd.min())/Nmcmc_pnts 
    else: 
        dt = Prot/4.
    tbin, fbin, efbin = llnl.boxcar(bjd,f,ef,dt=dt)
    theta_full = np.append(thetaGP, theta)
    sampler, samples = mcmcN.run_emcee(theta_full, tbin, fbin, efbin,
                                       initialize, nwalkers=nwalkers,
                                       burnin=burnin, nsteps=nsteps, a=a)
    results = mcmcN.get_results(samples)
    return sampler, samples, results


def save_fits(arr, fname):
    hdu = fits.PrimaryHDU(arr)
    hdu.writeto(fname, clobber=True)


def _optimize_GP(thetaGP, x, res, ey):
    '''Optimize the GP parameters of a mean-subtracted time-series and return 
    the mean GP model.'''
    assert len(thetaGP) == 4
    a, l, G, Pgp = np.exp(thetaGP)
    k1 = george.kernels.ExpSquaredKernel(l)
    k2 = george.kernels.ExpSine2Kernel(G,Pgp)
    gp = george.GP(a*(k1+k2))
    try:
	results = gp.optimize(x, res, ey)
	results = results[0]
    	gp.compute(x, ey)
    	mu, cov = gp.predict(res, x)
    	sig = np.sqrt(np.diag(cov))
    except (ValueError, np.linalg.LinAlgError):
	gp, results, mu, sig = None, np.zeros(len(thetaGP)), np.zeros(x.size), \
                               np.zeros(x.size)
    return gp, results, mu, sig


def _get_GP(thetaGP, x, res, ey):
    '''Compute the GP model.'''
    assert len(thetaGP) == 4
    a, l, G, Pgp = np.exp(thetaGP)
    k1 = george.kernels.ExpSquaredKernel(l)
    k2 = george.kernels.ExpSine2Kernel(G,Pgp)
    gp = george.GP(a*(k1+k2))
    try:
        gp.compute(x, ey)
        mu, cov = gp.predict(res, x)
        sig = np.sqrt(np.diag(cov))
        results = thetaGP
    except (ValueError, np.linalg.LinAlgError):
        gp, results, mu, sig = None, np.zeros(len(thetaGP)), np.zeros(x.size), \
                               np.zeros(x.size)
    return gp, results, mu, sig


def find_transits(self, bjd, f, ef, quarters, thetaGPs,
                  Npntsmin=5e2, Npntsmax=1e3, medkernel=99, Nsig=3,
		  Plims=(.5,1e2)):
    '''Search for periodic transit-like events.'''
    assert not np.all(ef == 0)
    assert np.unique(quarters).size == thetaGPs.shape[0]

    # "detrend" the lc
    detrend_LC(self, bjd, f, ef, quarters, thetaGPs, Npntsmin, Npntsmax,
               Nsig, medkernel)
    assert self.ef.mean() > 0

    # do linear search first
    print 'Computing lnL over transit times and durations...\n'
    self._pickleobject()
    bjd, fcorr, ef = self.bjd, self.fcorr, self.ef
    transit_times, durations, lnLs, depths = llnl.linear_search(bjd, fcorr, ef)
    self.transit_times, self.durations = transit_times, durations
    self.lnLs_linearsearch, self.depths_linearsearch = lnLs, depths
    self.SNRs_linearsearch = (lnLs-np.median(lnLs, axis=0)) / llnl.MAD2d(lnLs)
    self._pickleobject()
    
    # get transit candidates and initial parameters guesses
    print 'Computing lnL over periods and mid-transit times...\n'
    Ps, T0s, Ds, Zs, lnLs_transit = llnl.compute_transit_lnL(bjd, fcorr, ef,
                                                             transit_times,
                                                             durations, lnLs,
                                                             depths, SNRthresh)
    # set period limit
    Pmin, Pmax = Plims
    g = (Ps >= Pmin) & (Ps <= Pmax)
    Ps, T0s, Ds, Zs, lnLs_transit = Ps[g], T0s[g], Ds[g], Zs[g], lnLs_transit[g]

    # save here for debugging
    self.POIsorig, self.T0OIsorig, self.DOIsorig = Ps, T0s, Ds
    self.ZOIsorig, self.lnLOIsorig = Zs, lnLs_transit
    self._pickleobject()

    print 'Finding transit-like events and making transit parameter guesses...\n'
    POIs, T0OIs, DOIs, ZOIs, lnLOIs, params, EBparams, maybeEBparams = \
                        llnl.identify_transit_candidates(self, Ps, T0s, Ds, Zs,
                                                         lnLs_transit,
                                                         durations.size,
                                                         self.Rs, bjd, fcorr,
                                                         ef)
    self.POIs, self.T0OIs, self.DOIs, self.ZOIs, self.lnLOIs = POIs, T0OIs, \
                                                               DOIs, ZOIs, \
                                                               lnLOIs

    # return the parameters (P,T0,D,Z=depth) of the confirmed transit candidates
    return params, EBparams, maybeEBparams
 


def detrend_LC(self, bjd, f, ef, quarters, thetaGPs, Npntsmin, 
	       Npntsmax, Nsig, medkernel):
    assert thetaGPs.shape[1] == 4
    NGP = thetaGPs.shape[0]

    self.tbin, self.fbin, self.efbin = np.zeros(0), np.zeros(0), np.zeros(0) 
    mu, sig = np.zeros(bjd.size), np.zeros(bjd.size)
    fcorr = np.zeros(bjd.size)
    for i in range(NGP):
        
        Prot = np.exp(thetaGPs[i,3])
        Npnts_per_timescale, inds = 8., np.array([1,3])
        timescale_to_resolve = np.exp(thetaGPs[i,inds]).min() / \
                               Npnts_per_timescale

        # bin the light curve
        g1 = quarters == i
        Ttot = bjd[g1].max() - bjd[g1].min()
        if Ttot/timescale_to_resolve < Npntsmin:
            dt = Ttot / Npntsmin
        elif Ttot/timescale_to_resolve > Npntsmax:
            dt = Ttot / Npntsmax
        else:
            dt = timescale_to_resolve
            
        # trim outliers and median filter to avoid fitting deep transits
        g = abs(f[g1]-np.median(f[g1])) <= Nsig*np.std(f[g1])
        tbin, fbin, efbin = llnl.boxcar(bjd[g1][g], medfilt(f[g1][g],medkernel),
                                        ef[g1][g], dt=dt, include_edges=True,
                                        tfull=bjd[g1])
        self.tbin  = np.append(self.tbin, tbin)
        self.fbin  = np.append(self.fbin, fbin)
        self.efbin = np.append(self.efbin, efbin)

        # fit the GP model and save
        _,resultsGP, mubin, sigbin = _get_GP(thetaGPs[i], tbin, fbin, efbin)
        fintmu, fintsig = interp1d(tbin, mubin), interp1d(tbin, sigbin)
        mu[g1], sig[g1] = fintmu(bjd[g1]), fintsig(bjd[g1])
        fcorr[g1] = f[g1] - mu[g1] + 1 if mu[g1].sum() > 0 else f[g1] - mu[g1]
        ##ef[g1] = np.repeat(llnl.MAD1d(fcorr[g1]), fcorr[g1].size)

    # save
    self.bjd, self.f, self.ef = bjd, f, ef
    self.mu, self.sig, self.fcorr = mu, sig, fcorr
    self.resultsGP_detrend = thetaGPs



def estimate_box_transit_model(P, T0, Rs, t, f, ef):
    '''Estimate the transit depth and duration given P and T0. Return 
    the box transit model.'''
    phase = foldAt(t, P, T0)
    phase[phase>.5] -= 1
    intransit_approx = (phase*P <= 15./60/24) & (phase*P >= -15./60/24)
    depth = np.median(f[intransit_approx])
    duration = rvs.transit_width(P, Rs, Rs, 
				 rvs.m2Rearth(np.sqrt(depth)*rvs.Rsun2m(Rs)), 
				 0)
    model = llnl.box_transit_model((P,T0,depth,duration), t)
    return model


def is_good_star(hdr):
    if hdr['TESSMAG'] > 12:
	return False
    elif hdr['TEFF'] >= 4200:
	return False
    else:
	return True


def main(fname):
    # get fits file
    print 'Downloading fits file...\n'
    hdu = download_one_fits(fname)

    # get LC
    hdr, bjd, f, ef = get_lc(hdu)

    # only continue if it's a bright M-dwarf
    #if not is_good_star(hdr):  TEMP
#	raise ValueError('Not a star of interest.')

    # fit systematics with a GP
    print 'Fitting LC with GP alone...\n'
    fname_short = fname.replace('.fits','')
    self = Selfitivity(fname_short)
    samplerGP, samplesGP, resultsGP = do_mcmc_0(self, bjd, f, ef, fname_short)
    self.add_GPsamples(samplesGP)
    self.add_GPresults(resultsGP)
    ##save_fits(samplesGP, 'Results/%s/GP_samples'%fname_short)
    ##save_fits(resultsGP, 'Results/%s/GP_results'%fname_short)

    # search for transits in the corrected LC and get the transit parameters guesses
    print 'Searching for transit-like events...\n'
    params, EBparams = find_transits(self, bjd, f, ef, resultsGP[0], hdr, fname_short)
    self.params_guess = params
    self.EBparams_guess = EBparams

    # run full mcmc with GP+transit model
    # how to do the GP with the full LC??
    #if Ntransits > 0:
 	#print 'Fitting LC with GP + %i transit models'%Ntransits
    	#sampler, samples, results = do_mcmc_N(resultsGP[0], params, bjd, f, ef)
	#save_fits(samples, 'Results/%s/full_samples'%fname_short)


if __name__ == '__main__':
    fname = sys.argv[1]
    #main(fname)
