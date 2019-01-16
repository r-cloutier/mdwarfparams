from occurrencerateclass import *
import linear_lnlike as llnl
import rvs
from priors import *


def optimize_singletransit_params(params, bjd, fcorr, ef, Ms, Rs, u1, u2,
                                  pltt=True):
    '''Get best-fit parameters using the periodic fit parameters for 
    initialization (i.e. P_input < P_singletransit).'''
    assert params.shape == (4,)
    P, T0, depth, duration = params
    if depth >= .9:  # sometimes the dimming is passed instead of depth
        return np.nan, np.nan, np.nan, np.nan, \
            np.repeat(np.nan, bjd.size), np.repeat(np.nan, 7)
    # focus on data centered around T0
    g = (bjd >= T0-10*duration) & (bjd <= T0+10*duration)
    bjdred, fcorrred, efred = bjd[g], fcorr[g], ef[g]
    
    # initialize
    aRs = rvs.AU2m(rvs.semimajoraxis(P,Ms,0)) / rvs.Rsun2m(Rs)
    rpRs = np.sqrt(depth)
    p0 = P, T0, aRs, rpRs, 90.
    incs = np.array([float(rvs.inclination(P,Ms,Rs,1)),
                     float(rvs.inclination(P,Ms,Rs,-1))])
    bnds = ((0, T0-.2, 0, 0, incs.min()),
            (P*100, T0+.2, aRs*100, 1, incs.max()))
    try:
        popt,_ = curve_fit(llnl.transit_model_func_curve_fit(u1,u2),
                           bjdred, fcorrred, p0=p0, sigma=efred,
                           absolute_sigma=True, bounds=bnds)
        P, T0, aRs, rpRs, inc = popt
        depth = rpRs**2
        b = rvs.impactparam_inc(P, Ms, Rs, inc)
        duration = P/(np.pi*aRs) * np.sqrt((1+np.sqrt(depth))**2 - b*b)
        func = llnl.transit_model_func_curve_fit(u1, u2)
        fmodel = func(bjdred, P, T0, aRs, rpRs, inc)
        params = np.array([P,T0,aRs,rpRs,inc,u1,u2])
    except RuntimeError, ValueError:
        func = llnl.transit_model_func_curve_fit(u1, u2)
        fmodel = func(bjdred, P, T0, aRs, rpRs, 90.)
        P, T0, depth, duration = np.repeat(np.nan, 4)
        params = np.repeat(np.nan, 7)

    # plot
    if pltt:
        plt.plot(bjdred, fcorrred, '.', bjdred, fmodel, '-')
        plt.show()        
        
    return bjdred, fcorrred, efred, fmodel, params



def run_mcmc(params, bjd, fcorr, ef, Ms, eMs, Rs, eRs, Teff, eTeff, u1, u2,
             nwalkers=200, burnin=200, nsteps=400, pltt=True):
    assert params.shape == (4,)
    params_optimized = np.zeros(5)
    nwalkers, burnin, nsteps = int(nwalkers), int(burnin), int(nsteps)
    params_results = np.zeros((3,5))

    # get optimized parameters and get the LC around T0
    bjdred, fcorrred, efred, fmodel_opt, theta = \
                            optimize_singletransit_params(params, bjd, fcorr,
                                                          ef, Ms, Rs, u1, u2,
                                                          pltt=1)
    params_opt = theta[:5]

    # run MCMC on transit LC
    initialize = [1, 1e-3, 1, 1e-2, 1e-1]
    print 'Running MCMC on single-transit model'
    sampler,samples = run_emcee(params_opt, bjd, bjdred, fcorrred, efred,
                                initialize, u1, u2, Ms, Rs, a=2,
                                nwalkers=nwalkers, burnin=burnin, nsteps=nsteps,
                                zeroplanetmodel=False)
    results = get_results(samples)
    params_results = results
    func = llnl.transit_model_func_curve_fit(u1, u2)
    fmodel_mcmc = func(bjdred, *params_results[0])

    # estimate single transit P, aRs, rpRs, and inc
    Nsamp = samples.shape[0] 
    samp_Ms = np.random.randn(Nsamp)*eMs + Ms
    samp_Rs = np.random.randn(Nsamp)*eRs + Rs
    samp_rho = rvs.Msun2kg(samp_Ms) / rvs.Rsun2m(samp_Rs)**3
    samp_aRs = samples[:,2]
    v = np.percentile(samp_aRs, (16,50,84))
    aRs_est = [v[1], v[2]-v[1], v[1]-v[0]]
    samp_rpRs = samples[:,3]
    v = np.percentile(samp_rpRs, (16,50,84))
    rpRs_est = [v[1], v[2]-v[1], v[1]-v[0]]
    samp_inc = samples[:,4]
    v = np.percentile(samp_inc, (16,50,84))
    inc_est = [v[1], v[2]-v[1], v[1]-v[0]]
    #samp_P=rvs.sec2days(np.sqrt(4*np.pi*np.pi/(6.67e-11*samp_rho)*samp_aRs**3))
    samp_P = samples[:,0]
    v = np.percentile(samp_P, (16,50,84))
    P_est = [v[1], v[2]-v[1], v[1]-v[0]]
   
    # estimate F
    samp_Teff = np.random.randn(Nsamp)*eTeff + Teff
    samp_Ls = compute_Ls(samp_Rs, samp_Teff)
    samp_F = compute_F(samp_Ls, rvs.semimajoraxis(samp_P, samp_Ms, 0))
    v = np.percentile(samp_F, (16,50,84))
    F_est = [v[1], v[2]-v[1], v[1]-v[0]]
 
    # plotting
    if pltt:
        t0 = 2457000
        plt.figure(1)
        plt.plot(bjdred-t0, fcorrred, '.')
        plt.plot(bjdred-t0, fmodel_opt, '-', lw=3, label='optimized')
        plt.plot(bjdred-t0, fmodel_mcmc, '-', lw=2, label='MCMC model')
        plt.legend(loc='upper left')
        plt.figure(2)
        plt.hist(samp_P, bins=30)
        plt.show()  

    return bjd, fcorr, ef, params_opt, fmodel_opt, samples, params_results, fmodel_mcmc, samp_P, P_est, samp_F, F_est, aRs_est, rpRs_est, inc_est




def get_model1(theta, bjd, u1, u2):
    assert theta.size == 5
    P, T0, aRs, rpRs, inc = theta
    func = llnl.transit_model_func_curve_fit(u1, u2)
    fmodel = func(bjd, P, T0, aRs, rpRs, inc)
    return fmodel


def lnlike(theta, bjd, f, ef, u1, u2):
    fmodel = get_model1(theta, bjd, u1, u2)
    lnL = -.5*(np.sum((f-fmodel)**2/ef**2 - np.log(1./ef**2)))
    return lnL


def lnprior(theta, theta0, Plim, aRslim, inclims):
    P0, T00, aRs0, rpRs0,_ = theta0
    P, T0, aRs, rpRs, inc = theta
    lps = np.zeros(5)
    lps[0] = lnjeffreysprior(P, Plim, P0*100) if P > Plim else -np.inf
    lps[1] = lnuniform(T0, T00-.3, T00+.3)
    lps[2] = lnjeffreysprior(aRs, aRslim, aRs0*100)
    lps[3] = lnuniform(rpRs, 0, 1)
    lps[4] = lnuniform(inc, inclims[0], inclims[1])
    return lps.sum()


def lnprob(theta, theta0, bjd, f, ef, Plim, aRslim, inclims, u1, u2,
           zeroplanetmodel):
    if zeroplanetmodel:  # set rp/Rs to zero if considering a zero-planet model
	theta[3] = 0.
    lp = lnprior(theta, theta0, Plim, aRslim, inclims)
    if np.isfinite(lp):
        return lp + lnlike(theta, bjd, f, ef, u1, u2)
    else:
        return -np.inf


def run_emcee(theta, bjd, bjdred, fred, efred, initialize, u1, u2, Ms, Rs,
              nwalkers=200, burnin=200, nsteps=400, a=2, 
	      zeroplanetmodel=False):
    '''Run mcmc on an input light curve with no transit model.'''
    # get limits on P and aRs
    P, T0 = theta[:2]
    Plim = np.max([T0-bjd.min(), bjd.max()-T0])
    aRslim = rvs.AU2m(rvs.semimajoraxis(Plim,Ms,0)) / rvs.Rsun2m(Rs)
    print 'Plim = %.4f'%Plim
    print 'aRslim = %.4f'%aRslim
    
    # initialize chains
    assert len(theta) == len(initialize)
    assert len(theta) == 5
    theta[0], theta[2] = Plim+10, aRslim+10
    p0 = []
    for i in range(nwalkers):
    	p0.append(theta + initialize*np.random.randn(len(theta)))

    
    # initialize sampler
    P = theta[0]
    inclims = np.array([float(rvs.inclination(P,Ms,Rs,1)),
                        float(rvs.inclination(P,Ms,Rs,-1))])
    args = (theta, bjdred, fred, efred, Plim, aRslim, inclims, u1, u2,
            zeroplanetmodel)
    sampler = emcee.EnsembleSampler(nwalkers, len(theta), lnprob, args=args,
                                    a=a)

    # run burnin
    print 'Running burnin...'
    t0 = time.time()
    p0,_,_ = sampler.run_mcmc(p0, burnin)
    print 'Burnin acceptance fraction is %.4f'%np.mean(sampler.acceptance_fraction)
    print 'Burnin took %.4f minutes\n'%((time.time()-t0)/60.)
    sampler.reset()

    # run MCMC
    print 'Running full MCMC...'
    p0,_,_ = sampler.run_mcmc(p0, nsteps)
    samples = sampler.chain.reshape((-1, len(theta)))
    print "Mean acceptance fraction: %.4f"%np.mean(sampler.acceptance_fraction)
    print 'Full MCMC took %.4f minutes'%((time.time()-t0)/60.)

    return sampler, samples
