from imports import *
from priors import *
from linear_lnlike import *


def get_modelN(theta, t, f, ef, u1, u2):
    a, l, G, P = np.exp(theta[:4])
    s = theta[4]
    k1 = george.kernels.ExpSquaredKernel(l)
    k2 = george.kernels.ExpSine2Kernel(G,P)
    gp = george.GP(a*(k1+k2))
    try:
        gp.compute(t, np.sqrt(ef**2 + s**2))
    except (ValueError, np.linalg.LinAlgError):
        return -np.inf
    Ntransits = (theta.size-5) / 5
    assert Ntransits >= 1
    fmodel = np.zeros(t.size)
    for i in range(Ntransits):
	P, T0, aRs, rpRs, inc = theta[5+5*i:10+5*i]
        func = transit_model_func_curve_fit(u1, u2)
        fmodel += func(t, P, T0, aRs, rpRs, inc) - 1
    mu, cov = gp.predict(f-fmodel+1, t)
    sig = np.sqrt(np.diag(cov))
    return gp, mu, sig


def lnlike(theta, t, f, ef, u1, u2):
    gp,_,_ = get_modelN(theta, t, f, ef, u1, u2)
    lnL = gp.lnlikelihood(f, quiet=True)
    return lnL


def lnprior(theta, theta0, inclims):
    Ntransits = (theta.size-5) / 5
    P0, T00, aRs0, rpRs0, inc0 = theta0
    lna, lnl, lnG, lnP, s = theta[:5]
    lps = np.zeros(5+5*Ntransits)
    lps[0] = lnuniform(lna, np.log(1e-5), np.log(1))
    lps[1] = lnuniform(lnl, np.log(1), np.log(1e4))
    lps[2] = lnuniform(lnG, np.log(1e-2), np.log(1e2))
    lps[3] = lnuniform(lnP, np.log(1), np.log(5e2))
    lps[4] = lnuniform(s, 0, 1) if s > 0 else -np.inf
    for i in range(Ntransits):
        P, T0, aRs, rpRs, inc = theta[5+5*i:10+5*i]
	lps[5*i+5] = lnuniform(P, P0*.9, P0*1.1)
        lps[5*i+6] = lnuniform(T0, T00-P0*1.1, T00+P0*1.1)
        lps[5*i+7] = lnuniform(aRs, aRs0*.9, aRs0*1.1)
        lps[5*i+8] = lnuniform(rpRs, 0, 1)
        lps[5*i+9] = lnuniform(inc, inclims[0], inclims[1])
    return lps.sum()


def lnprob(theta, theta0, t, f, ef, inclims, u1, u2):
    lp = lnprior(theta, theta0, inclims)
    if np.isfinite(lp):
        return lp + lnlike(theta, t, f, ef, u1, u2)
    else:
        return -np.inf


def run_emcee(theta, t, f, ef, initialize, u1, u2, Ms, Rs,
              nwalkers=100, burnin=200, nsteps=400, a=2):
    '''Run mcmc on an input light curve with no transit model.'''
    # initialize chains
    assert len(theta) == len(initialize)
    assert len(theta) == 10
    p0 = []
    for i in range(nwalkers):
    	p0.append(theta + initialize*np.random.randn(len(theta)))
    
    # initialize sampler
    P = theta[5]
    inclims = np.array([float(rvs.inclination(P,Ms,Rs,1)),
                        float(rvs.inclination(P,Ms,Rs,-1))])
    args = (theta[5:10], t, f, ef, inclims, u1, u2)
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
    print 'Running MCMC (RVs)...'
    p0,_,_ = sampler.run_mcmc(p0, nsteps)
    samples = sampler.chain.reshape((-1, len(theta)))
    print "Mean acceptance fraction: %.4f"%np.mean(sampler.acceptance_fraction)
    print 'Full MCMC took %.4f minutes'%((time.time()-t0)/60.)

    return sampler, samples, get_results(samples)
