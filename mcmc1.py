from imports import *
from priors import *
from linear_lnlike import *


def get_model1_OLD(theta, bjd, f, ef, u1, u2):
    assert theta.size == 5
    gp = george.GP(0*george.kernels.ExpSquaredKernel(1))
    try:
        gp.compute(bjd)
    except (ValueError, np.linalg.LinAlgError):
        return -np.inf
    P, T0, aRs, rpRs, inc = theta
    func = transit_model_func_curve_fit(u1, u2)
    fmodel = func(bjd, P, T0, aRs, rpRs, inc)
    return gp, fmodel


def get_model1(theta, bjd, u1, u2):
    assert theta.size == 5
    P, T0, aRs, rpRs, inc = theta
    func = transit_model_func_curve_fit(u1, u2)
    fmodel = func(bjd, P, T0, aRs, rpRs, inc)
    return fmodel


def lnlike(theta, bjd, f, ef, u1, u2):
    #gp, fmodel = get_model1(theta, bjd, f, ef, u1, u2)
    #lnL = gp.lnlikelihood(f-fmodel, quiet=True)
    fmodel = get_model1(theta, bjd, u1, u2)
    lnL = -.5*(np.sum((f-fmodel)**2/ef**2 - np.log(1./ef**2)))
    return lnL


def lnprior(theta, theta0, inclims):
    P0, T00, aRs0, rpRs0,_ = theta0
    P, T0, aRs, rpRs, inc = theta
    lps = np.zeros(5)
    lps[0] = lnuniform(P, P0*.9, P0*1.1)
    lps[1] = lnuniform(T0, T00-P0*1.1, T00+P0*1.1)
    lps[2] = lnuniform(aRs, aRs0*.7, aRs0*1.3)
    lps[3] = lnuniform(rpRs, 0, 1)
    lps[4] = lnuniform(inc, inclims[0], inclims[1])
    return lps.sum()


def lnprob(theta, theta0, bjd, f, ef, inclims, u1, u2, zeroplanetmodel):
    if zeroplanetmodel:  # set rp/Rs to zero if considering a zero-planet model
	theta[3] = 0.
    lp = lnprior(theta, theta0, inclims)
    if np.isfinite(lp):
        return lp + lnlike(theta, bjd, f, ef, u1, u2)
    else:
        return -np.inf


def run_emcee(theta, bjd, f, ef, initialize, u1, u2, Ms, Rs,
              nwalkers=100, burnin=200, nsteps=400, a=2, 
	      zeroplanetmodel=False):
    '''Run mcmc on an input light curve with no transit model.'''
    # initialize chains
    assert len(theta) == len(initialize)
    assert len(theta) == 5
    p0 = []
    for i in range(nwalkers):
    	p0.append(theta + initialize*np.random.randn(len(theta)))
    
    # initialize sampler
    P = theta[0]
    inclims = np.array([float(rvs.inclination(P,Ms,Rs,1)),
                        float(rvs.inclination(P,Ms,Rs,-1))])
    args = (theta, bjd, f, ef, inclims, u1, u2, zeroplanetmodel)
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
