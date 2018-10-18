from imports import *


def compute_perrakis_FML(sampler, ndim, nbins=100):
    '''Compute the perrakis estimate of the FML from an emcee sampler object.'''
    samples = sampler.chain.reshape((-1, ndim))
    nsamples, nparams = samples.shape
    
    # estimate posterior density for each parameter
    marginal_posterior_density = np.zeros(samples.shape)
    for i in range(nparams):
        marginal_posterior_density[:,i] = estimate_density(samples[:,i],
                                                           nbins=nbins)
    
    prod_marginal_densities = marginal_posterior_density.prod(axis=1)

    # compute lnprob from the sampler
    #lnPs, lnLs = np.zeros(nsamples), np.zeros(nsamples)
    #for i in range(nsamples):
    #    lnPs[i] = lnprior(samples[i])
    #    lnLs[i] = lnlike(samples[i], bjd, rv, erv)
    #lnprobs = lnLs + lnPs
    lnprobs = sampler.lnprobability.flatten()    
    assert lnprobs.size == nsamples
    
    # Use identity for summation
    # http://en.wikipedia.org/wiki/List_of_logarithmic_identities#Summation.2Fsubtraction
    # ln(sum(x)) = ln(x[0]) + ln(1 + sum( exp( ln(x[1:]) - ln(x[0]) ) ) )
    # log_summands = log_likelihood[cond] + np.log(prior_probability[cond])
    #  - np.log(prod_marginal_densities[cond])
    log_summands = (lnprobs - np.log(prod_marginal_densities))
    perr = log_sum(log_summands) - np.log(len(log_summands))
    
    return perr


def compute_perrakis_FMLv2(samples, lnprobs, nbins=100):
    '''Compute the perrakis estimate of the FML from an emcee sampler object.'''
    nsamples, nparams = samples.shape
    
    # estimate posterior density for each parameter
    marginal_posterior_density = np.zeros(samples.shape)
    for i in range(nparams):
        marginal_posterior_density[:,i] = estimate_density(samples[:,i],
                                                           nbins=nbins)
    
    prod_marginal_densities = marginal_posterior_density.prod(axis=1)

    assert lnprobs.size == nsamples
    
    # Use identity for summation
    # http://en.wikipedia.org/wiki/List_of_logarithmic_identities#Summation.2Fsubtraction
    # ln(sum(x)) = ln(x[0]) + ln(1 + sum( exp( ln(x[1:]) - ln(x[0]) ) ) )
    # log_summands = log_likelihood[cond] + np.log(prior_probability[cond])
    #  - np.log(prod_marginal_densities[cond])
    log_summands = (lnprobs - np.log(prod_marginal_densities))
    perr = log_sum(log_summands) - np.log(len(log_summands))
    
    return perr


def estimate_density(samples1D, nbins=100):
    density, bin_edges = np.histogram(samples1D, nbins, density=True)

    # Find to which bin each element corresponds
    density_indexes = np.searchsorted(bin_edges, samples1D, side='left')

    # Correct to avoid index zero from being assiged to last element
    density_indexes = np.where(density_indexes > 0, density_indexes,
                               density_indexes + 1)

    return density[density_indexes - 1]


def log_sum(log_summands):
    a = np.inf
    x = log_summands.copy()
    while a == np.inf or a == -np.inf or np.isnan(a):
        a = x[0] + np.log(1 + np.sum(np.exp(x[1:] - x[0])))
        np.random.shuffle(x)
    return a
