import numpy as np
from scipy.stats import gaussian_kde
from scipy.ndimage import gaussian_filter1d


def lnjeffreysprior(p, pmin, pmax):
    '''See Gregory 2005, Eq 17'''
    pri = 1./(p*np.log(pmax/pmin))
    if np.isfinite(pri):
        return np.log(pri)
    else:
        return -np.inf

def lnmodjeffreysprior(p, pknee, pmax):
    '''See Gregory 2005, Eq 16'''
    pri = 1./(p+pknee) * 1./(np.log((pknee+pmax)/pknee))
    if np.isfinite(pri):
        return np.log(pri)
    else:
        return -np.inf

def lnuniform(p, pmin, pmax):
    '''See Gregory 2005, Eq 17'''
    if pmax < p or pmin > p:
        return -np.inf
    pri = 1./(pmax-pmin)
    if np.isfinite(pri):
        return np.log(pri)
    else:
        return -np.inf

def lngaussian(x, mu, sig):
    return -.5*((x-mu)/sig)**2


def get_results(samples, sigma=5):
    MAPs = np.zeros(samples.shape[1])
    plus_1sigs = np.zeros(samples.shape[1])
    min_1sigs = np.zeros(samples.shape[1])
    for i in range(MAPs.size):
	samp = samples[:,i]
	samp = samp[np.isfinite(samp)]
	try:
            kernel = gaussian_kde(samp)
            xarr = np.linspace(samp.min(), samp.max(), 1000)
            probs = kernel.pdf(xarr) / kernel.pdf(xarr).sum()
            probs = gaussian_filter1d(probs, sigma)
            MAPs[i] = xarr[probs==probs.max()][0]
	    # get percentiles
	    xarr = np.random.choice(xarr, 1000, p=probs/probs.sum())
	    v = np.percentile(xarr, (16,84))
	    plus_1sigs[i], min_1sigs[i] = v[1]-MAPs[i], MAPs[i]-v[0]
	except ValueError:
	    MAPs[i], plus_1sigs[i], min_1sigs[i] = np.repeat(np.nan,3)

    return np.array([MAPs, plus_1sigs, min_1sigs])
