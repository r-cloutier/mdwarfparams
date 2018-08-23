from imports import *

class GAIAMdwarfs():

    def __init__(self, fname):
	self.fname = fname


    def _pickleoject(self, fname=''):
	fname = self.fname if fname == '' else fname
        fObj = open('%s'%fname, 'wb')
        pickle.dump(self, fObj)
        fObj.close()


# general functions
def plot_error_scatter(x, ux, lx, y, uy, ly, g=[], xlim=(), ylim=(),
                       xlog=False, ylog=False):
    assert x.size == ux.size
    assert x.size == lx.size
    assert x.size == y.size
    assert x.size == uy.size
    assert x.size == ly.size
    if len(g) > 0:
        assert len(g) == x.size
        g = (np.isfinite(ux)) & (np.isfinite(lx)) & \
            (np.isfinite(uy)) & (np.isfinite(ly)) & g
    else:
        g = (np.isfinite(ux)) & (np.isfinite(lx)) & \
            (np.isfinite(uy)) & (np.isfinite(ly))
    plt.errorbar(x[g], y[g], xerr=[lx[g], ux[g]], yerr=[ly[g], uy[g]], fmt='k.')
    plt.plot([x[g].min(),x[g].max()], [x[g].min(),x[g].max()], 'b--')
    if len(xlim) == 2:
        plt.xlim(xlim)
    if len(ylim) == 2:
        plt.ylim(ylim)
    if xlog:
        plt.xscale('log')
    if ylog:
        plt.yscale('log')
    plt.show()
    

def compute_distance(parallax_mas, eparallax_mas):
    '''return distance in pc'''
    par = unp.uarray(parallax_mas, eparallax_mas)
    dist_pc = rvs.m2pc(rvs.AU2m(1)) / unp.radians(par*1e-3/3600)
    return unp.nominal_values(dist_pc), unp.std_devs(dist_pc)
    

def compute_MK(Kmag, eKmag, parallax, eparallax):
    dist_pc, edist_pc = compute_distance(parallax, eparallax)
    dist = unp.uarray(dist_pc, edist_pc)
    Kmag = unp.uarray(Kmag, eKmag)
    MK = Kmag - 5.*unp.log10(dist) + 5.
    return unp.nominal_values(MK), unp.std_devs(MK)
    

def compute_Ms_from_MK(MK, eMK):
    '''from Eq 11 and Table 13 in Benedict+2016 and 
    Eq on page 220 of Delfosse+2000'''
    # compute all masses using B16
    C0 = unp.uarray(.2311, .0004)
    C1 = unp.uarray(-.1352, .0007)
    C2 = unp.uarray(.04, .0005)
    C3 = unp.uarray(.0038, .0002)
    C4 = unp.uarray(-.0032, .0001)
    x0 = 7.5
    MK = unp.uarray(MK, eMK)
    Ms_B16 = C0 + C1*(MK-x0) + C2*(MK-x0)**2 + C3*(MK-x0)**3 + C4*(MK-x0)**4

    # compute all masses using D00
    coeffs = 1e-3 * np.array([0.37529, -6.2315, 13.205, 6.12, 1.8])
    p = np.poly1d(coeffs)
    Ms_D00 = 10**(p(MK))

    # select masses over valid ranges
    b = unp.nominal_values(MK) >= 5
    d = (unp.nominal_values(MK) >= 4.5) & (unp.nominal_values(MK) < 5)
    g = (unp.nominal_values(MK) < 4.5)
    assert b.sum() + d.sum() + g.sum() == MK.size
    Ms = unp.uarray(np.ones(MK.size) + np.nan, np.ones(MK.size) + np.nan)
    Ms[b] = Ms_B16[b]
    Ms[d] = Ms_D00[d]

    return unp.nominal_values(Ms), unp.std_devs(Ms)
    

def compute_Rs_from_Ms(Ms, eMs):
    '''from Eq 10 in Boyajian+2012'''
    a = unp.uarray(.32, .0165)
    b = unp.uarray(.6063, .0153)
    c = unp.uarray(.0906, .0027)
    Ms = unp.uarray(Ms, eMs)
    Rs = a*Ms*Ms + b*Ms + c
    return unp.nominal_values(Rs), unp.std_devs(Rs)


def compute_sma_from_Ms(Ps, ePs, Ms, eMs):
    Ps = unp.uarray(Ps, ePs)
    Ms = unp.uarray(Ms, eMs)
    sma = rvs.semimajoraxis(Ps, Ms, 0)
    return unp.nominal_values(sma), unp.std_devs(sma)


def compute_F_Rs_Ms(Ps, ePs, Ms, eMs, Rs, eRs, Teff, eTeff):
    sma, esma = compute_sma_from_Ms(Ps, ePs, Ms, eMs)
    sigma = 5.670367e-8
    sma = unp.uarray(sma, esma)
    Rs = unp.uarray(Rs, eRs)
    Teff = unp.uarray(Teff, eTeff)
    F = sigma * (rvs.Rsun2m(Rs)/rvs.AU2m(sma))**2 * Teff**4 / 1367.
    return unp.nominal_values(F), unp.std_devs(F)


def compute_mp_from_K(P, eP, Ms, eMs, K, eK):
    '''for use with the new GAIA-derived Ms.'''
    P = unp.uarray(P, eP)
    Ms = unp.uarray(Ms, eMs)
    K = unp.uarray(K, eK)
    mp = rvs.RV_mp(P, Ms, K)
    return unp.nominal_values(mp), unp.std_devs(mp)


def compute_rp_from_rpRS(rpRs, erpRs, Rs, eRs):
    '''for use with the new GAIA-derived Rs.'''
    rpRs = unp.uarray(rpRs, erpRs)
    Rs = unp.uarray(Rs, eRs)
    rp = rvs.m2Rearth(rvs.Rsun2m(rpRs * Rs))
    return unp.nominal_values(rp), unp.std_devs(rp)

    
def compute_planet_density(mp, emp, rp, erp):
    '''units: Mearth and Rearth'''
    mp = unp.uarray(mp, emp)
    rp = unp.uarray(rp, erp)
    rho = 5.51 * mps / rps**3  # g/cm3
    return unp.nominal_values(rho), unp.std_devs(rho)
    
        
def loadpickle(fname):
    fObj = open(fname, 'rb')
    self = pickle.load(fObj)
    fObj.close()
    return self
