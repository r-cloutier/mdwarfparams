from imports import *
from scipy.interpolate import LinearNDInterpolator as lint
import mwdust


class GAIAMdwarfs():

    def __init__(self, fname):
	self.fname = fname


    def make_cuts(self):
 	self.good = (self.RsN<1) & (self.loggN>3) & (self.MsN<1) & \
                    (self.parallax_reliable==1)
        

    def _pickleoject(self, fname=''):
	fname = self.fname if fname == '' else fname
        fObj = open('%s'%fname, 'wb')
        pickle.dump(self, fObj)
        fObj.close()


# general functions
def plot_error_scatter(x, ux, lx, y, uy, ly, g=[], c=[], xlim=(), ylim=(),
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
    if len(c) > 0:
        plt.errorbar(x[g], y[g], xerr=[lx[g], ux[g]], yerr=[ly[g], uy[g]],
                     fmt='k.', elinewidth=.5, ms=0)
        plt.scatter(x[g], y[g], s=20, c=c[g], cmap=plt.get_cmap('rainbow'))
        plt.colorbar()
    else:
        plt.errorbar(x[g], y[g], xerr=[lx[g], ux[g]], yerr=[ly[g], uy[g]],
                     fmt='k.', elinewidth=.5)
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
    

def compute_GAIAmag(Hmag, eHmag, Kmag, eKmag):
    '''using https://www.aanda.org/articles/aa/pdf/2018/08/aa32756-18.pdf
    Table A2'''
    Hmag = unp.uarray(Hmag, eHmag)
    Kmag = unp.uarray(Kmag, eKmag)
    H_K = Hmag - Kmag
    p = np.poly1d([-1.359, 12.073, 0.6613])
    GAIAmag = p(H_K) + Kmag
    GAIAmag = unp.uarray(unp.nominal_values(GAIAmag), \
                         np.sqrt(unp.std_devs(GAIAmag)**2 + .3692**2))
    return unp.nominal_values(GAIAmag), unp.std_devs(GAIAmag)


def offset_parallax(parallax_mas):
    '''Apply a systematic offset to the GAIA parallaxes in milli-arcseconds 
    from https://arxiv.org/pdf/1804.09366.pdf'''
    offset = .03
    return parallax_mas + offset


def compute_distance(parallax_mas, eparallax_mas):
    '''return distance in pc'''
    par = unp.uarray(parallax_mas, eparallax_mas)
    dist_pc = rvs.m2pc(rvs.AU2m(1)) / unp.radians(par*1e-3/3600)
    return unp.nominal_values(dist_pc), unp.std_devs(dist_pc)
    

def compute_AK_myself(EBV, eEBV):
    '''compute extinction coefficient using extinction vector following 
    what is done in Fulton+Petigura 2018 sect 3.3'''
    raise ValueError('this is old shit')
    b = .063
    RK = unp.uarray(.224,.224*.3)
    AK = unp.uarray(EBV, eEBV) * (RK+b)
    return unp.nominal_values(AK), unp.std_devs(AK)


def compute_AK_mwdust(ls, bs, dist, edist):
    '''Using the EB-V map from 2014MNRAS.443.2907S and the extinction vector
    RK = 0.31 from Schlafly and Finkbeiner 2011 (ApJ 737, 103)'''
    dustmap = mwdust.Combined15(filter='2MASS Ks')
    dist_kpc, edist_kpc = dist*1e-3, edist*1e-3
    AK, eAK = np.zeros(ls.size), np.zeros(ls.size)
    for i in range(ls.size):
        v = dustmap(ls[i], bs[i],
                    np.array([dist_kpc[i], dist_kpc[i]+edist_kpc[i]]))
        AK[i], eAK[i] = v[0], abs(np.diff(v))
    return AK, eAK


def compute_MK(Kmag, eKmag, dist_pc, edist_pc, AK, eAK):
    Kmag = unp.uarray(Kmag, eKmag)
    dist = unp.uarray(dist_pc, edist_pc)
    mu = 5.*unp.log10(dist) - 5.
    AK = unp.uarray(AK, eAK)
    MK = Kmag - mu - AK
    return unp.nominal_values(MK), unp.std_devs(MK)


def interpolate_BCK(Teff, eTeff, Ms, eMs, Rs, eRs, EBV, eEBV):
    '''interpolate the bolometric correction in the K-band from the MIST
    models: http://adsabs.harvard.edu/abs/2016ApJ...823..102C'''
    raise ValueError('this is old shit')
    # create arrays
    Teff = unp.uarray(Teff, eTeff)
    Ms = unp.uarray(Ms, eMs)
    Rs = unp.uarray(Rs, eRs)
    logg = unp.log10(6.67e-11*rvs.Msun2kg(Ms)*1e2 / rvs.Rsun2m(Rs)**2)
    EBV = unp.uarray(EBV, eEBV)
    AV = EBV*3.1
    
    # get MIST models
    TeffMISTfull, loggMISTfull, AVMISTfull, BC_KMISTfull = \
    np.loadtxt('input_data/UBVRIplus/fehp000.UBVRIplus',usecols=(0,1,3,12)).T
    lintBC = lint(np.array([TeffMISTfull, loggMISTfull, AVMISTfull]).T,
                  BC_KMISTfull)

    # interpolate over 1 sigma values for each star
    BCK, eBCK = np.zeros(Teff.size), np.zeros(Teff.size)
    for i in range(Teff.size):
        BCKtmp = np.zeros((3,3,3))
        for j in range(-1,2):
            for k in range(-1,2):
                for l in range(-1,2):
                    Teff_val = unp.nominal_values(Teff[i]) + \
                               j*unp.std_devs(Teff[i])
                    logg_val = unp.nominal_values(logg[i]) + \
                               k*unp.std_devs(logg[i])
                    AV_val = unp.nominal_values(AV[i]) + \
                             l*unp.std_devs(AV[i])
                    BCKtmp[j+1,k+1,l+1] = float(lintBC(Teff_val, logg_val,
                                                       AV_val))

        # get BCK
        BCK[i], eBCK[i] = np.median(BCKtmp), MAD(BCKtmp)

    return BCK, eBCK


def compute_BCK_from_V_J_MDWARFS(Vmag, eVmag, Jmag, eJmag):
    '''use cubic polynomial and coefficients from table 3 in 
    http://iopscience.iop.org/article/10.1088/0004-637X/804/1/64/pdf'''
    color2BCK = np.poly1d([6.263e-3, -0.09655, 0.6084, 1.421])
    V_J = unp.uarray(Vmag, eVmag) - unp.uarray(Jmag, eJmag)
    BCK = color2BCK(V_J)
    return unp.nominal_values(BCK), unp.std_devs(BCK)
    

def MAD(arr):
    return np.median(abs(arr-np.median(arr)))

    
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
    b = (unp.nominal_values(MK) >= 5) & (unp.nominal_values(MK) < 10)
    d = (unp.nominal_values(MK) >= 4.5) & (unp.nominal_values(MK) < 5)
    Ms = unp.uarray(np.ones(MK.size) + np.nan, np.ones(MK.size) + np.nan)
    Ms[b] = Ms_B16[b]
    Ms[d] = Ms_D00[d]
    #Ms = Ms_D00
    return unp.nominal_values(Ms), unp.std_devs(Ms)


def compute_Ms_from_MK_M15(MK, eMK):
    g = (MK >= 4.5) & (MK < 9.5)
    MK = unp.uarray(MK, eMK)
    p = np.poly1d([-2.7262e-4, .0106, -.1217, .3872, .5858])  
    Ms = unp.uarray(np.ones(MK.size) + np.nan, np.ones(MK.size) + np.nan)
    Ms[g] = p(MK[g])
    return unp.nominal_values(Ms), unp.std_devs(Ms)


def compute_Rs_from_Ms(Ms, eMs):
    '''from Eq 10 in Boyajian+2012'''
    a = unp.uarray(.32, .0165)
    b = unp.uarray(.6063, .0153)
    c = unp.uarray(.0906, .0027)
    Ms = unp.uarray(Ms, eMs)
    Rs = a*Ms*Ms + b*Ms + c
    return unp.nominal_values(Rs), unp.std_devs(Rs)


def compute_Mbol(BCK, eBCK, MK, eMK):
    BCK = unp.uarray(BCK, eBCK)
    MK = unp.uarray(MK, eMK)
    Mbol = BCK + MK
    return unp.nominal_values(Mbol), unp.std_devs(Mbol)


def compute_Rs_from_Mbol(Mbol, eMbol, Teff, eTeff):
    L0, sigma = 3.0128e28, 5.670367e-8
    Mbol = unp.uarray(Mbol, eMbol)
    Lbol = L0*10**(-.4*Mbol)
    Teff = unp.uarray(Teff, eTeff)
    Rs = rvs.m2Rsun(unp.sqrt(Lbol / (4*np.pi*sigma*Teff**4)))
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
