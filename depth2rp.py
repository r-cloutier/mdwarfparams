from imports import *
import rvs
from uncertainties import unumpy as unp


def depth2rp(P_days, depth, duration_days, Ms, Rs):
    '''Compute the planet radius from the transit depth and 
    duration using the analtyical treatment from Mandel & Agol 2002'''
    assert 0 < depth < 1

    # compute distance from centres at T0
    sma = rvs.semimajoraxis(P_days, Ms, 0)
    a_Rs = rvs.AU2m(sma) / rvs.Rsun2m(Rs)
    assert a_Rs > 1
    b = rvs.impactparam_T(P_days, Ms, Rs, duration_days)
    assert abs(b) <= 1
    inc = float(rvs.inclination(P_days, Ms, Rs, b))
    z = compute_distance_center(a_Rs, inc)

    # compute size ratio (p=rp/Rs)
    p_simple = unp.sqrt(depth)
    if z <= 1-p_simple:
        p = p_simple

    else:
        ps = np.logspace(-6,0,1000)
        depths = p2depth_grazing(ps, z)
        if (np.nanmax(depths) < z) or (np.nanmin(depths) > z):
            p = p_simple
        else:
            fint = interp1d(ps, depths)
            p = float(fint(depth))

    # compute planet radius
    rp = rvs.m2Rearth(rvs.Rsun2m(p*Rs))
    return rp
            

def compute_distance_center(a_Rs, inc_deg,
                            f=np.pi/2, ecc=0, omega=0):
    '''from batman paper
    http://iopscience.iop.org/article/10.1086/683602/pdf
    '''
    inc_rad = np.deg2rad(inc_deg)
    rt = np.sqrt(1. - np.sin(omega+f)**2 * np.sin(inc_rad)**2)
    z = a_Rs * (1-ecc**2) / (1. + ecc*np.cos(f)) * rt   # z = d/Rs
    return z


def p2depth_grazing(ps, z):
    kappa0 = np.arccos((ps*ps + z*z - 1.) / (2.*ps*z))
    kappa1 = np.arccos((1. - ps*ps + z*z) / (2.*z))
    bracket = ps*ps*kappa0 + kappa1 - np.sqrt(.25*(4*z*z - (1.+z*z-ps*ps)**2))
    depths = bracket / np.pi
    return depths
