from imports import *
from scipy.optimize import curve_fit
import vetting as vett
import rvs, batman
from scipy.interpolate import LinearNDInterpolator as lint

global dispersion_sig, depth_sig, bimodalfrac, T0tolerance, transitlikefrac, min_autocorr_coeff
#dispersion_sig, depth_sig, bimodalfrac = 3., 3., .5
#dispersion_sig, depth_sig, bimodalfrac = 2., 1.35, .5  # v3
#dispersion_sig, depth_sig, bimodalfrac = 1.6, 1., .5
# for real K2 LCs
#dispersion_sig, depth_sig, bimodalfrac, T0tolerance, transitlikefrac = \
                                                        #2., 1., .6, .1, .7
dispersion_sig, depth_sig, bimodalfrac, T0tolerance, transitlikefrac = \
							2.8, 8., .7, .1, .7
min_autocorr_coeff = .6


def lnlike(bjd, f, ef, fmodel):
    return -.5*(np.sum((f-fmodel)**2 / ef**2 - np.log(1./ef**2)))


def box_transit_model(theta, t):
    '''Return the box transmit model.'''
    fmodel = np.ones(t.size)
    P, T0, depth, duration = theta
    phase = foldAt(t, P, T0)
    phase[phase>.5] -= 1
    intransit = (phase*P <= .5*duration) & (phase*P >= -.5*duration)
    fmodel[intransit] = 1. - depth
    return fmodel


def box_transit_model_curve_fit(t, P, T0, depth, duration):
    '''Return the box transmit model.'''
    fmodel = np.ones(t.size)
    phase = foldAt(t, P, T0)
    phase[phase>.5] -= 1
    intransit = (phase*P <= .5*duration) & (phase*P >= -.5*duration)
    fmodel[intransit] = 1. - depth
    return fmodel


def box_transit_model_time(theta, t):
    '''Return the box transit model as a function of transit time 
    instead of P and T0.'''
    fmodel = np.ones(t.size)
    T, depth, duration = theta 
    intransit = (t >= T-.5*duration) & (t <= T+.5*duration)
    fmodel[intransit] = 1. - depth
    return fmodel


def get_depth_lnlike(theta, bjd, fcorr, ef, depthmax=.1, N=2e2):
    '''Given a transit time and duration, get the max lnL depth and its lnL.'''
    # sample depths
    depthmin, N = np.median(ef), int(N)
    depths = 10**(np.random.uniform(np.log10(depthmin), np.log10(depthmax), N))

    # compute lnLs of the depths given a transit time and duration
    T, duration = theta
    lnLs = np.zeros(N)
    for i in range(N):
    	fmodel = box_transit_model_time((T,depths[i],duration), bjd)
        lnLs[i] = lnlike(bjd, fcorr, ef, fmodel)

    # get maximum-likelihood depth
    g = lnLs == lnLs.max()
    if g.sum() == 1:
	lnL, D = float(lnLs[g]), float(depths[g])
    elif g.sum() == 0:
	lnL, D = np.nan, np.nan
    else:
	g = np.where(g)[0][0]
	lnL, D = lnLs[g], depths[g]
    return lnLs, depths, lnL, D


def linear_search(bjd, fcorr, ef):
    '''Evaluate the lnL as a function of transit duration and 
    transit time/epoch of a transit.'''
    # setup transit time and duration grids
    transit_times = np.arange(bjd.min(), bjd.max(), 30./60/24)
    durations = np.array([1.2,2.4,4.8]) / 24  # coarse grid in transit duration

    # get max lnL depth over the linear grid of transit times and durations
    NT, ND = transit_times.size, durations.size
    lnLs, depths = np.zeros((NT,ND)), np.zeros((NT,ND))
    for i in range(NT):
	print float(i) / NT
	for j in range(ND):
	    theta = transit_times[i], durations[j]
	    _,_,lnLs[i,j],depths[i,j] = get_depth_lnlike(theta, bjd, fcorr, ef)

    return transit_times, durations, lnLs, depths


def MAD1d(arr):
    assert len(arr.shape) == 1
    return np.median(abs(arr-np.median(arr)))

def MAD2d(arr):
    assert len(arr.shape) == 2
    mads = np.zeros(arr.shape[1])
    for i in range(mads.size):
	mads[i] = np.median(abs(arr[:,i] - np.median(arr[:,i])))
    return mads


def find_transit_parameters(bjd, fcorr, ef, 
			    transit_times, durations, lnLs, depths, 
			    SNRthresh):
    '''Find periodic transits and return initial guesses of P, T0, 
    duration, and depth.'''
    # find high S/N peaks in lnL as a function of the transit time
    SNRs = (lnLs-np.median(lnLs, axis=0)) / MAD2d(lnLs)
    gSNR = SNRs > SNRthresh

    # search over coarse duration grid
    NT, ND = transit_times.size, durations.size
    Ps_full, T0s_full, durations_full, depths_full, lnLs_full = np.zeros(0), \
						     		np.zeros(0), \
                                                     		np.zeros(0), \
                                                     		np.zeros(0), \
								np.zeros(0)
    for i in range(ND):

    	if gSNR[:,i].sum() == 0:  # no transit-like events
	    pass

   	else:  # some potential transit-like events
	    
	    # get unique approximate transit times of transit-like events
	    Ts_all = transit_times[gSNR[:,i]]
	    Ts_reduced = np.array([Ts_all[0]])
	    for j in range(1,Ts_all.size):
		if np.isclose(Ts_all[j]-Ts_reduced[-1], 0, atol=durations[i]*2):  # same transit
		    pass
		else:
		    Ts_reduced = np.append(Ts_reduced, Ts_all[j])

 	    # adjust to more accurate transit times
	    Ntransits = Ts_reduced.size
	    T0s = np.zeros(0)
	    for j in range(Ntransits):
		g = (bjd >= Ts_reduced[j]-2*durations[i]) & \
		    (bjd <= Ts_reduced[j]+2*durations[i])
		if g.sum() > 0:
		    fsmooth = gaussian_filter1d(fcorr[g], 5)
		    T0s = np.append(T0s, bjd[g][fsmooth == fsmooth.min()])
		    ##plt.plot(bjd, fcorr, '-', bjd[g], fcorr[g], 'o')
		    ##plt.plot(bjd[g], fsmooth, '-', lw=2)
		    ##plt.axvline(Ts_reduced[j]), plt.axvline(T0s[j], ls='--')
		    ##plt.show()
	    T0s = np.unique(T0s)

	    # search for all periodic transits (periods between transit events)
	    Ntransits = T0s.size
	    for j in range(Ntransits):
		for k in range(Ntransits):
		    if j != k:
			Ps_full = np.append(Ps_full, T0s[j]-T0s[k])
		      	T0s_full = np.append(T0s_full, T0s[j])
			durations_full = np.append(durations_full, durations[i])
		    	phase = foldAt(bjd, Ps_full[-1], T0s_full[-1])
		    	phase[phase > .5] -= 1
		    	intransit = (phase*Ps_full[-1] <= .5*durations_full[-1]) & \
				    (phase*Ps_full[-1] >= -.5*durations_full[-1])
                        depths_full = np.append(depths_full, 1-np.median(fcorr[intransit]))  # could be negative
			theta = Ps_full[-1], T0s_full[-1], depths_full[-1], durations_full[-1]
			fmodel = box_transit_model(theta, bjd)
			lnLs_full = np.append(lnLs_full, lnlike(bjd, fcorr, ef, fmodel))
			##plt.plot(phase, fcorr, '-'), plt.plot(phase[intransit], fcorr[intransit], 'o'), plt.show()

    # trim (must be greater than ~.5 and smaller than the baseline)
    g = (Ps_full >= .49) & (Ps_full <= (bjd.max()-bjd.min()))

    return Ps_full[g], T0s_full[g], durations_full[g], depths_full[g], lnLs_full[g]


def compute_P_lnL(bjd, fcorr, ef, theta, N=1e2):
    '''Compute lnL for each P marginalized over T0 then return the maxmimum 
    likelihood result.'''
    N = int(N)
    P, T00, Z, D = theta
    lnLs, thetas = np.zeros(N), np.zeros((N,4))
    for i in range(N):
	T0 = np.random.uniform(T00-.5*P,T00+.5*P)
	phase = foldAt(bjd, P, T0)
	phase[phase > .5] -= 1
	intransit = (phase*P <= .5*D) & (phase*P >= -.5*D)
	thetas[i] = P, T0, 1-np.median(fcorr[intransit]), D
	fmodel = box_transit_model(thetas[i], bjd)
	lnLs[i] = lnlike(bjd, fcorr, ef, fmodel)
    # get maximum lnL result
    g = lnLs == lnLs.max()
    return lnLs[g][0], thetas[g][0]


def compute_transit_lnL(bjd, fcorr, ef, transit_times, durations, lnLs,
                        depths, SNRthresh):
    '''Get transit parameters and compute the lnL to identify transits.'''
    # get transit parameters
    assert lnLs.shape == (transit_times.size, durations.size)
    assert depths.shape == (transit_times.size, durations.size)
    Ps, T0s, Ds, Zs, lnLs_transit = find_transit_parameters(bjd, fcorr, ef,
							    transit_times,
                                                            durations, 
							    lnLs, depths,
                                                            SNRthresh)
    assert Ps.size == T0s.size
    assert Ps.size == Ds.size
    assert Ps.size == Zs.size
    assert Ps.size == lnLs_transit.size
    return Ps, T0s, Ds, Zs, lnLs_transit


def remove_multiple_on_lnLs(bjd, ef, Ps, T0s, Ds, Zs, lnLs, rP=.01, rZ=.2):
    '''remove multiple orbital periods but dont assume the shortest one is
    correct, instead select the one with the highest lnL.'''
    assert Ps.size == T0s.size
    assert Ps.size == Ds.size
    assert Ps.size == Zs.size
    assert Ps.size == lnLs.size
    to_remove = np.zeros(0)
    for i in range(Ps.size):

        # check for positive multiples (eg 2P, 3P, 4P, ...)
	Ntransits = int((bjd.max()-bjd.min()) / Ps[i])
	lim = Ntransits+1 if Ntransits+1 > 2 else 3
	for j in range(2,lim):
	    # check positive multiples
	    isclose = np.isclose(Ps[i]*j, Ps, rtol=rP)
	    if np.any(isclose):
		# remove if nearby period has a lower lnL and has the same
                # depth within rZ (i.e. rZ=10%)
		iscloselnL = (lnLs[isclose] <= lnLs[i])#& \
                             #(abs(1-Zs[isclose]/Zs[i]) < SNRZ)
		to_remove = np.append(to_remove, Ps[isclose][iscloselnL])

        # check inverse multiples (eg P/2, P/3, ...)
        div = 2.
	Pmin = .1
        while Ps[i]/div >= Pmin:
            isclose = np.isclose(Ps, Ps[i]/div, rtol=rP)
            div += 1.
            if np.any(isclose):
                # remove if nearby period has a lower lnL and has the same
                # depth within rZ (i.e. rZ=10%)
                iscloselnL = (lnLs[isclose] <= lnLs[i])# & \
                             #(abs(Zs[isclose]-Zs[i])/ef[0] < SNRZ)
		to_remove = np.append(to_remove, Ps[isclose][iscloselnL])    
                
    to_remove = np.unique(to_remove)
    assert to_remove.size <= Ps.size
    to_remove_inds = np.where(np.in1d(Ps, to_remove))[0]
    Ps_final = np.delete(Ps, to_remove_inds)
    T0s_final = np.delete(T0s, to_remove_inds)
    Ds_final = np.delete(Ds, to_remove_inds)
    Zs_final = np.delete(Zs, to_remove_inds)
    lnLs_final = np.delete(lnLs, to_remove_inds)
    return Ps_final, T0s_final, Ds_final, Zs_final, lnLs_final


def remove_common_P(Ps, T0s, Ds, Zs, lnLs, rP=.2):
    assert Ps.size == T0s.size
    assert Ps.size == Ds.size
    assert Ps.size == Zs.size
    assert Ps.size == lnLs.size
    # remove common periods based on maximum likelihood
    sort = np.argsort(Ps)
    POIs, T0OIs, DOIs, ZOIs, lnLOIs = Ps[sort], T0s[sort], Ds[sort], \
                                      Zs[sort], lnLs[sort]
    POIs_red, T0OIs_red, DOIs_red, ZOIs_red, lnLOIs_red = np.zeros(0), \
                                                          np.zeros(0), \
							  np.zeros(0), \
                                                          np.zeros(0), \
							  np.zeros(0)
    for i in range(POIs.size):
        isclose = np.isclose(POIs, POIs[i], rtol=rP)
        if np.any(isclose):
            g = lnLOIs == lnLOIs[isclose].max()
            POIs_red = np.append(POIs_red, POIs[g])
            T0OIs_red = np.append(T0OIs_red, T0OIs[g])
            DOIs_red = np.append(DOIs_red, DOIs[g])
            ZOIs_red = np.append(ZOIs_red, ZOIs[g])
            lnLOIs_red = np.append(lnLOIs_red, lnLOIs[g])
    _,unique = np.unique(POIs_red, return_index=True)
    POIs_red, T0OIs_red, DOIs_red, ZOIs_red, lnLOIs_red = POIs_red[unique], \
							  T0OIs_red[unique], \
                                                          DOIs_red[unique], \
                                                          ZOIs_red[unique], \
                                                          lnLOIs_red[unique]
    return POIs_red, T0OIs_red, DOIs_red, ZOIs_red, lnLOIs_red

    

def consider_fractional_P(bjd, fcorr, ef, Ps, T0s, Ds, Zs, lnLs, Ms, Rs, Teff):
    '''given periods of interest, fit the transit model and compute lnL 
    for fractions of those periods as they may have been missed in the 
    linear search if 2 adjacent transits are not seen above the SNR 
    threshold.'''  
    assert Ps.size == T0s.size
    assert Ps.size == Ds.size
    assert Ps.size == Zs.size
    assert Ps.size == lnLs.size
    
    Ps2, T0s2, Ds2, Zs2, lnLs2 = np.zeros(0), np.zeros(0), np.zeros(0), \
                                 np.zeros(0), np.zeros(0)
    for i in range(Ps.size):
        Ps2 = np.append(Ps2, Ps[i])
        T0s2 = np.append(T0s2, T0s[i])
        Ds2 = np.append(Ds2, Ds[i])
        Zs2 = np.append(Zs2, Zs[i])
        lnLs2 = np.append(lnLs2, lnLs[i])
        div = 2.
        while Ps[i]/div >= .1:
            params = np.array([Ps[i]/div, T0s[i], Zs[i], Ds[i]])
            P, T0, Z, D, fmodel,_ = fit_params(params, bjd, fcorr, ef, 
                                               Ms, Rs, Teff)
            Ps2 = np.append(Ps2, P)
            T0s2 = np.append(T0s2, T0)
            Ds2 = np.append(Ds2, D)
            Zs2 = np.append(Zs2, Z)
	    lnLs2 = np.append(lnLs2, lnlike(bjd, fcorr, ef, fmodel))
            div += 1.
            
    return Ps2, T0s2, Ds2, Zs2, lnLs2


def identify_transit_candidates(self, Ps, T0s, Ds, Zs, lnLs, Ndurations, Rs,
                                bjd, fcorr, ef):
    '''Given the transit parameters and their lnLs, identify transit 
    candidates.'''
    assert Ps.size == T0s.size
    assert Ps.size == Ds.size
    assert Ps.size == Zs.size
    assert Ps.size == lnLs.size

    # remove negative entries (mostly negative depths) and only allow a
    # maximum number of periodicities for efficiency
    Nmax = 100
    g = (Ps>0) & (T0s>0) & (Ds>0) & (Zs>0) & (Zs<.9)
    if g.sum() > Nmax:
        Ps, T0s, Ds, Zs, lnLs = Ps[g], T0s[g], Ds[g], Zs[g], lnLs[g]
        s = np.argsort(lnLs)[::-1][:Nmax]
        Ps, T0s, Ds, Zs, lnLs = Ps[s], T0s[s], Ds[s], Zs[s], lnLs[s]
    else:
        Ps, T0s, Ds, Zs, lnLs = Ps[g], T0s[g], Ds[g], Zs[g], lnLs[g]
    
    # get optimized parameters to get more precise Ps and T0s which will help 
    # when removing multiples
    Ps2, T0s2, Ds2, Zs2 = np.zeros_like(Ps), np.zeros_like(Ps), \
		      	  np.zeros_like(Ps), np.zeros_like(Ps)
    lnLs2 = np.zeros_like(Ps)
    for i in range(Ps.size):
        params = np.array([Ps[i], T0s[i], Zs[i], Ds[i]])
        Ps2[i], T0s2[i], Zs2[i], Ds2[i], fmodel,_ = fit_params(params, bjd,
                                                               fcorr, ef,
						               self.Ms, self.Rs,
                                                               self.Teff)
        lnLs2[i] = lnlike(bjd, fcorr, ef, fmodel)

    # remove common periods based on maximum likelihood
    POIs1, T0OIs1, DOIs1, ZOIs1, lnLOIs1 = \
                                    remove_common_P(Ps2, T0s2, Ds2, Zs2, lnLs2)

    # update Z manually
    ZOIs2 = np.zeros(POIs1.size)
    for i in range(ZOIs2.size):
	phase = foldAt(bjd, POIs1[i], T0OIs1[i])
	phase[phase > .5] -= 1
	dur = .25*DOIs1[i]/POIs1[i]
	intransit = (phase >= -dur) & (phase <= dur)
	ZOIs2[i] = 1-np.median(fcorr[intransit])

    # remove multiple transits (i.e. 2P, 3P, 4P...)
    g = (ZOIs2>0) & (ZOIs2<.9)
    POIs2, T0OIs2, DOIs2, ZOIs2, lnLOIs2 = remove_multiple_on_lnLs(bjd, ef,
                                                                   POIs1[g],
                                                                   T0OIs1[g],
                                                                   DOIs1[g],
                                                                   ZOIs2[g],
                                                                   lnLOIs1[g])

    # do not consider too many planets to limit FPs
    g = ZOIs2 > 0
    params2 = np.array([POIs2, T0OIs2, ZOIs2, DOIs2]).T[g]
    params2, lnLOIs2 = trim_planets(params2, lnLOIs2[g])
    POIs2, T0OIs2, ZOIs2, DOIs2 = params2.T
    
    # consider integer fractions of Ps which may have been missed by
    # the linear search but should be 'detectable' when phase folded
    POIs3, T0OIs3, DOIs3, ZOIs3, lnLOIs3 = consider_fractional_P(bjd, fcorr, ef,
                                                                 POIs2, T0OIs2,
                                                                 DOIs2, ZOIs2,
                                                                 lnLOIs2,
                                                                 self.Ms,
                                                                 self.Rs,
                                                                 self.Teff)
    
    # remove duplicates and multiples
    POIs4, T0OIs4, DOIs4, ZOIs4, lnLOIs4 = remove_common_P(POIs3, T0OIs3, DOIs3,
                                                           ZOIs3, lnLOIs3)
    POIs5, T0OIs5, DOIs5, ZOIs5, lnLOIs5 = remove_multiple_on_lnLs(bjd, ef,
                                                                   POIs4,
                                                                   T0OIs4,
                                                                   DOIs4, ZOIs4,
                                                                   lnLOIs4)

    # do not consider too many planets to limit FPs
    g = ZOIs5 > 0
    params5 = np.array([POIs5, T0OIs5, ZOIs5, DOIs5]).T[g]
    params6, lnLOIs6 = trim_planets(params5, lnLOIs5[g])

    # identify bona-fide transit-like events
    self.params_guess_priorto_confirm, self.lnLOIs_priorto_confirm = params6, \
                                                                     lnLOIs6
    self._pickleobject()
    params6,Ntransits,lnLOIs6,cond_vals,conds = confirm_transits(params6, lnLOIs6,
                                                                 bjd, fcorr, ef,
                                                                 self.Ms, self.Rs,
                                                                 self.Teff)
    self.Ntransits = Ntransits
    self.transit_condition_free_params = np.array([dispersion_sig,
                                                   depth_sig,
                                                   bimodalfrac,
                                                   T0tolerance,
						   np.nan,
                                                   transitlikefrac,
                                                   min_autocorr_coeff])
    self.transit_condition_values = cond_vals
    self.transit_condition_bool = conds
    self.transit_condition_labels = np.array(['scatterin_gtr_scatterout',
                                              'depth_gtr_rms',
                                              'no_bimodal_flux_intransit',
                                              'flux_symmetric_in_time',
                                              'good_ephemeris',
                                              'indiv_transit_fraction',
                                              'not_autocorrelated_residuals'])

    # re-remove multiple transits based on refined parameters
    p,t0,d,z,lnLs = remove_common_P(params6[:,0], params6[:,1], params6[:,3],
                                    params6[:,2], lnLOIs6)
    p,t0,d,z,lnLs = remove_multiple_on_lnLs(bjd, ef, p, t0, d, z, lnLs)
    params = np.array([p,t0,z,d]).T

    # try to identify EBs
    params, EBparams, maybeEBparams, EBconditions, EBcondition_labels = \
                                    vett.identify_EBs(params, bjd, fcorr, ef,
                                                      self.Ms, self.Rs, 
						      self.Teff)
    self.EBconditions, self.EBcondition_labels = EBconditions, \
                                                 EBcondition_labels
    
    return params[:,0], params[:,1], params[:,3], params[:,2], lnLs, \
        params, EBparams, maybeEBparams


#def _optimize_box_transit(theta, bjd, fcorr, ef):
#    assert len(theta) == 4
#    #P, T0, depth, duration = theta
#    popt,_ = curve_fit(box_transit_model_curve_fit, bjd, fcorr, p0=theta, sigma=ef, absolute_sigma=True)
#    return popt


def get_LDcoeffs_Kepler(Ms, Rs, Teff, Z=0):
    '''Interpolate Claret+2012 grid of limb darkening coefficients to a
    given star.'''
    # get LD coefficient grid (Z is always 0 for some reason)
    clarlogg, clarTeff, clarZ, clar_a, clar_b = \
                                    np.loadtxt('LDcoeffs/claret12.tsv',
                                               delimiter=';', skiprows=40,
                                               usecols=(0,1,2,4,5)).T

    # interpolate to get the stellar LD coefficients
    logg = np.log10(6.67e-11*rvs.Msun2kg(Ms)*1e2 / rvs.Rsun2m(Rs)**2)
    lint_a = lint(np.array([clarTeff,clarlogg]).T, clar_a)
    lint_b = lint(np.array([clarTeff,clarlogg]).T, clar_b)

    return float(lint_a(Teff,logg)), float(lint_b(Teff,logg))


def get_LDcoeffs_TESS(Ms, Rs, Teff, Z=0):
    '''Interpolate Claret 2017 grid of limb darkening coefficients to a
    given star.'''
    # get LD coefficient grid (Z is always 0 for some reason)
    clarlogg, clarTeff, clarZ, clar_a, clar_b = \
                                    np.loadtxt('LDcoeffs/claret17.tsv',
                                               delimiter=';', skiprows=37).T

    # interpolate to get the stellar LD coefficients
    logg = np.log10(6.67e-11*rvs.Msun2kg(Ms)*1e2 / rvs.Rsun2m(Rs)**2)
    lint_a = lint(np.array([clarTeff,clarlogg]).T, clar_a)
    lint_b = lint(np.array([clarTeff,clarlogg]).T, clar_b)

    return float(lint_a(Teff,logg)), float(lint_b(Teff,logg))


def transit_model_func_curve_fit(u1, u2):
    def transit_model_func_in(bjd, P, T0, aRs, rpRs, inc):
        p = batman.TransitParams()
        p.t0, p.per, p.rp = 0., 1., rpRs
        p.a, p.inc, p.ecc = aRs, inc, 0.
        p.w, p.limb_dark, p.u = 90., 'quadratic', [u1,u2]
        phase = foldAt(bjd, P, T0)
        m = batman.TransitModel(p, phase)
        f = m.light_curve(p)
        return f
    return transit_model_func_in


def fit_params(params, bjd, fcorr, ef, Ms, Rs, Teff):
    '''Get best-fit parameters.'''
    assert params.shape == (4,)
    P, T0, depth, duration = params
    assert depth < .9  # sometimes the dimming is passed instead of depth 
    u1, u2 = get_LDcoeffs_Kepler(Ms, Rs, Teff)
    #u1, u2 = get_LDcoeffs_TESS(Ms, Rs, Teff)
    aRs = rvs.AU2m(rvs.semimajoraxis(P,Ms,0)) / rvs.Rsun2m(Rs)
    rpRs = np.sqrt(depth)
    p0 = P, T0, aRs, rpRs, 90.
    incs = np.array([float(rvs.inclination(P,Ms,Rs,1)), float(rvs.inclination(P,Ms,Rs,-1))])
    bnds = ((P*.9, T0-P*1.1, aRs*.9, 0, incs.min()),
            (P*1.1, T0+P*1.1, aRs*1.1, 1, incs.max()))
    try:
        popt,_ = curve_fit(transit_model_func_curve_fit(u1,u2),
                           bjd, fcorr, p0=p0, sigma=ef,
                           absolute_sigma=True, bounds=bnds)
        P, T0, aRs, rpRs, inc = popt
        depth = rpRs**2
        b = rvs.impactparam_inc(P, Ms, Rs, inc)
        duration = P/(np.pi*aRs) * np.sqrt((1+np.sqrt(depth))**2 - b*b)
	func = transit_model_func_curve_fit(u1, u2)
	fmodel = func(bjd, P, T0, aRs, rpRs, inc)
	params = np.array([P,T0,aRs,rpRs,inc,u1,u2])
    except RuntimeError, ValueError:
	func = transit_model_func_curve_fit(u1, u2)
	fmodel = func(bjd, P, T0, aRs, rpRs, 90.)
        P, T0, depth, duration = np.repeat(np.nan, 4)
	params = np.repeat(np.nan, 7)
    return P, T0, depth, duration, fmodel, params


def trim_planets(params, lnLOIs, Nplanetsmax=5):
    '''If there are too many planet candidates then remove some based on 
    their lnLs. Most are FPs.'''
    Nplanets = params.shape[0]
    assert params.shape == (Nplanets,4)
    assert lnLOIs.size == Nplanets
    if Nplanets > Nplanetsmax:
	tokeep = np.sort(np.argsort(lnLOIs)[::-1][:Nplanetsmax])  # retain the same ordering
	return params[tokeep], lnLOIs[tokeep]
    else:
	return params, lnLOIs


def compute_Ntransits(bjd, P, T0):
    possible_T0s = np.arange(-400,400)*P + T0
    i = 2
    while (possible_T0s.min() >= bjd.min()) or (possible_T0s.max() <= bjd.max()):
	possible_T0s = np.arange(-400*i,400*i)*P + T0
	i += 1
    Ntransits = np.sum((possible_T0s >= bjd.min()) & (possible_T0s <= bjd.max()))    
    return Ntransits

    
def is_not_autocorrelated(timeseries):
    x = np.ascontiguousarray(timeseries)
    n = x.size
    norm = (x - np.mean(x))
    result = np.correlate(norm, norm, mode='same')
    acorr = result[n//2 + 1:] / (x.var() * np.arange(n-1, n//2, -1))
    lag = np.abs(acorr).argmax() + 1
    r = acorr[lag-1]        
    if np.abs(r) <= min_autocorr_coeff:
        return True, np.abs(r)
    else: 
        return False, np.abs(r)
    
    
def confirm_transits(params, lnLs, bjd, fcorr, ef, Ms, Rs, Teff,
                     minNpnts_intransit=4):
    '''Look at proposed transits and confirm whether or not a significant 
    dimming is seen.'''
    Nplanets = params.shape[0]
    assert lnLs.size == Nplanets
    paramsout, to_remove_inds = np.zeros((Nplanets,4)), np.zeros(0)
    Ntransits = np.zeros((Nplanets))
    transit_condition_scatterin_val = np.zeros(Nplanets)
    transit_condition_scatterin_gtr_scatterout = np.zeros(Nplanets, dtype=bool)
    transit_condition_depth_val = np.zeros(Nplanets)
    transit_condition_depth_gtr_rms = np.zeros(Nplanets, dtype=bool)
    transit_condition_no_bimodal_val = np.zeros(Nplanets)
    transit_condition_no_bimodal_flux_intransit = np.zeros(Nplanets,dtype=bool)
    transit_condition_timesym_val = np.zeros(Nplanets)
    transit_condition_timesym = np.zeros(Nplanets,dtype=bool)
    transit_condition_indiv_transit_frac_val = np.zeros(Nplanets)
    transit_condition_indiv_transit_frac_gt_min = np.zeros(Nplanets, dtype=bool)
    transit_condition_ephemeris_fits_in_WF = np.zeros(Nplanets, dtype=bool)

    print 'Confirming proposed transits...'
    for i in range(Nplanets):
	print float(i) / Nplanets
	
	# try original first, then try optimized parameters if necessary
	# TEMP
	j = 0
	while j <= 1:
	    
	    if j == 0:
	    	P, T0, depth1, duration = params[i]
	    else:
		P, T0, depth1, duration,_,_ = fit_params(params[i], bjd, fcorr,
                                                         ef, Ms, Rs, Teff)

	    # get in and out of transit window
            phase = foldAt(bjd, P, T0)
            phase[phase > .5] -= 1
            # fraction of the duration in-transit
            # (should be <.5 to ignore ingress & egress)
            Dfrac = .25
            intransit = (phase*P >= -Dfrac*duration) & \
                        (phase*P <= Dfrac*duration)
            intransitfull = (phase*P >= -duration/2) & (phase*P <= duration/2)
	    Dfrac = .55
            outtransit = (phase*P <= -Dfrac*duration) | \
                         (phase*P >= Dfrac*duration)
            #plt.plot(phase, fcorr, 'ko', phase[intransit], fcorr[intransit],
            #         'bo'), plt.show()

            # calculate number of transits
            Ntransits[i] = compute_Ntransits(bjd, P, T0)
            
            # check scatter in and out of the proposed transit to see if the
            # transit is real
            cond1_val = (np.median(fcorr[outtransit]) - \
                         np.median(fcorr[intransit])) / \
                         MAD1d(fcorr[outtransit])
            cond1 = cond1_val > dispersion_sig
            transit_condition_scatterin_val[i] = cond1_val
            transit_condition_scatterin_gtr_scatterout[i] = cond1
	    # also check that the transit depth is significant relative to
            # the noise
            depth2 = 1-np.median(fcorr[intransit])
            sigdepth = np.std(fcorr[intransit])
            depth = depth2 #depth1
            cond5 = depth > 0
            cond2_val = depth / sigdepth * np.sqrt(Ntransits[i])
            cond2 = cond2_val > depth_sig
            transit_condition_depth_val[i] = cond2_val
	    transit_condition_depth_gtr_rms[i] = cond2
	    # ensure that the flux measurements intransit are not bimodal
            # (ie. f at depth and at f=1 which would indicate a 
	    # bad period and hence a FP
            Npnts_intransit = fcorr[intransit].size
            Npnts_lt_1sigdepth = (fcorr[intransit] <= 1-depth+sigdepth).sum()
	    if (Npnts_intransit >= minNpnts_intransit):
                cond3_val = float(Npnts_lt_1sigdepth) / Npnts_intransit
                cond3 = cond3_val >= bimodalfrac
            else:
            	cond3_val, cond3 = np.nan, False
	    transit_condition_no_bimodal_val[i] = cond3_val
	    transit_condition_no_bimodal_flux_intransit[i] = cond3

	    # ensure that the flux measurements over the full transit are
            # symmetrical in time
            Npnts_intransitfull = fcorr[intransitfull].size
            Npnts_lt_T0 = (phase[intransitfull] < 0).sum()
	    if Npnts_intransitfull >= minNpnts_intransit:
                cond4_val = float(Npnts_lt_T0) / Npnts_intransitfull
                cond4 = .5-T0tolerance <= cond4_val <= .5+T0tolerance
	    else:
            	cond4_val, cond4 = np.nan, False
            transit_condition_timesym_val[i] = cond4_val
	    transit_condition_timesym[i] = cond4

            # check that most individual transits look like transits
            Dfrac = .5
            events_BJD = np.arange(-1000,1000)*P + T0
            g = (events_BJD >= bjd.min()) & (events_BJD <= bjd.max())
            events_BJD = events_BJD[g]
            depth_gtr_sig = np.zeros(events_BJD.size)
            aretransits = np.zeros(events_BJD.size, dtype=bool)
            for k in range(depth_gtr_sig.size):
                t0 = events_BJD[k]
                phase = foldAt(bjd, P, t0)
                phase[phase>.5] -= 1
                intransit = (bjd >= t0-Dfrac*duration) & \
                            (bjd <= t0+Dfrac*duration)
                outtransit1 = (bjd >= t0-10*duration) & \
                              (bjd <= t0-8*duration)
                outtransit2 = (bjd <= t0+10*duration) & \
                              (bjd >= t0+8*duration)
                ##plt.plot(bjd, fcorr, 'o', bjd[outtransit1],
                ##         fcorr[outtransit1], 'o'), plt.show()
                fin = np.median(fcorr[intransit])
                sigdepth = np.std(fcorr[intransit]) if intransit.sum() > 1 \
                           else np.median(ef[intransit])
                fout = np.median(np.append(fcorr[outtransit1],
                                           fcorr[outtransit2]))
                aretransits[k] = fin+sigdepth < fout

            cond6_val = aretransits.sum() / float(aretransits.size)
            cond6 = cond6_val >= transitlikefrac
            transit_condition_indiv_transit_frac_val[i] = cond6_val
            transit_condition_indiv_transit_frac_gt_min[i] = cond6
    
            # ensure that at least two transits will fit within the observing
            # window otherwise its just a
            # single transit-like event
            cond7 = ((T0-P >= bjd.min()) | (T0+P <= bjd.max())) & \
                    (T0 >= bjd.min()) & (T0 <= bjd.max()) & \
                    (P < bjd.max()-bjd.min())
            transit_condition_ephemeris_fits_in_WF[i] = cond7
            paramsout[i] = P, T0, depth, duration
            if cond1 and cond2 and cond3 and cond4 and cond5 and cond6 and \
               cond7:
	        j += 2
	        pass
            else:
                to_remove_inds = np.append(to_remove_inds, i)
                j += 1

    # only remove parameter sets that fail both the original and optimized tests
    to_remove_inds, counts = np.unique(to_remove_inds, return_counts=True)
    to_remove_inds = to_remove_inds[counts==2]

    # remove false transits
    paramsout = np.delete(paramsout, to_remove_inds, 0)
    lnLsout = np.delete(lnLs, to_remove_inds)

    # check autocorrelation after removing planets (is autocorrelated if the
    # systematic correction is bad -> planets detections are not robust)
    fmodel_tot = np.zeros(bjd.size)
    for i in range(paramsout.shape[0]):
        _,_,_,_,fmodel,_ = fit_params(paramsout[i], bjd, fcorr,
                                      ef, Ms, Rs, Teff)
        fmodel_tot += fmodel - 1.
    cond_autocorrs, autocorr_coeffs = np.zeros(2, dtype=bool), np.zeros(2)
    cond_autocorrs[0], autocorr_coeffs[0] = is_not_autocorrelated(fcorr - \
                                                                  fmodel_tot)
    cond_autocorrs[1], autocorr_coeffs[1] = is_not_autocorrelated(fcorr)
    transit_condition_autocorr_leq_max = np.repeat(np.any(cond_autocorrs),
                                                   Nplanets)
    transit_condition_autocorr_val = np.repeat(autocorr_coeffs.min(), Nplanets)
    if np.all(cond_autocorrs == False):
        paramsout, lnLsout = np.zeros((0,4)), np.zeros(0)
 
    # combine conditions
    cond_vals = np.array([transit_condition_scatterin_val, \
                          transit_condition_depth_val, \
                          transit_condition_no_bimodal_val, \
                          transit_condition_timesym_val, \
                          transit_condition_ephemeris_fits_in_WF, \
                          transit_condition_indiv_transit_frac_val, \
                          transit_condition_autocorr_val]).T
    
    cond_bool = np.array([transit_condition_scatterin_gtr_scatterout, \
                          transit_condition_depth_gtr_rms, \
                          transit_condition_no_bimodal_flux_intransit, \
                          transit_condition_timesym, \
                          transit_condition_ephemeris_fits_in_WF, \
                          transit_condition_indiv_transit_frac_gt_min, \
                          transit_condition_autocorr_leq_max]).T

    return paramsout, Ntransits, lnLsout, cond_vals, cond_bool


def identify_EBs(params, bjd, fcorr, ef, Rs, SNRthresh=3., rpmax=30):
    '''For each proposed planet in params, check if there is a clearly-defined 
    secondary eclipse as is indicative of a secondary eclipse. Could also have 
    a V-shaped "transit" but so highly inclined transiting planets.'''
    Nplanets = params.shape[0]
    notEB = np.ones(Nplanets)
    for i in range(Nplanets):
	# check for significant eclipse depths at two times in case T0 is near an edge of the WF
	P, T0, depth, duration = params[i]
	eclipse1 = (bjd >= T0+.5*P-.5*duration) & (bjd <= T0+.5*P+.5*duration)
	eclipse2 = (bjd >= T0-.5*P-.5*duration) & (bjd <= T0-.5*P+.5*duration)	
	outeclipse1 = (bjd >= T0+.5*P-2*duration) & (bjd <= T0+.5*P-duration)
        outeclipse2 = (bjd >= T0-.5*P-2*duration) & (bjd <= T0-.5*P-duration)

	rms_ineclipse1  = fcorr[eclipse1].std()
	rms_outeclipse1 = fcorr[outeclipse1].std()
        rms_ineclipse2  = fcorr[eclipse2].std()
        rms_outeclipse2 = fcorr[outeclipse2].std()
        if (rms_ineclipse1 >= SNRthresh*rms_outeclipse1) or (rms_ineclipse2 >= SNRthresh*rms_outeclipse2):
            notEB[i] = 0.

    # check that planets are not too big
    rpRss = np.sqrt(params[:,2])
    rps = rvs.m2Rearth(rvs.Rsun2m(rpRss*Rs))
    notEB[(rpRss > .5) | (rps > rpmax)] = 0

    # save planet and EB parameters
    notEB = notEB.astype(bool)
    params, EBparams = params[notEB], params[np.invert(notEB)]
    return params, EBparams 
