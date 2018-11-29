from LCclass import *
from planetsearch import *
import linear_lnlike as llnl


def transit_model_func_in(bjd, P, T0, aRs, rpRs, inc, u1, u2):
    p = batman.TransitParams()
    p.t0, p.per, p.rp = 0., 1., rpRs
    p.a, p.inc, p.ecc = aRs, inc, 0.
    p.w, p.limb_dark, p.u = 90., 'quadratic', [u1,u2]
    phase = foldAt(bjd, P, T0)
    m = batman.TransitModel(p, phase)
    f = m.light_curve(p)
    return f


def remove_detected_planets(folder, IDnum, prefix, bjd, f):
    '''Remove planet detections using their optimized parameters to clean the 
    light curve before searching for injected planets.'''
    # get planet search results
    try:
        d = loadpickle('%s/%s_%i/LC_-00099'%(folder,prefix,IDnum))
    except IOError:
        raise ValueError('initial planet search has not been run.')

    # get planet transit models
    fmodel = np.ones(bjd.size)
    for i in range(d.Ndet):
        func = llnl.transit_model_func_curve_fit(d.u1, d.u2)
        P, T0, aRs, rpRs, inc = d.params_optimized[i]
        fmodel *= func(bjd, P, T0, aRs, rpRs, inc)

    # remove planets to clean the light curve 
    f_noplanets = f / fmodel
    return f_noplanets


def sample_planets_uniform(bjd, Ms, Rs, Teff, Plims=(.5,200), rplims=(.5,10)):
    '''Sample M dwarf planets over a log uniform grid.'''
    Nplanets = 0
    while Nplanets < 1:
        Nplanets = int(np.round(np.random.normal(2.5,.2))) # from DC15

    Pmin, Pmax = Plims
    rpmin, rpmax = rplims
    assert Pmin < Pmax
    assert rpmin < rpmax

    # get planets
    Ptrue = np.ones(Nplanets)
    while np.any(np.diff(Ptrue) / Ptrue[1:] <= .1) or not np.any(Ptrue < 80):
        Ptrue = np.sort(10**np.random.uniform(np.log10(Pmin), np.log10(Pmax),
                                              Nplanets))
    T0true = np.zeros(Nplanets)
    while np.any((T0true <= bjd.min()) | (T0true >= bjd.max())):
        T0true = np.random.choice(bjd, Nplanets) + np.random.randn(Nplanets)*.5
    rptrue = 10**np.random.uniform(np.log10(rpmin), np.log10(rpmax), Nplanets)
    rpRs = rvs.Rearth2m(rptrue) / rvs.Rsun2m(Rs)
    aRs = rvs.AU2m(rvs.semimajoraxis(Ptrue, Ms, 0)) / rvs.Rsun2m(Rs)
    assert np.all(aRs > 1)
    bs = np.random.uniform(-.7, .7, Nplanets)
    incs = rvs.inclination(Ptrue, Ms, Rs, bs)
    u1, u2 = llnl.get_LDcoeffs_Kepler(Ms, Rs, Teff)
    
    # compute transit model(s)
    fmodel = np.ones(bjd.size)
    depthtrue, durationtrue = np.zeros(Nplanets), np.zeros(Nplanets)
    for i in range(Nplanets):
        fmodeltmp = transit_model_func_in(bjd, Ptrue[i], T0true[i], aRs[i],
                                          rpRs[i], incs[i], u1, u2)
        fmodel *= fmodeltmp
        depthtrue[i] = 1 - abs(fmodeltmp.min())
        durationtrue[i] = rvs.transit_width(Ptrue[i], Ms, Rs, rptrue[i], bs[i])
        
    return Ptrue, T0true, depthtrue, durationtrue, rptrue, fmodel



def injected_planet_search(folder, IDnum, index, K2=False, Kep=False, TESS=False):
    '''Inject planets into a K2 light curve using the pipeline defined in
    planet_search to search for planets.'''

    # get data and only run the star if it is of interest
    if not is_star_of_interest(IDnum, Kep=Kep, K2=K2):
        return None

    if K2:
        name, star_dict, bjd, f, ef, quarters = read_K2_data(IDnum)
	prefix, Nopt = 'EPIC', 10
    elif Kep:
	name, star_dict, bjd, f, ef, quarters = read_Kepler_data(IDnum)
	prefix, Nopt = 'KepID', 5
    elif TESS:
	name, star_dict, bjd, f, ef, quarters = read_TESS_data(IDnum)
	prefix, Nopt = 'TIC', 10
    else:
	return None
    self = LCclass(folder, name, index)
    self.bjd, self.f_orig, self.ef, self.quarters = bjd, np.copy(f), ef, quarters
    for attr in star_dict.keys():
        setattr(self, attr, star_dict[attr])
    self.DONE = False
    self._pickleobject()

    # remove planets detected by the planet search which should already
    # have been run
    self.f_noplanets = remove_detected_planets(folder, IDnum, prefix, self.bjd, self.f_orig)
    
    # sample and inject planet(s)
    Ptrue, T0true, depthtrue, durationtrue, rptrue, fmodel = \
                                    sample_planets_uniform(bjd, self.Ms,
                                                           self.Rs, self.Teff)
    self.params_true = np.array([Ptrue, T0true, depthtrue, durationtrue]).T
    self.Ptrue, self.rptrue = Ptrue, rptrue
    self.f = self.f_noplanets * fmodel
    self._pickleobject()
    
    # read-in initial GP
    d = loadpickle('%s/%s/LC_-00099'%(folder, name))
    self.thetaGPin, self.thetaGPout = d.thetaGPin, d.thetaGPout
    self._pickleobject()
    
    # search for transits in the corrected LC and get the transit parameters
    # guesses
    print 'Searching for transit-like events...\n'
    if K2 or Kep:
	Kep, TESS = True, False
    else:
	Kep, TESS = False, True
    params, EBparams, maybeEBparams = find_transits(self, self.bjd, self.f,
                                                    self.ef, self.quarters,
						    self.thetaGPout, Kep=Kep
						    TESS=TESS)
    self.params_guess = params
    self.params_guess_labels = np.array(['Ps', 'T0s', 'depths [Z]', \
                                         'durations [D]'])
    self.EBparams_guess, self.maybeEBparams_guess = EBparams, maybeEBparams
    self._pickleobject()

    # save depth and SNR for recovered planets plus the SNR of injected
    p = compute_SNRtransit(self.bjd, self.fcorr, self.params_guess)
    self.SNRtransits, self.depths, self.sigtransit = p
    Ntransits = np.array([llnl.compute_Ntransits(self.bjd,Ptrue[i],T0true[i])
                          for i in range(Ptrue.size)])
    self.SNRtrue = (depthtrue / self.sigtransit) * np.sqrt(Ntransits)
    assert self.SNRtrue.size == self.Ptrue.size
    
    # check if planets are detected
    self.is_detected = np.array([int(np.any(np.isclose(params[:,0], Ptrue[i],
                                                       rtol=.02)))
                                 for i in range(Ptrue.size)]).astype(bool)
    self.is_FP = np.array([int(np.invert(np.any(np.isclose(Ptrue, params[i,0],
                                                           rtol=.02))))
                           for i in range(params.shape[0])]).astype(bool)
    self.DONE = True
    self._pickleobject()



def do_i_run_this_sim(folder, IDnum, prefix, index):
    # check if planet search has already been run
    fname = '%s/%s_%i/LC_-00099'%(folder, prefix, IDnum)
    if os.path.exists(fname) and loadpickle(fname).DONE:
        
        # check if star is already done
        fname = '%s/%s_%i/LC_%.5d'%(folder, prefix, IDnum, index)
        if os.path.exists(fname):
            return not loadpickle(fname).DONE
        else:
            return True

    else:
        return False
        

if __name__ == '__main__':
    startind = int(sys.argv[1])
    Nstars = int(sys.argv[2])
    Njobs = int(sys.argv[3])
    Nsystems_per_star = int(sys.argv[4])
    folder = sys.argv[5]
    f = glob.glob('%s/*_detectionsfirst.csv'%folder)
    ids = np.loadtxt(f[0])
    if 'EPIC' in folder:
	prefix, K2, Kep, TESS = 'EPIC', True, False, False
    if 'KepID' in folder:
	prefix, K2, Kep, TESS = 'KepID', False, True, False
    if 'TIC' in folder:
        prefix, K2, Kep, TESS = 'TIC', False, False, True
        
    for i in range(startind, Nstars, Njobs):
        print ids[i]
        for j in range(Nsystems_per_star):
            index = j + 0#starting_index
            if do_i_run_this_sim(folder, ids[i], prefix, index):
                injected_planet_search(folder, ids[i], index,
                                       Kep=Kep, K2=K2, TESS=TESS)
