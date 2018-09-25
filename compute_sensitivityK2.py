from K2LCclass import *
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


def sample_planets_uniform(bjd, Ms, Rs, Teff, Plims=(.1,30), rplims=(.5,4)):
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
    while np.any(np.diff(Ptrue) / Ptrue[1:] <= .1):
        Ptrue = np.sort(10**np.random.uniform(np.log10(Pmin), np.log10(Pmax),
                                              Nplanets))
    T0true = np.zeros(Nplanets)
    while np.any((T0true <= bjd.min()) | (T0true >= bjd.max())):
        T0true = np.random.choice(bjd, Nplanets) + np.random.randn(Nplanets)*.5
    rptrue = np.random.uniform(rpmin, rpmax, Nplanets)
    rpRs = rvs.Rearth2m(rptrue) / rvs.Rsun2m(Rs)
    aRs = rvs.AU2m(rvs.semimajoraxis(Ptrue, Ms, 0)) / rvs.Rsun2m(Rs)
    assert np.all(aRs > 1)
    bs = np.random.uniform(-1, 1, Nplanets)
    incs = rvs.inclination(Ptrue, Ms, Rs, bs)
    u1, u2 = llnl.get_LDcoeffs_Kepler(Ms, Rs, Teff)

    print '\n', Ptrue, T0true, aRs, rpRs, incs, u1, u2
    
    # compute transit model(s)
    fmodel = np.ones(bjd.size)
    depthtrue, durationtrue = np.zeros(Nplanets), np.zeros(Nplanets)
    for i in range(Nplanets):
        fmodeltmp = transit_model_func_in(bjd, Ptrue[i], T0true[i], aRs[i],
                                          rpRs[i], incs[i], u1, u2) - 1.
        fmodel += fmodeltmp
        depthtrue[i] = abs(fmodeltmp.min())
        #prior2transit = (bjd >= T0true[i]-Ptrue[i]*.5) & (bjd < T0true[i])
        #aftertransit = (bjd <= T0true[i]+Ptrue[i]*.5) & (bjd > T0true[i])
        #T1 = bjd[prior2transit][fmodeltmp[prior2transit]==0][-1]
        #T4 = bjd[aftertransit][fmodeltmp[aftertransit]==0][0]        
        #durationtrue[i] = T4-T1
        durationtrue[i] = rvs.transit_width(Ptrue[i], Ms, Rs, rptrue[i], bs[i])
        
    return Ptrue, T0true, depthtrue, durationtrue, rptrue, fmodel
    

def injected_planet_search(epicnum, index):
    '''Inject planets into a K2 light curve using the pipeline defined in
    planet_search to search for planets.'''

    # get data and only run the star if it is of interest
    if not is_star_of_interest(epicnum):
        return None
    name, Kepmag, logg, Ms, Rs, Teff, bjd, f, ef = read_K2_data(epicnum)
    self = K2LC('%s_%.4d'%(name, index))
    self.bjd, self.f, self.ef = bjd, f, ef
    self.Kepmag, self.logg, self.Ms, self.Rs, self.Teff = Kepmag, logg, Ms, \
                                                          Rs, Teff
    self.DONE = False
    self._pickleobject()

    # sample and inject planet(s)
    Ptrue, T0true, depthtrue, durationtrue, rptrue, fmodel = \
                                    sample_planets_uniform(bjd, Ms, Rs, Teff)
    self.params_true = np.array([Ptrue, T0true, depthtrue, durationtrue]).T
    self.Ptrue, self.rptrue = Ptrue, rptrue
    f += fmodel - 1.
    self.f = f
    self._pickleobject()
    
    # fit initial GP
    thetaGPall, resultsGPall, thetaGPin, thetaGPout = do_optimize_0(bjd, f, ef)
    self.thetaGPall, self.resultsGPall = thetaGPall, resultsGPall
    self.thetaGPin, self.thetaGPout = thetaGPin, thetaGPout
    self._pickleobject()

    # search for transits in the corrected LC and get the transit parameters
    # guesses
    print 'Searching for transit-like events...\n'
    params, EBparams, maybeEBparams = find_transits(self, bjd, f, ef,
                                                    thetaGPout)
    self.params_guess = params
    self.params_guess_labels = np.array(['Ps', 'T0s', 'depths [Z]', \
                                         'durations [D]'])
    self.EBparams_guess, self.maybeEBparams_guess = EBparams, maybeEBparams
    self._pickleobject()

    # check if planets are detected
    self.is_detected = np.array([int(np.any(np.isclose(params[:,0], Ptrue[i],
                                                       rtol=.02)))
                                 for i in range(Ptrue.size)]).astype(bool)
    self.DONE = True
    self._pickleobject()



def do_i_run_this_sim(epicnum, index):
    # check if star is already done
    fname = 'PipelineResults/EPIC_%i_%.4d/K2LC'%(epicnum, index)
    if os.path.exists(fname):
        return not loadpickle(fname).DONE
    else:
        return True


if __name__ == '__main__':
    startind = int(sys.argv[1])
    endind = int(sys.argv[2])
    Nsystems = int(sys.argv[3])
    epics= np.loadtxt(K2Mdwarffile, delimiter=',')[:,0]
    for i in range(startind, endind):
        print epics[i]
        for j in range(Nsystems):
            if do_i_run_this_sim(epics[i], j):
                injected_planet_search(epics[i], j)
