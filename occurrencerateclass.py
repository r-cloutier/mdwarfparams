from LCclass import *
from sensclass import *
from uncertainties import unumpy as unp
import rvs, time
from scipy.ndimage.filters import gaussian_filter # for map smoothing if desired
from scipy.interpolate import LinearNDInterpolator as lint
from priors import get_results
from estimate_PDF import *


class OccurrenceRateclass:

    def __init__(self, folder, prefix, startstarind, endstarind,
                 xlen=20, ylen=12, compute_detections=False,
                 compute_sens=False, compute_occurrence_rate=False,
                 fine_factor=6, Plims=(.5,1e2), Flims=(.1,4e2),
                 smalims=(5e-3,.5), rplims=(.5,10)):
        self.folder, self.prefix = folder, prefix
        if 'EPIC' in self.folder:
            self.Kep, self.K2, self.TESS = True, False, False
        if 'KepID' in self.folder:
            self.Kep, self.K2, self.TESS = False, True, False
        if 'TIC' in self.folder: 
            self.Kep, self.K2, self.TESS = False, False, True
        self._startstarind, self._endstarind = int(startstarind), int(endstarind)
        assert self._startstarind < self._endstarind
        self.fname_out = '%s/%s_results_%i_%i'%(self.folder, self.prefix,
                                                self._startstarind,
                                                self._endstarind)
        
        self._xlen, self._ylen = int(xlen), int(ylen)
	self._fine_factor = int(fine_factor)
        self.Plims, self.Flims, self.smalims = Plims, Flims, smalims
        self.rplims = rplims

        self._t0 = time.time()
        
        if compute_detections:
	    self.fname_out += '_det'
            self.get_planetsearch_results()
            self.save_ids_detectionsfirst()
            self._pickleobject()
            
        if compute_sens:
	    self.fname_out += '_sens'
            self.get_simulation_results()
            self._pickleobject()

        if compute_occurrence_rate:
	    self.fname_out += '_occ'
            self.compute_occurrence_rate()
	    self._pickleobject()


    def get_planetsearch_results(self):
        '''Get the results from the planetsearch, i.e. the detected planets 
        and stellar properties.'''
        fs = np.array(glob.glob('%s/%s_*/LC_-00099'%(self.folder, self.prefix)))
        fs = fs[self._startstarind:self._endstarind]
        if fs.size == 0:
            return None
	self.fs_planetsearch, self.names_planetsearch = [], np.zeros(0)

        # detected planet params
        self.Ndetected, self.params_guess = np.zeros(0), np.zeros((0,4))
	self.params_optimized = np.zeros((0,5))
        self.Ps, self.e_Ps = np.zeros(0), np.zeros(0)
        self.T0s, self.e_T0s = np.zeros(0), np.zeros(0)
        self.aRss, self.ehi_aRss, self.elo_aRss = np.zeros(0), np.zeros(0), \
                                                  np.zeros(0)
        self.rpRss, self.ehi_rpRss, self.elo_rpRss = np.zeros(0), np.zeros(0), \
                                                     np.zeros(0)
        self.incs, self.ehi_incs, self.elo_incs = np.zeros(0), np.zeros(0), \
                                                  np.zeros(0)
        self.rps, self.ehi_rps, self.elo_rps = np.zeros(0), np.zeros(0), \
                                               np.zeros(0)
        self.smas, self.ehi_smas, self.elo_smas = np.zeros(0), np.zeros(0), \
                                                  np.zeros(0)
        self.Fs, self.ehi_Fs, self.elo_Fs = np.zeros(0),np.zeros(0),np.zeros(0)
        self.CDPPs, self.depths, self.SNRtransits = np.zeros(0), np.zeros(0), \
                                                    np.zeros(0)
        # POI params
        self.cond_vals, self.cond_free_params = np.zeros((0,9)), np.zeros((0,11))

        # are EBs detected
        self.isEB = np.zeros(0, dtype=int)
        
        # stellar params
        if self.Kep or self.K2:
            self.Kepmags = np.zeros(0)
        if self.TESS:
            self.TESSmags = np.zeros(0)

        self.GBPmags, self.e_GBPmags = np.zeros(0), np.zeros(0)
        self.GRPmags, self.e_GRPmags = np.zeros(0), np.zeros(0)
        self.Jmags, self.e_Jmags = np.zeros(0), np.zeros(0)
        self.Hmags, self.e_Hmags = np.zeros(0), np.zeros(0)
        self.Kmags, self.e_Kmags = np.zeros(0), np.zeros(0)
        self.efs = np.zeros(0)
	self.ras, self.decs = np.zeros(0), np.zeros(0)
        self.pars,self.e_pars = np.zeros(0),np.zeros(0)
        self.dists,self.ehi_dists,self.elo_dists = np.zeros(0),np.zeros(0),np.zeros(0)
        self.MKs,self.ehi_MKs,self.elo_MKs = np.zeros(0),np.zeros(0),np.zeros(0)
        self.Mss,self.ehi_Mss,self.elo_Mss = np.zeros(0),np.zeros(0),np.zeros(0)
        self.Rss,self.ehi_Rss,self.elo_Rss = np.zeros(0),np.zeros(0),np.zeros(0)
        self.Teffs, self.ehi_Teffs, self.elo_Teffs = np.zeros(0), np.zeros(0), \
                                                     np.zeros(0)
        self.loggs, self.ehi_loggs, self.elo_loggs = np.zeros(0), np.zeros(0), \
                                                     np.zeros(0)
        self.Lss,self.ehi_Lss,self.elo_Lss = np.zeros(0),np.zeros(0),np.zeros(0)
       	for i in range(fs.size):

            print float(i)/fs.size, fs[i]
	    self._conditional_pickleobject(self._t0, time.time())
            try:
             	d = loadpickle(fs[i])
            except (EOFError, ValueError):
                pass
	    if d.DONE:
                
		for j in range(d.Ndet+1):

                    self.fs_planetsearch.append(fs[i])
                    self.names_planetsearch = \
                                        np.append(self.names_planetsearch,
                                                  d.object_name)
                    self.Ndetected = np.append(self.Ndetected, d.Ndet)
                    if self.Kep or self.K2:
                        self.Kepmags = np.append(self.Kepmags, d.Kepmag)
                    if self.TESS:
                        self.TESSmags = np.append(self.TESSmags, d.TESSmag)

                    self.GBPmags = np.append(self.GBPmags, d.GBPmag)
                    self.e_GBPmags = np.append(self.e_GBPmags, d.e_GBPmag)
                    self.GRPmags = np.append(self.GRPmags, d.GRPmag)
                    self.e_GRPmags = np.append(self.e_GRPmags, d.e_GRPmag)
                    self.Jmags = np.append(self.Jmags, d.Jmag)
                    self.e_Jmags = np.append(self.e_Jmags, d.e_Jmag)
                    self.Hmags = np.append(self.Hmags, d.Hmag)
                    self.e_Hmags = np.append(self.e_Hmags, d.e_Hmag)
                    self.Kmags = np.append(self.Kmags, d.Kmag)
                    self.e_Kmags = np.append(self.e_Kmags, d.e_Kmag)
                    
                    self.efs = np.append(self.efs, d.ef.mean())
		    self.ras = np.append(self.ras, d.ra)
		    self.decs = np.append(self.decs, d.dec)
                    self.pars = np.append(self.pars, d.par)
                    self.e_pars = np.append(self.e_pars, d.e_par)
                    self.dists = np.append(self.dists, d.dist)
                    self.ehi_dists = np.append(self.ehi_dists, d.ehi_dist)
                    self.elo_dists = np.append(self.elo_dists, d.elo_dist)
                    self.MKs = np.append(self.MKs, d.MK)
                    self.ehi_MKs = np.append(self.ehi_MKs, d.ehi_MK)
                    self.elo_MKs = np.append(self.elo_MKs, d.elo_MK)
                    self.Mss = np.append(self.Mss, d.Ms)
                    self.ehi_Mss = np.append(self.ehi_Mss, d.ehi_Ms)
                    self.elo_Mss = np.append(self.elo_Mss, d.elo_Ms)
                    self.Rss = np.append(self.Rss, d.Rs)
                    self.ehi_Rss = np.append(self.ehi_Rss, d.ehi_Rs)
                    self.elo_Rss = np.append(self.elo_Rss, d.elo_Rs)
                    self.Teffs = np.append(self.Teffs, d.Teff)
                    self.ehi_Teffs = np.append(self.ehi_Teffs, d.ehi_Teff)
                    self.elo_Teffs = np.append(self.elo_Teffs, d.elo_Teff)
                    self.loggs = np.append(self.loggs, d.logg)
                    self.ehi_loggs = np.append(self.ehi_loggs, d.ehi_logg)
                    self.elo_loggs = np.append(self.elo_loggs, d.elo_logg)
                    _,_,samp_Ls = sample_Ls(d.object_name, self.prefix) 
                    Lss = get_results(samp_Ls.reshape(samp_Ls.size,1))
                    self.Lss = np.append(self.Lss, Lss[0])
                    self.ehi_Lss = np.append(self.ehi_Lss, Lss[1])
                    self.elo_Lss = np.append(self.elo_Lss, Lss[2])

                    isEB = int(np.all(d.EBconditions[j-1])) if j > 0 \
                           else np.nan
                    self.isEB = np.append(self.isEB, isEB)
                
                    # save parameter guesses of detected planets
                    params = d.params_guess[j-1] if j > 0 \
                             else np.repeat(np.nan,4)
		    self.params_guess = np.append(self.params_guess,
                                                  params.reshape(1,4), axis=0)

                    # save optimized parameters
                    params_opt = d.params_optimized[j-1] if j > 0 \
                                 else np.repeat(np.nan,5)
                    params_res = d.params_results[j-1] if j > 0 \
                                 else np.repeat(np.nan,15).reshape(3,5)
                    assert params_res.shape == (3,5)
                    self.params_optimized = np.append(self.params_optimized,
                                                      params_opt.reshape(1,5),
                                                      axis=0)

                    # save transit vetting results
                    P = params[0]
                    Pss = d.params_guess_priorto_confirm[:,0]
                    if Pss.size > 0:
                        g = abs(Pss-P) == np.min(abs(Pss-P))
                        assert g.sum() in range(2)
                        cond_vals = d.transit_condition_values[g] if j > 0 \
                                    else np.repeat(np.nan,9)
                    else:
                        cond_vals = np.repeat(np.nan, 9)
		    self.cond_vals = np.append(self.cond_vals,
                                               cond_vals.reshape(1,9), axis=0)
                    self.cond_free_params = np.append(self.cond_free_params,
                                d.transit_condition_free_params.reshape(1,11),
                                                      axis=0)

                    # save planet params
                    self.Ps = np.append(self.Ps, params_opt[0])
                    self.e_Ps = np.append(self.e_Ps,
                                          get_1sigma(params_res[:,0]))
                    self.T0s = np.append(self.T0s, params_opt[1])
                    self.e_T0s = np.append(self.e_T0s,
                                           get_1sigma(params_res[:,1]))
                    if j > 0:
		 	_,_,samp_aRs = get_samples_from_percentiles(*params_res[:,2])
                    	_,_,samp_rpRs = get_samples_from_percentiles(*params_res[:,3])
                    	_,_,samp_inc = get_samples_from_percentiles(*params_res[:,4])
		    else:
			samp_aRs = np.repeat(np.nan, 1000)
			samp_rpRs = np.repeat(np.nan, 1000)
                        samp_inc = np.repeat(np.nan, 1000)
                    v = np.percentile(samp_aRs, (16,50,84))
                    self.aRss = np.append(self.aRss, v[1])
                    self.ehi_aRss = np.append(self.ehi_aRss, v[2]-v[1])
                    self.elo_aRss = np.append(self.elo_aRss, v[1]-v[0])
                    v = np.percentile(samp_rpRs, (16,50,84))
                    self.rpRss = np.append(self.rpRss, v[1])
                    self.ehi_rpRss = np.append(self.ehi_rpRss, v[2]-v[1])
                    self.elo_rpRss = np.append(self.elo_rpRss, v[1]-v[0])
                    v = np.percentile(samp_inc, (16,50,84))
                    self.incs = np.append(self.incs, v[1])
                    self.ehi_incs = np.append(self.ehi_incs, v[2]-v[1])
                    self.elo_incs = np.append(self.elo_incs, v[1]-v[0])
                    samp_Rs = np.random.randn(samp_aRs.size)*d.ehi_Rs + d.Rs
                    samp_rp = rvs.m2Rearth(rvs.Rsun2m(samp_Rs*samp_rpRs))
                    v = np.percentile(samp_rp, (16,50,84))
                    self.rps = np.append(self.rps, v[1])
                    self.ehi_rps = np.append(self.ehi_rps, v[2]-v[1])
                    self.elo_rps = np.append(self.elo_rps, v[1]-v[0])
                    samp_Ms = np.random.randn(samp_aRs.size)*d.ehi_Ms + d.Ms
                    samp_sma = rvs.semimajoraxis(np.random.randn(samp_aRs.size)*self.e_Ps[-1] + self.Ps[-1], samp_Ms, 0)
                    v = np.percentile(samp_sma, (16,50,84))
                    self.smas = np.append(self.smas, v[1])
                    self.ehi_smas = np.append(self.ehi_smas, v[2]-v[1])
                    self.elo_smas = np.append(self.elo_smas, v[1]-v[0])
                    samp_Teff = np.random.randn(samp_aRs.size)*d.ehi_Teff + d.Teff
                    samp_Ls = compute_Ls(samp_Rs, samp_Teff)
                    samp_F = compute_F(samp_Ls, samp_sma)
                    v = np.percentile(samp_F, (16,50,84))
                    self.Fs = np.append(self.Fs, v[1])
                    self.ehi_Fs = np.append(self.ehi_Fs, v[2]-v[1])
                    self.elo_Fs = np.append(self.elo_Fs, v[1]-v[0])
                                        
                    #samp_rp = sample_rp(d.object_name, self.prefix, params_opt[3],
                    #                    get_1sigma(params_res[:,3]))
                    #rp,ehi_rp,elo_rp = get_results(samp_rp.reshape(samp_rp.size,
                    #                                               1))
                    #self.rps = np.append(self.rps, rp)
                    #self.ehi_rps = np.append(self.ehi_rps, ehi_rp)
                    #self.elo_rps = np.append(self.elo_rps, elo_rp)
                    #samp_sma = sample_sma(d.object_name, self.prefix, self.Ps[-1],
                    #                      self.e_Ps[-1])
                    #smas = get_results(samp_sma.reshape(samp_sma.size,1))
                    #self.smas = np.append(self.smas, smas[0])
                    #self.ehi_smas = np.append(self.ehi_smas, smas[1])
                    #self.elo_smas = np.append(self.elo_smas, smas[2])
                    #samp_F = compute_F(samp_Ls, samp_sma)
                    #Fs = get_results(samp_F.reshape(samp_F.size,1))
                    #self.Fs = np.append(self.Fs, Fs[0])
                    #self.ehi_Fs = np.append(self.ehi_Fs, Fs[1])
                    #self.elo_Fs = np.append(self.elo_Fs, Fs[2])

                    if j > 0:
                        depth = d.depths[j-1]
                        cdpp = d.CDPPs[j-1]
                        SNRtransit = d.SNRtransits[j-1]
                    else:
                        depth, cdpp, SNRtransit = np.repeat(np.nan, 3)
                    self.depths = np.append(self.depths, depth)
                    self.CDPPs = np.append(self.CDPPs, cdpp)
                    self.SNRtransits = np.append(self.SNRtransits, SNRtransit)
                    
                    # save planet samples for MC sampling
                    if j > 0:
                        samp_P = np.random.normal(self.Ps[-1], self.e_Ps[-1],
                                                  samp_rp.size)
                        outarr = np.array([samp_P, samp_sma, samp_F, samp_rp]).T
                        hdr = 'P_days,sma_AU,F_Fearth,rp_Rearth'
                        fname = 'Gaia-DR2-distances_custom/DistancePosteriors/'
                        fname += '%s_planetpost_%i'%(d.object_name, j)
                        np.savetxt(fname, outarr, header=hdr, delimiter=',',
                                   fmt='%.8e')

        # get confirmed planets from the NASA exoplanet archive
        #self.flag_confirmed_planets_KepK2()
                        
        # save stuff
        _, self.unique_inds = np.unique(self.names_planetsearch,
                                        return_index=True)
        self.Nstars = self.unique_inds.size
        self.Nplanets_detected = self.Ndetected[self.unique_inds].sum()
        self.fs_planetsearch = np.array(self.fs_planetsearch)
	self.SNRtransit_cond = self.cond_vals[:,2]

        # compute map of planet detections via MC simulations
        self.compute_Ndet_map()



    def flag_confirmed_planets_KepK2(self):
        '''Flag planets that are found in the planet search, with confirmed 
        exoplanets from Kepler or K2 based on the NASA exoplanet archive.'''
        # get confirmed planets
        if self.Kep:
            fname = 'input_data/Keplertargets/NASAarchive_confirmed_and_candidates_Kepler_Mdwarf_planets.csv'
            confirmed_IDs, confirmed_Ps = np.loadtxt(fname, delimiter=',', usecols=(1,4),
                                                     skiprows=45).T
        elif self.K2:
            fname = 'input_data/K2targets/NASAarchive_confirmed_and_candidates_K2_Mdwarf_planets.csv'
            epics, Ps = np.genfromtxt(fname, delimiter=',', usecols=(1,5), skip_header=46,
                                      dtype='|S20').T
            confirmed_IDs = np.array([int(s.split(' ')[-1]) for s in epics])
            confirmed_Ps = Ps.astype(float)
        else:
            raise ValueError('Planet search must come from Kepler or K2.')

        # flag confirmed planets in the planetsearch results
        planetsearch_IDs = np.array([int(s.split('_')[-1])
                                     for s in self.names_planetsearch])
        N = planetsearch_IDs.size
        self.confirmed_planet = np.zeros(N) + np.nan
        for i in range(N):
            if self.Ndetected[i] > 0:
                g = confirmed_IDs == planetsearch_IDs[i]
                if g.sum() > 0:  # there is a confirmed planet around this star
                    isconfirmed = np.any(np.isclose(confirmed_Ps[g], self.Ps[i],
                                                    rtol=.05))
                    self.confirmed_planet[i] = float(isconfirmed) if np.isfinite(self.Ps[i]) else np.nan
                else:  # there is no confirmed planet around this star (i.e. FP)
                    self.confirmed_planet[i] = 0.
        
        # revise planet detections
        self.Ndetected_putative = np.copy(self.Ndetected)
        self.Ndetected = something
        
            
    
    def compute_Ndet_map(self, Ntrials=1e3):
        '''Compute the map of planet detections over P and rp using MC 
        simulations over the measurement uncertainties.'''
        if not hasattr(self, 'fs_planetsearch'):
            return None
        N, Ntrials = self.fs_planetsearch.size, int(Ntrials)
        planet_indices = np.zeros(N)
        planets_considered = []
        self.Ps_MC = np.zeros((N, Ntrials))
        self.Fs_MC = np.zeros((N, Ntrials))
        self.as_MC = np.zeros((N, Ntrials))
        self.rps_MC = np.zeros((N, Ntrials))        
        
        # for each detected planet around each star, compute realizations
        # over the uncertainties
        for i in range(N):

            self._conditional_pickleobject(self._t0, time.time())

            if (self.Ndetected[i] > 0) & (np.isfinite(self.Ps[i])):
                
                # compute MC realizations of this planet around this star
                name = self.names_planetsearch[i]
                planet_index = planet_indices[self.names_planetsearch == name].sum() + 1
                planet_indices[i] = 1
                PDF_fname = 'Gaia-DR2-distances_custom/DistancePosteriors/'
                PDF_fname += '%s_planetpost_%i'%(name, planet_index)
                assert PDF_fname not in planets_considered
                planets_considered.append(PDF_fname)
                self.Ps_MC[i],self.Fs_MC[i],self.as_MC[i],self.rps_MC[i] = \
                                    sample_planets(PDF_fname, Ntrials)
                
            else:
                self.Ps_MC[i]  = np.repeat(np.nan, Ntrials)
                self.Fs_MC[i]  = np.repeat(np.nan, Ntrials)
                self.as_MC[i]  = np.repeat(np.nan, Ntrials)
                self.rps_MC[i] = np.repeat(np.nan, Ntrials)

        # compute Ndet map over P and rp for each star
        self.logPgrid = np.logspace(np.log10(self.Plims[0]),
                                    np.log10(self.Plims[1]), self._xlen+1)
        self.logFgrid = np.logspace(np.log10(self.Flims[0]),
                                    np.log10(self.Flims[1]), self._xlen+1)
        self.logagrid = np.logspace(np.log10(self.smalims[0]),
                                    np.log10(self.smalims[1]), self._xlen+1)
        self.logrpgrid = np.logspace(np.log10(self.rplims[0]),
                                     np.log10(self.rplims[1]), self._ylen+1)
        self.NdetP_i = np.zeros((self.Nstars, self._xlen, self._ylen))
        self.NdetF_i = np.zeros((self.Nstars, self._xlen, self._ylen))
        self.Ndeta_i = np.zeros((self.Nstars, self._xlen, self._ylen))
        for i in range(self.Nstars):

            self._conditional_pickleobject(self._t0, time.time())
            name =  self.names_planetsearch[self.unique_inds[i]]
            g1 = self.names_planetsearch == name
            assert g1.sum() > 0

            for j in range(self._xlen):
                for k in range(self._ylen):

                    g2 = (self.Ps_MC[g1] >= self.logPgrid[j]) & \
                         (self.Ps_MC[g1] <= self.logPgrid[j+1]) & \
                         (self.rps_MC[g1] >= self.logrpgrid[k]) & \
                         (self.rps_MC[g1] <= self.logrpgrid[k+1])
                    self.NdetP_i[i,j,k] = g2.sum() / float(Ntrials)

                    g3 = (self.Fs_MC[g1] >= self.logFgrid[j]) & \
                         (self.Fs_MC[g1] <= self.logFgrid[j+1]) & \
                         (self.rps_MC[g1] >= self.logrpgrid[k]) & \
                         (self.rps_MC[g1] <= self.logrpgrid[k+1])
                    self.NdetF_i[i,j,k] = g3.sum() / float(Ntrials)

                    g4 = (self.as_MC[g1] >= self.logagrid[j]) & \
                         (self.as_MC[g1] <= self.logagrid[j+1]) & \
                         (self.rps_MC[g1] >= self.logrpgrid[k]) & \
                         (self.rps_MC[g1] <= self.logrpgrid[k+1])
                    self.Ndeta_i[i,j,k] = g4.sum() / float(Ntrials)
                    
	self.NdetP_tot = fill_map_nans(np.nansum(self.NdetP_i, axis=0))
	self.NdetF_tot = fill_map_nans(np.nansum(self.NdetF_i, axis=0))
	self.Ndeta_tot = fill_map_nans(np.nansum(self.Ndeta_i, axis=0))

        # interpolate onto a fine grid
        xlen, ylen = self._xlen*self._fine_factor, self._ylen*self._fine_factor
        self.logPgrid_fine = np.logspace(np.log10(self.Plims[0]),
                                         np.log10(self.Plims[1]), xlen)
        self.logFgrid_fine = np.logspace(np.log10(self.Flims[0]),
                                         np.log10(self.Flims[1]), xlen)
        self.logagrid_fine = np.logspace(np.log10(self.smalims[0]),
                                         np.log10(self.smalims[1]), xlen)
        self.logrpgrid_fine = np.logspace(np.log10(self.rplims[0]),
                                          np.log10(self.rplims[1]), ylen)
        Pgrid  = np.repeat(self.logPgrid_fine, ylen)
        Fgrid  = np.repeat(self.logFgrid_fine, ylen)
        agrid  = np.repeat(self.logagrid_fine, ylen)
        rpgrid = np.array(list(self.logrpgrid_fine)*xlen)
        self.NdetP_tot_fine_grid = \
                    interpolate_grid(self.logPgrid,self.logrpgrid,
                                     gaussian_filter(self.NdetP_tot,0),
                                     Pgrid,rpgrid).reshape(xlen,ylen)
        self.NdetF_tot_fine_grid = \
                    interpolate_grid(self.logFgrid, self.logrpgrid,
                                     gaussian_filter(self.NdetF_tot,0),
                                     Fgrid, rpgrid).reshape(xlen,ylen)
        self.Ndeta_tot_fine_grid = \
                    interpolate_grid(self.logagrid, self.logrpgrid,
                                     gaussian_filter(self.Ndeta_tot,0),
                                     agrid, rpgrid).reshape(xlen,ylen)


    def save_ids_detectionsfirst(self):
        '''Save a text file of ID numbers ordered by number of detected 
        planet candidates. Also useful for seeing what stars are included in 
        the planet search.'''
        assert self.Ndetected.size == self.names_planetsearch.size
        s = np.argsort(self.Ndetected[self.unique_inds])[::-1]
        IDnums = np.array([int(name.split('_')[-1])
                           for name in self.names_planetsearch[self.unique_inds][s]])
        fname = '%s/%s_detectionsfirst.csv'%(self.folder, self.prefix)
        np.savetxt(fname, IDnums, fmt='%i')
        

    def get_simulation_results(self):
        '''Get data from the simulated planetary systems to compute the 
        sensitivity and FP corrections for each star with a detected planet
        candidate.'''
        # get results from injection/recovery
        fs = glob.glob('%s/%s_*/LC_0*'%(self.folder, self.prefix))
        self.names_simulated = np.unique([i.split('/')[1] for i in fs])
        self.names_simulated = self.names_simulated[self._startstarind:self._endstarind]
        self.Nstars_simulated = self.names_simulated.size
        if self.Nstars_simulated == 0:
            return None

        # stellar params
        self.Nsims = np.zeros(self.Nstars_simulated)
        if self.Kep or self.K2:
            self.Kepmags = np.zeros(self.Nstars_simulated)
        if self.TESS:
            self.TESSmags = np.zeros(self.Nstars_simulated)
        self.efs = np.zeros(self.Nstars_simulated)
        self.Mss = np.zeros(self.Nstars_simulated)
        self.ehi_Mss = np.zeros(self.Nstars_simulated)
        self.elo_Mss = np.zeros(self.Nstars_simulated)
        self.Rss = np.zeros(self.Nstars_simulated)
        self.ehi_Rss = np.zeros(self.Nstars_simulated)
        self.elo_Rss = np.zeros(self.Nstars_simulated)
        self.Teffs = np.zeros(self.Nstars_simulated)
        self.ehi_Teffs = np.zeros(self.Nstars_simulated)
        self.elo_Teffs = np.zeros(self.Nstars_simulated)
        self.loggs = np.zeros(self.Nstars_simulated)
        self.ehi_loggs = np.zeros(self.Nstars_simulated)
        self.elo_loggs = np.zeros(self.Nstars_simulated)
        self.Lss = np.zeros(self.Nstars_simulated)
        self.ehi_Lss = np.zeros(self.Nstars_simulated)
        self.elo_Lss = np.zeros(self.Nstars_simulated)

        Nmaxfs, NmaxPs = 700, 20
        self.Nplanets_inj = np.zeros((self.Nstars_simulated, Nmaxfs)) + np.nan
        self.Nplanets_rec = np.zeros((self.Nstars_simulated, Nmaxfs)) + np.nan
        self.cond_vals_inj = np.zeros((self.Nstars_simulated, Nmaxfs, NmaxPs, 9))+np.nan
        self.cond_free_params_inj = np.zeros((self.Nstars_simulated, Nmaxfs, NmaxPs, 11))+np.nan
        self.Ps_inj = np.zeros((self.Nstars_simulated, Nmaxfs, NmaxPs)) + np.nan
        self.Fs_inj = np.zeros((self.Nstars_simulated, Nmaxfs, NmaxPs)) + np.nan
        self.as_inj = np.zeros((self.Nstars_simulated, Nmaxfs, NmaxPs)) + np.nan
        self.rps_inj = np.zeros((self.Nstars_simulated, Nmaxfs, NmaxPs)) + np.nan
        self.is_rec = np.zeros((self.Nstars_simulated, Nmaxfs, NmaxPs)) + np.nan
        self.Ps_rec = np.zeros((self.Nstars_simulated, Nmaxfs, NmaxPs)) + np.nan
        self.Fs_rec = np.zeros((self.Nstars_simulated, Nmaxfs, NmaxPs)) + np.nan
        self.as_rec = np.zeros((self.Nstars_simulated, Nmaxfs, NmaxPs)) + np.nan
        self.rps_rec = np.zeros((self.Nstars_simulated, Nmaxfs, NmaxPs)) + np.nan
        self.is_FP = np.zeros((self.Nstars_simulated, Nmaxfs, NmaxPs)) + np.nan
        self.depths_inj = np.zeros((self.Nstars_simulated, Nmaxfs, NmaxPs)) + np.nan
        self.CDPPs_inj = np.zeros((self.Nstars_simulated, Nmaxfs, NmaxPs)) + np.nan
        self.Ntransits_inj = np.zeros((self.Nstars_simulated, Nmaxfs, NmaxPs)) + np.nan
        self.SNRtransits_inj = np.zeros((self.Nstars_simulated, Nmaxfs, NmaxPs)) + np.nan
        self.depths_rec = np.zeros((self.Nstars_simulated, Nmaxfs, NmaxPs)) + np.nan
        self.CDPPs_rec = np.zeros((self.Nstars_simulated, Nmaxfs, NmaxPs)) + np.nan
        self.Ntransits_rec = np.zeros((self.Nstars_simulated, Nmaxfs, NmaxPs)) + np.nan
        self.SNRtransits_rec = np.zeros((self.Nstars_simulated, Nmaxfs, NmaxPs)) + np.nan
        
        for i in range(self.Nstars_simulated):

            name = self.names_simulated[i]
            print float(i)/self.Nstars_simulated, name
            self._conditional_pickleobject(self._t0, time.time())
            fs = np.array(glob.glob('%s/%s/LC*'%(self.folder, name)))

	    # remove planet search result (i.e. with index -99)
            g = np.in1d(fs, '%s/%s/LC_-00099'%(self.folder, name))
	    if np.any(g):
	        fs = np.delete(fs, np.where(g)[0][0])
            if fs.size == 0:
	        pass

            # get params of injected and recovered planets for this star
            for j in range(fs.size):

                print float(j) / fs.size
                try:
                    d = loadpickle(fs[j])
                except (EOFError, ValueError):
                    pass

                # save stellar params
                if j == 0:
                    if self.Kep or self.K2:
                        self.Kepmags[i] = d.Kepmag
                    if self.TESS:
                        self.TESSmags[i] = d.TESSmag
		    #self.sigtransits[i] = d.sigtransit
                    self.efs[i] = d.ef.mean()
                    self.Mss[i] = d.Ms
                    self.ehi_Mss[i] = d.ehi_Ms
                    self.elo_Mss[i] = d.elo_Ms
                    self.Rss[i] = d.Rs
                    self.ehi_Rss[i] = d.ehi_Rs
                    self.elo_Rss[i] = d.elo_Rs
                    self.Teffs[i] = d.Teff
                    self.ehi_Teffs[i] = d.ehi_Teff
                    self.elo_Teffs[i] = d.elo_Teff
                    self.loggs[i] = d.logg
                    self.ehi_loggs[i] = d.ehi_logg
                    self.elo_loggs[i] = d.elo_logg
                    _,_,samp_Ls = sample_Ls(d.object_name, self.prefix) 
                    Lss = get_results(samp_Ls.reshape(samp_Ls.size,1))
                    self.Lss[i], self.ehi_Lss, self.elo_Lss = Lss
                
                if d.DONE:
                    self.Nsims[i] += 1
                    self.Nplanets_inj[i,j] = d.Ptrue.size
                    filler = np.repeat(np.nan, NmaxPs-self.Nplanets_inj[i,j])
	            self.Ps_inj[i,j] = np.append(d.Ptrue, filler)
                    sma = rvs.semimajoraxis(d.Ptrue, d.Ms, 0)
                    self.as_inj[i,j] = np.append(sma, filler)
                    F = compute_F(compute_Ls(d.Rs, d.Teff), sma)
                    self.Fs_inj[i,j]  = np.append(F, filler) 
            	    self.rps_inj[i,j] = np.append(d.rptrue, filler)
            	    self.is_rec[i,j]  = np.append(d.is_detected, filler)
                    self.depths_inj[i,j] = np.append(d.depths_inj, filler)
                    self.CDPPs_inj[i,j] = np.append(d.CDPPs_inj, filler)
                    self.Ntransits_inj[i,j] = np.append(d.Ntransits_inj, filler)
                    self.SNRtransits_inj[i,j] = np.append(d.SNRtransits_inj, filler)
                    
                    # get false positives 
            	    params = d.params_guess
		    self.Nplanets_rec[i,j] = params.shape[0]
		    filler2 = np.repeat(np.nan, NmaxPs-self.Nplanets_rec[i,j])
                    rp = rvs.m2Rearth(rvs.Rsun2m(np.sqrt(params[:,2])*d.Rs))
                    self.Ps_rec[i,j] = np.append(params[:,0], filler2)
                    sma = rvs.semimajoraxis(params[:,0], d.Ms, 0)
                    self.as_rec[i,j] = np.append(sma, filler2)
                    F = compute_F(compute_Ls(d.Rs, d.Teff), sma)
                    self.Fs_rec[i,j]  = np.append(F, filler2)
                    self.rps_rec[i,j] = np.append(rp, filler2)
            	    self.is_FP[i,j]   = np.append(d.is_FP, filler2)
                    self.depths_rec[i,j] = np.append(d.depths_rec, filler2)
                    self.CDPPs_rec[i,j] = np.append(d.CDPPs_rec, filler2)
                    self.Ntransits_rec[i,j] = np.append(d.Ntransits_rec, filler2)
                    self.SNRtransits_rec[i,j] = np.append(d.SNRtransits_rec, filler2)
                    
                    # save vetting results for diagnostic purposes
                    NPOIs = d.params_guess_priorto_confirm.shape[0]
                    filler3 = np.repeat(np.nan, 9*(NmaxPs-NPOIs)).reshape(NmaxPs-NPOIs,9)
                    self.cond_vals_inj[i,j] = np.append(d.transit_condition_values, filler3, 0)
                    free_params = np.array(list(d.transit_condition_free_params)*NPOIs).reshape(NPOIs,11)
                    filler4 = np.repeat(np.nan, 11*(NmaxPs-NPOIs)).reshape(NmaxPs-NPOIs,11)
                    self.cond_free_params_inj[i,j] = np.append(free_params, filler4, 0)


        # compute sensitivity and transit probability maps
        self.compute_sens_maps()
        #self.compute_transitprob_maps()

        
        
    def compute_sens_maps(self):
        '''Get all the simulations for all stars and compute the
        sensitivity and the number of FPs as functions of P and rp.'''
        self.logPgrid = np.logspace(np.log10(self.Plims[0]),
                                    np.log10(self.Plims[1]), self._xlen+1)
        self.logFgrid = np.logspace(np.log10(self.Flims[0]),
                                    np.log10(self.Flims[1]), self._xlen+1)
        self.logagrid = np.logspace(np.log10(self.smalims[0]),
                                    np.log10(self.smalims[1]), self._xlen+1)
        self.logrpgrid = np.logspace(np.log10(self.rplims[0]),
                                     np.log10(self.rplims[1]), self._ylen+1)

        self.NrecP_i = np.zeros((self.Nstars_simulated, self._xlen, self._ylen))
        self.NinjP_i = np.zeros((self.Nstars_simulated, self._xlen, self._ylen))
        self.NFPP_i  = np.zeros((self.Nstars_simulated, self._xlen, self._ylen))

        self.NrecF_i = np.zeros((self.Nstars_simulated, self._xlen, self._ylen))
        self.NinjF_i = np.zeros((self.Nstars_simulated, self._xlen, self._ylen))
        self.NFPF_i  = np.zeros((self.Nstars_simulated, self._xlen, self._ylen))
        
        self.Nreca_i = np.zeros((self.Nstars_simulated, self._xlen, self._ylen))
        self.Ninja_i = np.zeros((self.Nstars_simulated, self._xlen, self._ylen))
        self.NFPa_i  = np.zeros((self.Nstars_simulated, self._xlen, self._ylen))

	print 'Compute N_rec and N_inj maps (coarse grid)...'
        for i in range(self.Nstars_simulated):
            self._conditional_pickleobject(self._t0, time.time())
            for j in range(self._xlen):
                for k in range(self._ylen):

                    g = (self.Ps_inj[i] >= self.logPgrid[j]) & \
                        (self.Ps_inj[i] <= self.logPgrid[j+1]) & \
                        (self.rps_inj[i] >= self.logrpgrid[k]) & \
                        (self.rps_inj[i] <= self.logrpgrid[k+1])
                    self.NrecP_i[i,j,k] = self.is_rec[i,g].sum()
                    self.NinjP_i[i,j,k] = self.is_rec[i,g].size

		    g = (self.Ps_rec[i] >= self.logPgrid[j]) & \
                        (self.Ps_rec[i] <= self.logPgrid[j+1]) & \
                        (self.rps_rec[i] >= self.logrpgrid[k]) & \
                        (self.rps_rec[i] <= self.logrpgrid[k+1])
                    self.NFPP_i[i,j,k] = self.is_FP[i,g].sum()

                    g = (self.Fs_inj[i] >= self.logFgrid[j]) & \
                        (self.Fs_inj[i] <= self.logFgrid[j+1]) & \
                        (self.rps_inj[i] >= self.logrpgrid[k]) & \
                        (self.rps_inj[i] <= self.logrpgrid[k+1])
                    self.NrecF_i[i,j,k] = self.is_rec[i,g].sum()
                    self.NinjF_i[i,j,k] = self.is_rec[i,g].size

		    g = (self.Fs_rec[i] >= self.logFgrid[j]) & \
                        (self.Fs_rec[i] <= self.logFgrid[j+1]) & \
                        (self.rps_rec[i] >= self.logrpgrid[k]) & \
                        (self.rps_rec[i] <= self.logrpgrid[k+1])
                    self.NFPF_i[i,j,k] = self.is_FP[i,g].sum()

                    g = (self.as_inj[i] >= self.logagrid[j]) & \
                        (self.as_inj[i] <= self.logagrid[j+1]) & \
                        (self.rps_inj[i] >= self.logrpgrid[k]) & \
                        (self.rps_inj[i] <= self.logrpgrid[k+1])
                    self.Nreca_i[i,j,k] = self.is_rec[i,g].sum()
                    self.Ninja_i[i,j,k] = self.is_rec[i,g].size

		    g = (self.as_rec[i] >= self.logagrid[j]) & \
                        (self.as_rec[i] <= self.logagrid[j+1]) & \
                        (self.rps_rec[i] >= self.logrpgrid[k]) & \
                        (self.rps_rec[i] <= self.logrpgrid[k+1])
                    self.NFPa_i[i,j,k] = self.is_FP[i,g].sum()

        # compute sensitivity
        self.sensP_i = self.NrecP_i / self.NinjP_i.astype(float)
        self.e_sensP_i = np.sqrt(self.NrecP_i) / self.NinjP_i.astype(float)
	self.sensP_avg = fill_map_nans(np.nanmean(self.sensP_i, axis=0))
        self.sensF_i = self.NrecF_i / self.NinjF_i.astype(float)
        self.e_sensF_i = np.sqrt(self.NrecF_i) / self.NinjF_i.astype(float)
	self.sensF_avg = fill_map_nans(np.nanmean(self.sensF_i, axis=0))
        self.sensa_i = self.Nreca_i / self.Ninja_i.astype(float)
        self.e_sensa_i = np.sqrt(self.Nreca_i) / self.Ninja_i.astype(float)
	self.sensa_avg = fill_map_nans(np.nanmean(self.sensa_i, axis=0))

        # compute yield correction to multiply Ndet by
        self.yield_corrP_i = 1 - self.NFPP_i / \
                             (self.NrecP_i + self.NFPP_i.astype(float))
        self.e_yield_corrP_i = np.sqrt(self.NFPP_i) / \
                               (self.NrecP_i + self.NFPP_i.astype(float))
	self.yield_corrP_avg = fill_map_nans(np.nanmean(self.yield_corrP_i,
                                                        axis=0))
        self.yield_corrF_i = 1 - self.NFPF_i / \
                             (self.NrecF_i + self.NFPF_i.astype(float))
        self.e_yield_corrF_i = np.sqrt(self.NFPF_i) / \
                               (self.NrecF_i + self.NFPF_i.astype(float))
	self.yield_corrF_avg = fill_map_nans(np.nanmean(self.yield_corrF_i,
                                                        axis=0))
        self.yield_corra_i = 1 - self.NFPa_i / \
                             (self.Nreca_i + self.NFPa_i.astype(float))
        self.e_yield_corra_i = np.sqrt(self.NFPa_i) / \
                               (self.Nreca_i + self.NFPa_i.astype(float))
	self.yield_corra_avg = fill_map_nans(np.nanmean(self.yield_corra_i,
                                                        axis=0))

        # interpolate onto a fine grid
	print 'Compute N_rec and N_inj maps (fine grid)...'
	return None  # TEMP
        xlen, ylen = self._xlen*self._fine_factor, self._ylen*self._fine_factor
        self.logPgrid_fine = np.logspace(np.log10(self.Plims[0]),
                                         np.log10(self.Plims[1]), xlen)
        self.logFgrid_fine = np.logspace(np.log10(self.Flims[0]),
                                         np.log10(self.Flims[1]), xlen)
        self.logagrid_fine = np.logspace(np.log10(self.smalims[0]),
                                         np.log10(self.smalims[1]), xlen)
        self.logrpgrid_fine = np.logspace(np.log10(self.rplims[0]),
                                          np.log10(self.rplims[1]), ylen)
        Pgrid  = np.repeat(self.logPgrid_fine, ylen)
        Fgrid  = np.repeat(self.logFgrid_fine, ylen)
        agrid  = np.repeat(self.logagrid_fine, ylen)
        rpgrid = np.array(list(self.logrpgrid_fine)*xlen)
        self.sensP_avg_fine_grid = \
                    interpolate_grid(self.logPgrid, self.logrpgrid,
                                     gaussian_filter(self.sensP_avg,.1),
                                     Pgrid, rpgrid).reshape(xlen,ylen)
        self.sensF_avg_fine_grid = \
                    interpolate_grid(self.logFgrid, self.logrpgrid,
                                     gaussian_filter(self.sensF_avg,.1),
                                     Fgrid, rpgrid).reshape(xlen,ylen)
        self.sensa_avg_fine_grid = \
                    interpolate_grid(self.logagrid, self.logrpgrid,
                                     gaussian_filter(self.sensa_avg,.1),
                                     agrid, rpgrid).reshape(xlen,ylen)

        self.yield_corrP_avg_fine_grid = \
                    interpolate_grid(self.logPgrid, self.logrpgrid,
                                     gaussian_filter(self.yield_corrP_avg,.1),
                                     Pgrid, rpgrid).reshape(xlen,ylen)
        self.yield_corrF_avg_fine_grid = \
                    interpolate_grid(self.logFgrid, self.logrpgrid,
                                     gaussian_filter(self.yield_corrF_avg,.1),
                                     Fgrid, rpgrid).reshape(xlen,ylen)
        self.yield_corra_avg_fine_grid = \
                    interpolate_grid(self.logagrid, self.logrpgrid,
                                     gaussian_filter(self.yield_corra_avg,.1),
                                     agrid, rpgrid).reshape(xlen,ylen)



    def compute_transitprob_maps(self):
        '''Compute the transiting probability maps for each star with a
        detected planet candidate.'''
        self.transit_probP_i = np.zeros_like(self.sensP_i)
        self.e_transit_probP_i = np.zeros_like(self.sensP_i)
        self.transit_probF_i = np.zeros_like(self.sensP_i)
        self.e_transit_probF_i = np.zeros_like(self.sensP_i)
        self.transit_proba_i = np.zeros_like(self.sensP_i)
        self.e_transit_proba_i = np.zeros_like(self.sensP_i)

        print 'Compute transit probability maps (coarse grid)...'
        for i in range(self.Nstars_simulated):

            self._conditional_pickleobject(self._t0, time.time())
            # get parameters pdfs
            name = self.names_simulated[i]
            ID = int(name.split('_')[-1])
            fname = 'Gaia-DR2-distances_custom/DistancePosteriors/'
            fname += '%s_allpost_%i'%(self.prefix, ID)
            inds = np.array([9,11])
            samp_Rs,samp_Ms = np.loadtxt(fname, delimiter=',')[:,inds].T
            samp_Rs,_,samp_Ls = sample_Ls(name, self.prefix)
            samp_Ms = resample_PDF(samp_Ms, samp_Rs.size, 1e-3)
            
            for j in range(self._xlen):
                for k in range(self._ylen):

                    # get central values
                    Pmid = 10**(np.log10(self.logPgrid[j]) + \
                                np.diff(np.log10(self.logPgrid[:2])/2))
                    Fmid = 10**(np.log10(self.logFgrid[j]) + \
                                np.diff(np.log10(self.logFgrid[:2])/2))
                    amid = 10**(np.log10(self.logagrid[j]) + \
                                np.diff(np.log10(self.logagrid[:2])/2))
                    rpmid = 10**(np.log10(self.logrpgrid[k]) + \
                                 np.diff(np.log10(self.logrpgrid[:2])/2))

                    # compute transit prob vs P,F,& sma
                    samp_smaP = rvs.AU2m(rvs.semimajoraxis(Pmid, samp_Ms, 0))
                    samp_probP = (rvs.Rsun2m(samp_Rs) + rvs.Rearth2m(rpmid)) / \
                                 samp_smaP
                    probP,_,_ = get_median_results(samp_probP)
                    self.transit_probP_i[i,j,k] = probP
                    #self.e_transit_probP_i[i,j,k] = unp.std_devs(probP)
                    
                    samp_smaF = sma_from_F(Fmid, samp_Ls)
                    samp_probF = (rvs.Rsun2m(samp_Rs) + rvs.Rearth2m(rpmid)) / \
                                 samp_smaF
                    probF,_,_ = get_median_results(samp_probF)
                    self.transit_probF_i[i,j,k] = probF
                    #self.e_transit_probF_i[i,j,k] = unp.std_devs(probF)
                    
                    samp_proba = (rvs.Rsun2m(samp_Rs) + rvs.Rearth2m(rpmid)) / \
                                 amid
                    proba,_,_ = get_median_results(samp_proba)
                    self.transit_proba_i[i,j,k] = proba
                    #self.e_transit_proba_i[i,j,k] = unp.std_devs(proba)


	# correction from beta distribution fit (Kipping 2013)
	self.transit_prob_factor = 1.08
	self.transit_probP_i *= self.transit_prob_factor
	self.e_transit_probP_i *= self.transit_prob_factor
	self.transit_probP_avg = fill_map_nans(np.nanmean(self.transit_probP_i,
                                                          axis=0))
	self.transit_probF_i *= self.transit_prob_factor
	self.e_transit_probF_i *= self.transit_prob_factor
	self.transit_probF_avg = fill_map_nans(np.nanmean(self.transit_probF_i,
                                                          axis=0))
	self.transit_proba_i *= self.transit_prob_factor
	self.e_transit_proba_i *= self.transit_prob_factor
	self.transit_proba_avg = fill_map_nans(np.nanmean(self.transit_proba_i,
                                                          axis=0))

        # interpolate onto a fine grid
        print 'Compute transit probability maps (fine grid)...'
        return None # TEMP
        xlen, ylen = self.logPgrid_fine.size, self.logrpgrid_fine.size
        Pgrid  = np.repeat(self.logPgrid_fine, ylen)
        Fgrid  = np.repeat(self.logFgrid_fine, ylen)
        agrid  = np.repeat(self.logagrid_fine, ylen)
        rpgrid = np.array(list(self.logrpgrid_fine)*xlen)
        self.transit_probP_avg_fine_grid = \
                    interpolate_grid(self.logPgrid, self.logrpgrid,
                                     gaussian_filter(self.transit_probP_avg,0),
                                     Pgrid, rpgrid).reshape(xlen,ylen)
        self.transit_probF_avg_fine_grid = \
                    interpolate_grid(self.logFgrid, self.logrpgrid,
                                     gaussian_filter(self.transit_probF_avg,0),
                                     Fgrid, rpgrid).reshape(xlen,ylen)
        self.transit_proba_avg_fine_grid = \
                    interpolate_grid(self.logagrid, self.logrpgrid,
                                     gaussian_filter(self.transit_proba_avg,0),
                                     agrid, rpgrid).reshape(xlen,ylen)

        
        
    def compute_occurrence_rate(self):
        '''Use the maps of Ndet, yield_corr, sensitivity, and transit 
        probability to compute the occurrence rate over P/F and rp.'''
        assert self.NdetP_i.shape[0] == self.Nstars
        assert self.sensP_i.shape[0] == self.Nstars_simulated
        assert self.NdetP_i.shape[1:] == self.sensP_i.shape[1:]
        assert self.NdetF_i.shape[0] == self.Nstars
        assert self.sensF_i.shape[0] == self.Nstars_simulated
        assert self.NdetF_i.shape[1:] == self.sensF_i.shape[1:]
        assert self.Ndeta_i.shape[0] == self.Nstars
        assert self.sensa_i.shape[0] == self.Nstars_simulated
        assert self.Ndeta_i.shape[1:] == self.sensa_i.shape[1:]
        
        self.occurrence_rateP_i = np.zeros_like(self.sens_i)
        self.occurrence_rateF_i = np.zeros_like(self.sens_i)
        self.occurrence_ratea_i = np.zeros_like(self.sens_i)
        for i in range(self.Nstars_simulated):

	    print float(i) / self.Nstars_simulated
            self._conditional_pickleobject(self._t0, time.time())
            g = self.names_planetsearch[self.unique_inds] == \
                self.names_simulated[i]
            
            NdetP_i = self.NdetP_i[g].reshape(self._xlen, self._ylen)
            SP_i = fill_map_nans(self.sensP_i[i])
            CP_i = fill_map_nans(self.yield_corrP_i[i])
            tP_i = self.transit_probP_i[i]
            self.occurrence_rateP_i[i] = NdetP_i * CP_i / \
                                         (SP_i * tP_i) / self.Nstars

            NdetF_i = self.NdetF_i[g].reshape(self._xlen, self._ylen)
            SF_i = fill_map_nans(self.sensF_i[i])
            CF_i = fill_map_nans(self.yield_corrF_i[i])
            tF_i = self.transit_probF_i[i]
            self.occurrence_rateF_i[i] = NdetF_i * CF_i / \
                                         (SF_i * tF_i) / self.Nstars

            Ndeta_i = self.Ndeta_i[g].reshape(self._xlen, self._ylen)
            Sa_i = fill_map_nans(self.sensa_i[i])
            Ca_i = fill_map_nans(self.yield_corra_i[i])
            ta_i = self.transit_proba_i[i]
            self.occurrence_ratea_i[i] = Ndeta_i * Ca_i / \
                                         (Sa_i * ta_i) / self.Nstars
            
        # compute two verions of the occurrence rate 
        self.occurrence_rateP_v1 = np.nanmean(self.occurrence_rateP_i, axis=0)
        self.occurrence_rateF_v1 = np.nanmean(self.occurrence_rateF_i, axis=0)
        self.occurrence_ratea_v1 = np.nanmean(self.occurrence_ratea_i, axis=0)
        self.occurrence_rateP_v2 = self.NdetP_tot * self.yield_corrP_avg / \
                                   (self.sensP_avg * self.transit_probP_avg) / \
                                   self.Nstars
        self.occurrence_rateF_v2 = self.NdetF_tot * self.yield_corrF_avg / \
                                   (self.sensF_avg * self.transit_probF_avg) / \
                                   self.Nstars
        self.occurrence_ratea_v2 = self.Ndeta_tot * self.yield_corra_avg / \
                                   (self.sensa_avg * self.transit_proba_avg) / \
                                   self.Nstars

        # interpolate onto a fine grid
        xlen, ylen = self.logPgrid_fine.size, self.logrpgrid_fine.size
        Pgrid  = np.repeat(self.logPgrid_fine, ylen)
        Fgrid  = np.repeat(self.logFgrid_fine, ylen)
        agrid  = np.repeat(self.logagrid_fine, ylen)
        rpgrid = np.array(list(self.logrpgrid_fine)*xlen)
        self.occurrence_rateP_fine_grid_v1 = \
                    interpolate_grid(self.logPgrid, self.logrpgrid,
                                     gaussian_filter(self.occurrence_rateP_v1,
                                                     .1), Pgrid,
                                     rpgrid).reshape(xlen,ylen)
        self.occurrence_rateP_fine_grid_v2 = \
                    interpolate_grid(self.logPgrid, self.logrpgrid,
                                     gaussian_filter(self.occurrence_rateP_v2,
                                                     .05), Pgrid,
                                     rpgrid).reshape(xlen,ylen)

        self.occurrence_rateF_fine_grid_v1 = \
                    interpolate_grid(self.logFgrid, self.logrpgrid,
                                     gaussian_filter(self.occurrence_rateF_v1,
                                                     .1), Fgrid,
                                     rpgrid).reshape(xlen,ylen)
        self.occurrence_rateF_fine_grid_v2 = \
                    interpolate_grid(self.logFgrid, self.logrpgrid,
                                     gaussian_filter(self.occurrence_rateF_v2,
                                                     .05), Fgrid,
                                     rpgrid).reshape(xlen,ylen)
        self.occurrence_ratea_fine_grid_v1 = \
                    interpolate_grid(self.logagrid, self.logrpgrid,
                                     gaussian_filter(self.occurrence_ratea_v1,
                                                     .1), agrid,
                                     rpgrid).reshape(xlen,ylen)
        self.occurrence_ratea_fine_grid_v2 = \
                    interpolate_grid(self.logagrid, self.logrpgrid,
                                     gaussian_filter(self.occurrence_ratea_v2,
                                                     .05), agrid,
                                     rpgrid).reshape(xlen,ylen)


    def _conditional_pickleobject(self, t0, tnow):
	'''Pickle object if enough time has passed. Useful for long loops. 
	Input times should be obtained from time.time()'''
        assert tnow > t0
        t_elapsed_hrs = (tnow - t0) / 36e2
	t_tolerance_min = 3
	t_tolerance_hrs = float(t_tolerance_min) / 60
        if np.isclose(t_elapsed_hrs, 24, atol=t_tolerance_hrs) or \
           np.isclose(t_elapsed_hrs, 48, atol=t_tolerance_hrs):
            self._pickleobject()


    def _pickleobject(self):
        fObj = open(self.fname_out, 'wb')
        pickle.dump(self, fObj)
        fObj.close()



def compute_Ls(Rs_Sun, Teff_K):
    return Rs_Sun**2 * (Teff_K / 5772.)**4


def resample_PDF(pdf, Nsamp, sig=1e-3):
    pdf_resamp = np.random.choice(pdf, int(Nsamp)) + np.random.randn(int(Nsamp))*sig
    return pdf_resamp


def sample_Ls(object_name, prefix, Nsamp=1e4):
    '''Sample Ls PDF from the Rs and Teff PDFs for this star.'''
    # get Rs and Teff pdfs
    Nsamp = int(Nsamp)
    try:
        ID = int(object_name.split('_')[-1])
        fname = 'Gaia-DR2-distances_custom/DistancePosteriors/'
        fname += '%s_allpost_%i'%(prefix, ID)
        inds = np.arange(9,11)
        samp_Rs_tmp, samp_Teff_tmp = np.loadtxt(fname, delimiter=',')[:,inds].T
        samp_Rs = resample_PDF(samp_Rs_tmp, Nsamp, 1e-3)
        samp_Teff = resample_PDF(samp_Teff_tmp, Nsamp, 1e-1)
    except IOError:
        samp_Rs = np.repeat(np.nan, Nsamp)
        samp_Teff = np.repeat(np.nan, Nsamp)

    # compute Ls distribution
    return samp_Rs, samp_Teff, samp_Rs**2 * (samp_Teff / 5772.)**4
            

def compute_F(Ls_Sun, smas_AU):
    return Ls_Sun / smas_AU**2


def sma_from_F(F_Sun, samp_Ls_Sun):
    samp_sma = np.sqrt(samp_Ls_Sun / F_Sun)
    return samp_sma


def get_median_results(arr):
    v = np.percentile(arr, (16,50,84))
    return v[1], v[2]-v[1], v[1]-v[0]


def sample_rp(object_name, prefix, rpRs, e_rpRs, Nsamp=1e4):
    '''Sample rp PDF from the Rs PDF for this star.'''
    # get Rs pdf
    Nsamp = int(Nsamp)
    try:
        ID = int(object_name.split('_')[-1])
        fname = 'Gaia-DR2-distances_custom/DistancePosteriors/'
        fname += '%s_allpost_%i'%(prefix, ID)
        samp_Rs_tmp = np.loadtxt(fname, delimiter=',')[:,9]
        samp_Rs = resample_PDF(samp_Rs_tmp, Nsamp, 1e-3)
    except IOError:
        samp_Rs = np.repeat(np.nan, Nsamp)

    # create gaussian rpRs ditribution
    samp_rpRs = np.random.normal(rpRs, e_rpRs, Nsamp)
    
    # compute rp distribution
    return rvs.m2Rearth(rvs.Rsun2m(samp_rpRs*samp_Rs))


def sample_sma(object_name, prefix, P, e_P, Nsamp=1e4):
    '''Sample sma PDF from the Ms PDF for this star.'''
    # get Ms pdf
    Nsamp = int(Nsamp)
    try:
        ID = int(object_name.split('_')[-1])
        fname = 'Gaia-DR2-distances_custom/DistancePosteriors/'
        fname += '%s_allpost_%i'%(prefix, ID)
        samp_Ms_tmp = np.loadtxt(fname, delimiter=',')[:,11]
        samp_Ms = resample_PDF(samp_Ms_tmp, Nsamp, 1e-3)
    except IOError:
        samp_Ms = np.repeat(np.nan, Nsamp)

    # create gaussian P ditribution
    samp_P = np.random.normal(P, e_P, Nsamp)
    
    # compute sma distribution
    return rvs.semimajoraxis(samp_P, samp_Ms, 0)

    
def get_1sigma(results):
    '''results = med, plus_1sig, minus_1sig'''
    assert results.size == 3
    return np.mean(results[1:])


def rpRs2rp(rpRs, Rs):
    rp = rvs.m2Rearth(rvs.Rsun2m(rpRs*Rs))
    return unp.nominal_values(rp), unp.std_devs(rp)


def sample_planets(PDF_fname, N):
    '''Sample planets P, F, and radius from their saved posterior PDFs.'''
    samp_P, samp_sma, samp_F, samp_rp = np.loadtxt(PDF_fname, delimiter=',').T
    N = int(N)
    samp_P = resample_PDF(samp_P, N, 1e-3)
    samp_F = resample_PDF(samp_F, N, 1e-3)
    samp_sma = resample_PDF(samp_sma, N, 1e-3)
    samp_rp = resample_PDF(samp_rp, N, 1e-2)
    return samp_P, samp_F, samp_sma, samp_rp


def fill_map_nans(arr):
    '''Fill nan cells by averging nearby finite cells.'''
    assert len(arr.shape) == 2
    xlen, ylen = arr.shape
    arrv2 = np.copy(arr)

    if np.any(np.isnan(arrv2)):
        # find number of finite cells surrounding nan cells
        while np.any(np.isnan(arrv2)):
            for i in range(xlen):
                for j in range(ylen):
                
                    # index surrounding pixels without going outside the grid
                    surrounding_x_inds_tmp = i+np.delete(range(-1,2)*3, 4)
                    surrounding_y_inds_tmp = j+np.delete(np.append(np.append( \
                                                                np.zeros(3)-1,
                                                                np.zeros(3)),
                                                                np.ones(3)), 4)
                    bad = (surrounding_x_inds_tmp < 0) | \
                          (surrounding_x_inds_tmp >= xlen) | \
                          (surrounding_y_inds_tmp < 0) | \
                          (surrounding_y_inds_tmp >= ylen)
                    surrounding_x_inds = np.delete(surrounding_x_inds_tmp,
                                                   np.where(bad)[0]).astype(int)
                    surrounding_y_inds = np.delete(surrounding_y_inds_tmp,
                                                   np.where(bad)[0]).astype(int)

                    # update the nan cells
                    if np.isnan(arrv2[i,j]):
                        arrv2[i,j] = np.nanmean(arrv2[surrounding_x_inds,
                                                      surrounding_y_inds])

    return arrv2
                        
        

def plot_map(xarr, yarr, zmap, zlabel='', avgtitle=False,
             sumtitle=False, issens=False, pltt=True, label=''):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    img = ax.pcolormesh(xarr, yarr, zmap.T, cmap=plt.get_cmap('hot_r'))
    cbar_axes = fig.add_axes([.1,.1,.87,.04])
    cbar = fig.colorbar(img, cax=cbar_axes, orientation='horizontal')
    cbar.set_label(zlabel)
    if not np.all(np.isclose(np.diff(xarr), np.diff(xarr)[0])):
        ax.set_xscale('log')
    if not np.all(np.isclose(np.diff(yarr), np.diff(yarr)[0])):
        ax.set_yscale('log')

    ax.set_xlabel('Period [days]')
    ax.set_ylabel('Planet Radius [R$_{\oplus}$]')
    if sumtitle:
	ax.set_title('Total = %i'%np.nansum(zmap), fontsize=12)
    if avgtitle:
	ax.set_title('Average = %.3f'%np.nanmean(zmap), fontsize=12)

    # fill nans
    g = np.where(np.isnan(zmap))
    for i in range(g[0].size):
	x1, x2 = xarr[g[0][i]], xarr[g[0][i]+1]
	y1, y2 = yarr[g[1][i]], yarr[g[1][i]+1]
	ax.fill([x1,x2,x2,x1], [y1,y1,y2,y2], fill=False, hatch='\\')

    # fill low sens if plotting a sensitivity map
    if issens:
	g = np.where(zmap<.15)
        for i in range(g[0].size):
            x1, x2 = xarr[g[0][i]], xarr[g[0][i]+1]
            y1, y2 = yarr[g[1][i]], yarr[g[1][i]+1]
            ax.fill([x1,x2,x2,x1], [y1,y1,y2,y2], fill=False, hatch='//')

    fig.subplots_adjust(bottom=.24, top=.95)
    if label != '':
        plt.savefig(label)
    if pltt:
	plt.show()
    plt.close('all')


    
def plot_rp_hist(selfsens, selfdet, pltt=True, label=False):
    '''Plot the 1D planetary radius distribution to visualize the radius gap.'''
    # get 1D rp occurrence rate
    bins_rp = np.logspace(np.log10(.5), np.log10(10), 30)
    Nrec1D = np.zeros(bins_rp.size-1)
    Ninj1D = np.zeros(bins_rp.size-1)
    NFP1D = np.zeros(bins_rp.size-1)
    self = selfsens
    for i in range(bins_rp.size-1):
        g = (self.rps_inj >= bins_rp[i]) & (self.rps_inj <= bins_rp[i+1])
        Nrec1D[i] = self.is_rec[g].sum()
        Ninj1D[i] = self.is_rec[g].size
        g = (self.rps_rec >= bins_rp[i]) & (self.rps_rec <= bins_rp[i+1])
        NFP1D[i] = self.is_FP[g].sum()
    sens1D = Nrec1D / Ninj1D
    yield_corr1D = 1. - NFP1D / (NFP1D + Nrec1D)
    transit_prob = self.transit_probP_avg.mean() # fix prob since the dependence on rp is weak
        
    fig = plt.figure()
    ax = fig.add_subplot(111)
    self = selfdet
    g = np.isfinite(self.rps)
    Ndet, rp_edges = np.histogram(self.rps[g], bins=bins_rp)
    eNdet = np.sqrt(Ndet)
    rparr = 10**(np.log10(rp_edges[:-1]) + np.diff(np.log10(rp_edges))[0]/2)
    occ = rparr*yield_corr1D / sens1D / transit_prob / self.Nstars
    ax.errorbar(rparr, occ, eNdet*0, lw=2, elinewidth=1, capsize=2,
                color='k', drawstyle='steps-mid')
    ax.set_xscale('log')
    if pltt:
        plt.show()
    plt.close('all')



def interpolate_grid(logxarr, logyarr, zarr, xval, yval):
    '''Interpolate over x and y to get the value at the 2d grid z.'''
    assert len(zarr.shape) == 2
    assert logxarr.size == zarr.shape[0]+1
    assert logyarr.size == zarr.shape[1]+1
    
    # create new arrays for interpolation
    logxarrv2 = 10**(np.log10(logxarr[1:]) - np.diff(np.log10(logxarr))[0]/2)
    logxarrv3 = np.repeat(logxarrv2, logyarr.size-1)
    logyarrv2 = 10**(np.log10(logyarr[1:]) - np.diff(np.log10(logyarr))[0]/2)
    logyarrv3 = np.array(list(logyarrv2)*(logxarr.size-1))

    # interpolate
    zarrv2 = zarr.reshape(logxarrv2.size * logyarrv2.size)
    lintz = lint(np.array([logxarrv3, logyarrv3]).T, zarrv2)
    zarrout = lintz(xval, yval)
    return float(zarrout) if zarrout.size == 1 else zarrout


def combine_individual_sens(folder, prefix, label='_0_70'):
    '''sensitivity calculations are done for subsets of all stars to 
    save time. Combine them here.'''
    # get sensitivity objects
    fs = np.sort(np.array(glob.glob('%s/%s_results_*_sens'%(folder, prefix))))
    Nsens = fs.size
    print Nsens
    self = loadpickle(fs[0])
    assert self.Nstars_simulated > 0
    assert label in self.fname_out
    self.fname_out = self.fname_out.replace(label,'')
    
    # append all objects together
    for i in range(1,Nsens):
        
        print float(i) / Nsens
        d = loadpickle(fs[i])
	if d.Nstars_simulated != self.Nstars_simulated:
	    continue

        self.as_inj = np.append(self.as_inj, d.as_inj, 0)
        self.as_rec = np.append(self.as_rec, d.as_rec, 0)
        self.SNRtransits_rec = np.append(self.SNRtransits_rec, d.SNRtransits_rec, 0)
	self.SNRtransits_inj = np.append(self.SNRtransits_inj, d.SNRtransits_inj, 0)
        self.cond_free_params_inj = np.append(self.cond_free_params_inj,
                                              d.cond_free_params_inj, 0)
        self.cond_vals_inj = np.append(self.cond_vals_inj, 
				       d.cond_vals_inj, 0)
        self.e_sensa_i = np.append(self.e_sensa_i, d.e_sensa_i, 0)
        self.e_sensF_i = np.append(self.e_sensF_i, d.e_sensF_i, 0)
        self.e_sensP_i = np.append(self.e_sensP_i, d.e_sensP_i, 0)
        self.e_yield_corra_i = np.append(self.e_yield_corra_i,
                                         d.e_yield_corra_i, 0)
        self.e_yield_corrF_i = np.append(self.e_yield_corrF_i,
                                         d.e_yield_corrF_i, 0)
        self.e_yield_corrP_i = np.append(self.e_yield_corrP_i,
                                         d.e_yield_corrP_i, 0)
        self.efs = np.append(self.efs, d.efs)
        self.ehi_loggs = np.append(self.ehi_loggs, d.ehi_loggs)
        self.ehi_Lss = np.append(self.ehi_Lss, d.ehi_Lss)
        self.ehi_Mss = np.append(self.ehi_Mss, d.ehi_Mss)
        self.ehi_Rss = np.append(self.ehi_Rss, d.ehi_Rss)
        self.ehi_Teffs = np.append(self.ehi_Teffs, d.ehi_Teffs)
        self.elo_loggs = np.append(self.elo_loggs, d.elo_loggs)
        self.elo_Lss = np.append(self.elo_Lss, d.elo_Lss)
        self.elo_Mss = np.append(self.elo_Mss, d.elo_Mss)
        self.elo_Rss = np.append(self.elo_Rss, d.elo_Rss)
        self.elo_Teffs = np.append(self.elo_Teffs, d.elo_Teffs)
        self.Fs_inj = np.append(self.Fs_inj, d.Fs_inj, 0)
        self.Fs_rec = np.append(self.Fs_rec, d.Fs_rec, 0)
        self.is_FP = np.append(self.is_FP, d.is_FP, 0)
        self.is_rec = np.append(self.is_rec, d.is_rec, 0)
        if self.Kep or self.K2:
            self.Kepmags = np.append(self.Kepmags, d.Kepmags)
	if self.TESS:
            self.TESSmags = np.append(self.TESSmags, d.TESSmags)
	#self.sigtransits = np.append(self.sigtransits, d.sigtransits)
        self.Lss = np.append(self.Lss, d.Lss)
        self.loggs = np.append(self.loggs, d.loggs)
        self.Mss = np.append(self.Mss, d.Mss)
        self.names_simulated = np.append(self.names_simulated,
                                         d.names_simulated)
        self.NFPa_i = np.append(self.NFPa_i, d.NFPa_i, 0)
        self.NFPF_i = np.append(self.NFPF_i, d.NFPF_i, 0)
        self.NFPP_i = np.append(self.NFPP_i, d.NFPP_i, 0)
        self.Ninja_i = np.append(self.Ninja_i, d.Ninja_i, 0)
        self.NinjF_i = np.append(self.NinjF_i, d.NinjF_i, 0)
        self.NinjP_i = np.append(self.NinjP_i, d.NinjP_i, 0)
        self.Nplanets_inj = np.append(self.Nplanets_inj, d.Nplanets_inj, 0)
        self.Nplanets_rec = np.append(self.Nplanets_rec, d.Nplanets_rec, 0)
        self.Nreca_i = np.append(self.Nreca_i, d.Nreca_i, 0)
        self.NrecF_i = np.append(self.NrecF_i, d.NrecF_i, 0)
        self.NrecP_i = np.append(self.NrecP_i, d.NrecP_i, 0)
        self.Nsims = np.append(self.Nsims, d.Nsims)
        self.Nstars_simulated += d.Nstars_simulated
        self.Ps_inj = np.append(self.Ps_inj, d.Ps_inj, 0)
        self.Ps_rec = np.append(self.Ps_rec, d.Ps_rec, 0)
        self.Rss = np.append(self.Rss, d.Rss)      
        self.rps_inj = np.append(self.rps_inj, d.rps_inj, 0)
        self.rps_rec = np.append(self.rps_rec, d.rps_rec, 0)
        self.sensa_i = np.append(self.sensa_i, d.sensa_i, 0)
        self.sensF_i = np.append(self.sensF_i, d.sensF_i, 0)
        self.sensP_i = np.append(self.sensP_i, d.sensP_i, 0)
        self.Teffs = np.append(self.Teffs, d.Teffs)
        self.yield_corra_i = np.append(self.yield_corra_i, d.yield_corra_i, 0)
        self.yield_corrF_i = np.append(self.yield_corrF_i, d.yield_corrF_i, 0)
        self.yield_corrP_i = np.append(self.yield_corrP_i, d.yield_corrP_i, 0)
        
    # trim excess planets
    endfs = int(np.nanmax(self.Nsims))
    endPs = int(np.nanmax(self.Nplanets_inj))
    self.cond_vals_inj = self.cond_vals_inj[:,:endfs,:endPs,:]
    self.cond_free_params_inj = self.cond_free_params_inj[:,:endfs,:endPs,:]
    self.SNRtransits_inj = self.SNRtransits_inj[:,:endfs,:endPs]
    self.Ps_inj = self.Ps_inj[:,:endfs,:endPs]
    self.Fs_inj = self.Fs_inj[:,:endfs,:endPs]
    self.as_inj = self.as_inj[:,:endfs,:endPs]
    self.rps_inj = self.rps_inj[:,:endfs,:endPs]
    self.is_rec = self.is_rec[:,:endfs,:endPs]
    endPs = int(np.nanmax(self.Nplanets_rec))
    self.Ps_rec = self.Ps_rec[:,:endfs,:endPs]
    self.Fs_rec = self.Fs_rec[:,:endfs,:endPs]
    self.as_rec = self.as_rec[:,:endfs,:endPs]
    self.SNRtransits_rec  = self.SNRtransits_rec[:,:endfs,:endPs]
    self.rps_rec = self.rps_rec[:,:endfs,:endPs]
    self.is_FP = self.is_FP[:,:endfs,:endPs]
    
    self._pickleobject()

        

if __name__ == '__main__':
    folder = sys.argv[1]  # 'PipelineResults_KepID'
    prefix = sys.argv[2]  # 'KepID'
    startstarind = int(sys.argv[3])  # 0
    endstarind = int(sys.argv[4])    # 10
    #self = OccurrenceRateclass(folder, prefix, startstarind, endstarind,
    #                           compute_detections=False,
    #                           compute_sens=True,
    #                           compute_occurrence_rate=False)
    combine_individual_sens(folder, prefix, '_0_81')
