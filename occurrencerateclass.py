from LCclass import *
from K2sensclass import *
from uncertainties import unumpy as unp
import rvs
from scipy.ndimage.filters import gaussian_filter # for map smoothing if desired
from scipy.interpolate import LinearNDInterpolator as lint


class OccurrenceRateclass:

    def __init__(self, folder, prefix, xlen=20, ylen=12, 
		 compute_detections=False, compute_sens=False, 
		 compute_occurrence_rate=False,
                 Plims=(.5,80), Flims=(.1,4e2), rplims=(.5,10)):
        self.folder, self.prefix = folder, prefix
	self.fname_out = '%s/%s_K2results'%(self.folder, self.prefix)
        self._xlen, self._ylen = int(xlen), int(ylen)
        self.Plims, self.Flims, self.rplims = Plims, Flims, rplims
        
        if compute_detections:
            self.get_planetsearch_results()
            self._pickleobject()
            
        if compute_sens:
            self.get_simulation_results()
            self._pickleobject()

        if compute_occurrence_rate:
            self.compute_occurrence_rate()
	    self._pickleobject()


    def get_planetsearch_results(self):
        '''Get the results from the planetsearch, i.e. the detected planets 
        and stellar properties.'''
        fs = np.array(glob.glob('%s/%s_*/LC_-00099'%(self.folder, self.prefix))
        if fs.size == 0:
            return None
	self.fs_planetsearch, self.epicnums_planetsearch = [], np.zeros(0)

        # detected planet params
        self.Ndetected, self.params_guess = np.zeros(0), np.zeros((0,4))
        self.params_optimized = np.zeros((0,5))
        self.Ps, self.e_Ps = np.zeros(0), np.zeros(0)
        self.rps, self.e_rps = np.zeros(0), np.zeros(0)        

        # POI params
        self.cond_vals, self.cond_free_params = np.zeros((0,7)), np.zeros((0,7))

        # stellar params
        self.Kepmags, self.efs = np.zeros(0), np.zeros(0)
        self.Mss, self.e_Mss = np.zeros(0), np.zeros(0)
        self.Rss, self.e_Rss = np.zeros(0), np.zeros(0)
        self.Teffs, self.e_Teffs = np.zeros(0), np.zeros(0)
        self.loggs, self.e_loggs = np.zeros(0), np.zeros(0)

        for i in range(fs.size):

            print float(i)/fs.size, fs[i]
            d = loadpickle(fs[i])
	    if d.DONE:
                
		for j in range(d.Ndet+1):

                    self.fs_planetsearch.append(fs[i])
                    self.epicnums_planetsearch = \
                                        np.append(self.epicnums_planetsearch,
                                                  d.epicnum)
                    self.Ndetected = np.append(self.Ndetected, d.Ndet)
                    self.Kepmags = np.append(self.Kepmags, d.Kepmag)
                    self.efs = np.append(self.efs, d.ef.mean())
                    self.Mss = np.append(self.Mss, d.Ms)
                    self.e_Mss = np.append(self.e_Mss, d.e_Ms)
                    self.Rss = np.append(self.Rss, d.Rs)
                    self.e_Rss = np.append(self.e_Rss, d.e_Rs)
                    self.Teffs = np.append(self.Teffs, d.Teff)
                    self.e_Teffs = np.append(self.e_Teffs, d.e_Teff)
                    self.loggs = np.append(self.loggs, d.logg)
                    self.e_loggs = np.append(self.e_loggs, d.e_logg)

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
                                    else np.repeat(np.nan,7)
                    else:
                        cond_vals = np.repeat(np.nan, 7)
		    self.cond_vals = np.append(self.cond_vals,
                                               cond_vals.reshape(1,7), axis=0)
                    self.cond_free_params = np.append(self.cond_free_params,
                                d.transit_condition_free_params.reshape(1,7),
                                                      axis=0)

                    # save periods and planets radii
                    self.Ps = np.append(self.Ps, params_opt[0])
                    self.e_Ps = np.append(self.e_Ps,
                                          get_1sigma(params_res[:,0]))
                    rpRs = unp.uarray(params_opt[3],
                                      get_1sigma(params_res[:,3]))
                    Rs = unp.uarray(d.Rs, d.e_Rs)
                    rp, e_rp = rpRs2rp(rpRs, Rs)
                    self.rps = np.append(self.rps, rp)
                    self.e_rps = np.append(self.e_rps, e_rp)

        # save stellar luminosities and planet sma and F
        Rss   = unp.uarray(self.Rss, self.e_Rss)
        Teffs = unp.uarray(self.Teffs, self.e_Teffs)
        Lss = compute_Ls(Rss, Teffs)
        self.Lss, self.e_Lss = unp.nominal_values(Lss), unp.std_devs(Lss)
        Ps = unp.uarray(self.Ps, self.e_Ps)
        Mss   = unp.uarray(self.Mss, self.e_Mss)
        smas = rvs.semimajoraxis(Ps, Mss, 0)
        self.smas, self.e_smas = unp.nominal_values(smas), unp.std_devs(smas)
        Fs = compute_F(Lss, smas)
        self.Fs, self.e_Fs = unp.nominal_values(Fs), unp.std_devs(Fs)
        
        # save stuff
        _, self.unique_inds = np.unique(self.epicnums_planetsearch,
                                        return_index=True)
        self.Nstars = self.unique_inds.size
        self.Nplanets_detected = self.Ndetected[self.unique_inds].sum()
        self.fs_planetsearch = np.array(self.fs_planetsearch)

        # compute map of planet detections via MC simulations
        self.compute_Ndet_map()

        
    def compute_Ndet_map(self, Ntrials=1e3):
        '''Compute the map of planet detections over P and rp using MC 
        simulations over the measurement uncertainties.'''
        if not hasattr(self, 'fs_planetsearch'):
            return None
        N, Ntrials = self.fs_planetsearch.size, int(Ntrials)
        self.Ps_MC = np.zeros((N, Ntrials))
        self.Fs_MC = np.zeros((N, Ntrials))
        self.rps_MC = np.zeros((N, Ntrials))        
        
        # for each detected planet around each star, compute realizations
        # over the uncertainties
        for i in range(N):

            if self.Ndetected[i] > 0:
                
                # compute MC realizations of this planet
                P, eP = self.Ps[i], self.e_Ps[i]
                rp, erp = self.rps[i], self.e_rps[i]
                Ms, eMs = self.Mss[i], self.e_Mss[i]
                Ls, eLs = self.Lss[i], self.e_Lss[i]
                self.Ps_MC[i],self.Fs_MC[i],self.rps_MC[i] = \
                                                    sample_planets(P, eP,
                                                                   rp, erp,
                                                                   Ms, eMs,
                                                                   Ls, eLs,
                                                                   Ntrials)
            else:
                self.Ps_MC[i]  = np.repeat(np.nan, Ntrials)
                self.Fs_MC[i]  = np.repeat(np.nan, Ntrials)
                self.rps_MC[i] = np.repeat(np.nan, Ntrials)

        # compute Ndet map over P and rp for each star
        self.logPgrid = np.logspace(np.log10(self.Plims[0]),
                                    np.log10(self.Plims[1]), self._xlen+1)
        self.logFgrid = np.logspace(np.log10(self.Flims[0]),
                                    np.log10(self.Flims[1]), self._xlen+1)
        self.logrpgrid = np.logspace(np.log10(self.rplims[0]),
                                     np.log10(self.rplims[1]), self._ylen+1)
        self.NdetP_i = np.zeros((self.Nstars, self._xlen, self._ylen))
        self.NdetF_i = np.zeros((self.Nstars, self._xlen, self._ylen))
        for i in range(self.Nstars):

            epicnum =  self.epicnums_planetsearch[self.unique_inds[i]]
            g1 = self.epicnums_planetsearch == epicnum
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
                    
	self.NdetP_tot = fill_map_nans(np.nansum(self.NdetP_i, axis=0))
	self.NdetF_tot = fill_map_nans(np.nansum(self.NdetF_i, axis=0))

        # interpolate onto a fine grid
        xlen, ylen = self._xlen*6, self._ylen*6
        self.logPgrid_fine = np.logspace(np.log10(self.Plims[0]),
                                         np.log10(self.Plims[1]), xlen)
        self.logFgrid_fine = np.logspace(np.log10(self.Flims[0]),
                                         np.log10(self.Flims[1]), xlen)
        self.logrpgrid_fine = np.logspace(np.log10(self.rplims[0]),
                                          np.log10(self.rplims[1]), ylen)
        Pgrid  = np.repeat(self.logPgrid_fine, ylen)
        Fgrid  = np.repeat(self.logFgrid_fine, ylen)
        rpgrid = np.array(list(self.logrpgrid_fine)*xlen)
        self.NdetP_tot_fine_grid = \
                    interpolate_grid(self.logPgrid,self.logrpgrid,
                                     gaussian_filter(self.NdetP_tot,0),
                                     Pgrid,rpgrid).reshape(xlen,ylen)
        self.NdetF_tot_fine_grid = \
                    interpolate_grid(self.logFgrid, self.logrpgrid,
                                     gaussian_filter(self.NdetF_tot,0),
                                     Fgrid, rpgrid).reshape(xlen,ylen)
                                                    

    def save_stars_with_detections(self):
        epic_tmp = self.epicnums_planetsearch[self.unique_inds]
        epicnum2save = epic_tmp[self.Ndetected[self.unique_inds] > 0]
	self.Nstars_wdet = epicnum2save.size
        np.savetxt('input_data/K2targets/K2Mdwarfs_withdetections.csv',
                   epicnum2save, delimiter=',', fmt='%i')


    def get_simulation_results(self):
        '''Get data from the simulated planetary systems to compute the 
        sensitivity and FP corrections for each star with a detected planet
        candidate.'''
        # get results from injection/recovery
        self.epicnums_wdet = np.loadtxt('input_data/K2targets/' + \
                                        'K2Mdwarfs_withdetections.csv',
                                        delimiter=',')
        self.Nstars_wdet = self.epicnums_wdet.size
        self.Nsims = np.zeros(self.Nstars_wdet)
        Nmaxfs, NmaxPs = 700, 20
        self.Nplanets_inj = np.zeros((self.Nstars_wdet, Nmaxfs)) + np.nan
        self.Nplanets_rec = np.zeros((self.Nstars_wdet, Nmaxfs)) + np.nan
        self.Ps_inj = np.zeros((self.Nstars_wdet, Nmaxfs, NmaxPs)) + np.nan
        self.Fs_inj = np.zeros((self.Nstars_wdet, Nmaxfs, NmaxPs)) + np.nan
        self.rps_inj = np.zeros((self.Nstars_wdet, Nmaxfs, NmaxPs)) + np.nan
        self.is_rec = np.zeros((self.Nstars_wdet, Nmaxfs, NmaxPs)) + np.nan
        self.Ps_rec = np.zeros((self.Nstars_wdet, Nmaxfs, NmaxPs)) + np.nan
        self.Fs_rec = np.zeros((self.Nstars_wdet, Nmaxfs, NmaxPs)) + np.nan
        self.rps_rec = np.zeros((self.Nstars_wdet, Nmaxfs, NmaxPs)) + np.nan
        self.is_FP = np.zeros((self.Nstars_wdet, Nmaxfs, NmaxPs)) + np.nan
        
        for i in range(self.Nstars_wdet):

            epicnum = self.epicnums_wdet[i]
            print i, epicnum
            fs = np.array(glob.glob('%s/%s_%i/LC*'%(self.folder, self.prefix, epicnum)))

	    # remove planet search result (i.e. with index -99)
            g = np.in1d(fs, '%s/%s_%i/LC_-00099'%(self.folder, self.prefix,
						  epicnum))
	    if np.any(g):
	        fs = np.delete(fs, np.where(g)[0][0])
            if fs.size == 0:
	        pass

            # get params of injected and recovered planets for this star
            for j in range(fs.size):

                print float(j) / fs.size
                d = loadpickle(fs[j])
                if d.DONE:
                    self.Nsims[i] += 1
                    self.Nplanets_inj[i,j] = d.Ptrue.size
                    filler = np.repeat(np.nan, NmaxPs-self.Nplanets_inj[i,j])
	            self.Ps_inj[i,j]  = np.append(d.Ptrue, filler)
                    F = compute_F(compute_Ls(d.Rs, d.Teff),
                                  rvs.semimajoraxis(d.Ptrue, d.Ms, 0))
                    self.Fs_inj[i,j]  = np.append(F, filler) 
            	    self.rps_inj[i,j] = np.append(d.rptrue, filler)
            	    self.is_rec[i,j]  = np.append(d.is_detected, filler)

                    # get false positives 
            	    params = d.params_guess
		    self.Nplanets_rec[i,j] = params.shape[0]
		    filler2 = np.repeat(np.nan, NmaxPs-self.Nplanets_rec[i,j])
                    rp = rvs.m2Rearth(rvs.Rsun2m(np.sqrt(params[:,2])*d.Rs))
                    self.Ps_rec[i,j]  = np.append(params[:,0], filler2)
                    F = compute_F(compute_Ls(d.Rs, d.Teff),
                                  rvs.semimajoraxis(params[:,0], d.Ms, 0))
                    self.Fs_rec[i,j]  = np.append(F, filler2)
                    self.rps_rec[i,j] = np.append(rp, filler2)
            	    self.is_FP[i,j]   = np.append(d.is_FP, filler2)

        # trim excess planets
        endfs = int(np.nanmax(self.Nsims))
        endPs = int(np.nanmax(self.Nplanets_inj))
        self.Ps_inj = self.Ps_inj[:,:endfs,:endPs]
        self.Fs_inj = self.Fs_inj[:,:endfs,:endPs]
        self.rps_inj = self.rps_inj[:,:endfs,:endPs]
        self.is_rec = self.is_rec[:,:endfs,:endPs]

        endPs = int(np.nanmax(self.Nplanets_rec))
        self.Ps_rec = self.Ps_rec[:,:endfs,:endPs]
        self.Fs_rec = self.Fs_rec[:,:endfs,:endPs]
        self.rps_rec = self.rps_rec[:,:endfs,:endPs]
        self.is_FP = self.is_FP[:,:endfs,:endPs]

        # compute sensitivity and transit probability maps
        self.compute_sens_maps()
        self.compute_transitprob_maps()

        
        
    def compute_sens_maps(self):
        '''Get all the simulations for all stars and compute the
        sensitivity and the number of FPs as functions of P and rp.'''
        self.logPgrid = np.logspace(np.log10(self.Plims[0]),
                                    np.log10(self.Plims[1]), self._xlen+1)
        self.logFgrid = np.logspace(np.log10(self.Flims[0]),
                                    np.log10(self.Flims[1]), self._xlen+1)
        self.logrpgrid = np.logspace(np.log10(self.rplims[0]),
                                     np.log10(self.rplims[1]), self._ylen+1)
        self.NrecP_i = np.zeros((self.Nstars_wdet, self._xlen, self._ylen))
        self.NinjP_i = np.zeros((self.Nstars_wdet, self._xlen, self._ylen))
        self.NFPP_i  = np.zeros((self.Nstars_wdet, self._xlen, self._ylen))

        self.NrecF_i = np.zeros((self.Nstars_wdet, self._xlen, self._ylen))
        self.NinjF_i = np.zeros((self.Nstars_wdet, self._xlen, self._ylen))
        self.NFPF_i  = np.zeros((self.Nstars_wdet, self._xlen, self._ylen))

        for i in range(self.Nstars_wdet):
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

        # compute sensitivity
        self.sensP_i = self.NrecP_i / self.NinjP_i.astype(float)
        self.e_sensP_i = np.sqrt(self.NrecP_i) / self.NinjP_i.astype(float)
	self.sensP_avg = fill_map_nans(np.nanmean(self.sensP_i, axis=0))
        self.sensF_i = self.NrecF_i / self.NinjF_i.astype(float)
        self.e_sensF_i = np.sqrt(self.NrecF_i) / self.NinjF_i.astype(float)
	self.sensF_avg = fill_map_nans(np.nanmean(self.sensF_i, axis=0))

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

        # interpolate onto a fine grid
        xlen, ylen = self.logPgrid_fine.size, self.logrpgrid_fine.size
        Pgrid  = np.repeat(self.logPgrid_fine, ylen)
        Fgrid  = np.repeat(self.logFgrid_fine, ylen)
        rpgrid = np.array(list(self.logrpgrid_fine)*xlen)
        self.sensP_avg_fine_grid = \
                    interpolate_grid(self.logPgrid, self.logrpgrid,
                                     gaussian_filter(self.sensP_avg,.1),
                                     Pgrid, rpgrid).reshape(xlen,ylen)
        self.sensF_avg_fine_grid = \
                    interpolate_grid(self.logFgrid, self.logrpgrid,
                                     gaussian_filter(self.sensF_avg,.1),
                                     Fgrid, rpgrid).reshape(xlen,ylen)

        self.yield_corrP_avg_fine_grid = \
                    interpolate_grid(self.logPgrid, self.logrpgrid,
                                     gaussian_filter(self.yield_corrP_avg,.1),
                                     Pgrid, rpgrid).reshape(xlen,ylen)
        self.yield_corrF_avg_fine_grid = \
                    interpolate_grid(self.logFgrid, self.logrpgrid,
                                     gaussian_filter(self.yield_corrF_avg,.1),
                                     Fgrid, rpgrid).reshape(xlen,ylen)

        

    def compute_transitprob_maps(self):
        '''Compute the transiting probability maps for each star with a
        detected planet candidate.'''
        self.transit_probP_i = np.zeros_like(self.sensP_i)
        self.e_transit_probP_i = np.zeros_like(self.sensP_i)
        self.transit_probF_i = np.zeros_like(self.sensP_i)
        self.e_transit_probF_i = np.zeros_like(self.sensP_i)

        for i in range(self.Nstars_wdet):
            for j in range(self._xlen):
                for k in range(self._ylen):
                    
                    Pmid = 10**(np.log10(self.logPgrid[j]) + \
                                np.diff(np.log10(self.logPgrid[:2])/2))
                    Fmid = 10**(np.log10(self.logFgrid[j]) + \
                                np.diff(np.log10(self.logFgrid[:2])/2))
                    rpmid = 10**(np.log10(self.logrpgrid[k]) + \
                                 np.diff(np.log10(self.logrpgrid[:2])/2))
                    
                    epicnum = self.epicnums_wdet[i]
                    g = np.where(self.epicnums_planetsearch == epicnum)[0][0]
                    Ms = unp.uarray(self.Mss[g], self.e_Mss[g])
                    smaP = rvs.AU2m(rvs.semimajoraxis(Pmid, Ms, 0))
                    Ls = unp.uarray(self.Lss[g], self.e_Lss[g])
                    smaF = sma_from_F(Fmid, Ls)
                    Rs = unp.uarray(self.Rss[g], self.e_Rss[g])

                    probP = (rvs.Rsun2m(Rs) + rvs.Rearth2m(rpmid)) / smaP
                    self.transit_probP_i[i,j,k] = unp.nominal_values(probP)
                    self.e_transit_probP_i[i,j,k] = unp.std_devs(probP)

                    probF = (rvs.Rsun2m(Rs) + rvs.Rearth2m(rpmid)) / smaF
                    self.transit_probF_i[i,j,k] = unp.nominal_values(probF)
                    self.e_transit_probF_i[i,j,k] = unp.std_devs(probF)
                    

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

        # interpolate onto a fine grid
        xlen, ylen = self.logPgrid_fine.size, self.logrpgrid_fine.size
        Pgrid  = np.repeat(self.logPgrid_fine, ylen)
        Fgrid  = np.repeat(self.logFgrid_fine, ylen)
        rpgrid = np.array(list(self.logrpgrid_fine)*xlen)
        self.transit_probP_avg_fine_grid = \
                    interpolate_grid(self.logPgrid, self.logrpgrid,
                                     gaussian_filter(self.transit_probP_avg,0),
                                     Pgrid, rpgrid).reshape(xlen,ylen)
        self.transit_probF_avg_fine_grid = \
                    interpolate_grid(self.logFgrid, self.logrpgrid,
                                     gaussian_filter(self.transit_probF_avg,0),
                                     Fgrid, rpgrid).reshape(xlen,ylen)

        
        
    def compute_occurrence_rate(self):
        '''Use the maps of Ndet, yield_corr, sensitivity, and transit 
        probability to compute the occurrence rate over P/F and rp.'''
        assert self.NdetP_i.shape[0] == self.Nstars
        assert self.sensP_i.shape[0] == self.Nstars_wdet
        assert self.NdetP_i.shape[1:] == self.sensP_i.shape[1:]
        assert self.NdetF_i.shape[0] == self.Nstars
        assert self.sensF_i.shape[0] == self.Nstars_wdet
        assert self.NdetF_i.shape[1:] == self.sensF_i.shape[1:]
        
        self.occurrence_rateP_i = np.zeros_like(self.sens_i)
        self.occurrence_rateF_i = np.zeros_like(self.sens_i)
        for i in range(self.Nstars_wdet):

	    print float(i) / self.Nstars_wdet
            g = self.epicnums_planetsearch[self.unique_inds] == \
                self.epicnums_wdet[i]
            
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
            
        # compute two verions of the occurrence rate 
        self.occurrence_rateP_v1 = np.nanmean(self.occurrence_rateP_i, axis=0)
        self.occurrence_rateF_v1 = np.nanmean(self.occurrence_rateF_i, axis=0)
        self.occurrence_rateP_v2 = self.NdetP_tot * self.yield_corrP_avg / \
                                   (self.sensP_avg * self.transit_probP_avg) / \
                                   self.Nstars
        self.occurrence_rateF_v2 = self.NdetF_tot * self.yield_corrF_avg / \
                                   (self.sensF_avg * self.transit_probF_avg) / \
                                   self.Nstars

        # interpolate onto a fine grid
        xlen, ylen = self.logPgrid_fine.size, self.logrpgrid_fine.size
        Pgrid  = np.repeat(self.logPgrid_fine, ylen)
        Fgrid  = np.repeat(self.logFgrid_fine, ylen)
        rpgrid = np.array(list(self.logrpgrid_fine)*xlen)
        self.occurrence_rateP_fine_grid_v1 = \
                    interpolate_grid(self.logPgrid, self.logrpgrid,
                                     gaussian_filter(self.occurrence_rateP_v1,
                                                     .1), Pgrid,
                                     rpgrid).reshape(xlen,ylen)
        self.transit_probF_fine_grid_v1 = \
                    interpolate_grid(self.logFgrid, self.logrpgrid,
                                     gaussian_filter(self.occurrence_rateF_v1,
                                                     .1), Fgrid,
                                     rpgrid).reshape(xlen,ylen)
        self.occurrence_rateP_fine_grid_v2 = \
                    interpolate_grid(self.logPgrid, self.logrpgrid,
                                     gaussian_filter(self.occurrence_rateP_v2,
                                                     .05), Pgrid,
                                     rpgrid).reshape(xlen,ylen)
        self.transit_probF_fine_grid_v2 = \
                    interpolate_grid(self.logFgrid, self.logrpgrid,
                                     gaussian_filter(self.occurrence_rateF_v2,
                                                     .05), Fgrid,
                                     rpgrid).reshape(xlen,ylen)

        


    def _pickleobject(self):
        fObj = open(self.fname_out, 'wb')
        pickle.dump(self, fObj)
        fObj.close()



def compute_Ls(Rs_Sun, Teff_K):
    return Rs_Sun**2 * (Teff_K / 5772.)**4


def compute_F(Ls_Sun, smas_AU):
    return Ls_Sun / smas_AU**2


def sma_from_F(F_Sun, Ls_Sun):
    return unp.sqrt(Ls_Sun / F_sun)
    

def get_1sigma(results):
    '''results = med, plus_1sig, minus_1sig'''
    assert results.size == 3
    return np.mean(results[1:])


def rpRs2rp(rpRs, Rs):
    rp = rvs.m2Rearth(rvs.Rsun2m(rpRs*Rs))
    return unp.nominal_values(rp), unp.std_devs(rp)


def sample_planets(P, eP, rp, erp, Ms, eMs, Ls, eLs, N):
    Psout = np.random.normal(P, eP, int(N))
    Ps = unp.uarray(Psout, np.repeat(eP, int(N)))
    smas = rvs.semimajoraxis(Ps, unp.uarray(Ms,eMs), 0)
    Fsout = unp.nominal_values(compute_F(unp.uarray(Ls,eLs), smas))
    rpsout = np.random.normal(rp, erp, int(N))
    return Psout, Fsout, rpsout


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



if __name__ == '__main__':
    folder = sys.argv[1]
    self = OccurrenceRateclass(folder, compute_detections=True, compute_sens=True, compute_occurrence_rate=True)
