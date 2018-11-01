from K2LCclass import *
from K2sensclass import *
from uncertainties import unumpy as unp
import rvs
from scipy.ndimage.filters import gaussian_filter # for map smoothing if desired


class K2occurrencerate:

    def __init__(self, folder, xlen=20, ylen=12, compute_detections=False,
                 compute_sens=False, comupute_occurrence_rate=False,
                 Plims=(.5,80), rplims=(.5,10)):
        self.folder = folder
        self._xlen, self._ylen = int(xlen), int(ylen)
        self.Plims, self.rplims = Plims, rplims
        
        if compute_detections:
            self.get_planetsearch_results()
            self.fname_planetsearch_results = '%s/EPIC_K2results'%self.folder
            self._pickleobject(self.fname_planetsearch_results)
            
        if compute_sens:
            self.get_simulation_results()
            self.fname_sens_results = '%s/EPIC_K2sens'%self.folder
            self._pickleobject(self.fname_sens_results)

        if comupute_occurrence_rate:
            self.compute_occurrence_rate()


    def get_planetsearch_results(self):
        '''Get the results from the planetsearch, i.e. the detected planets 
        and stellar properties.'''
        fs = np.array(glob.glob('%s/EPIC_*/K2LC_-00099'%self.folder))
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
        self.rps_MC = np.zeros((N, Ntrials))

        # for each detected planet around each star, compute realizations
        # over the uncertainties
        for i in range(N):

            if self.Ndetected[i] > 0:
                
                # compute MC realizations of this planet
                P, eP = self.Ps[i], self.e_Ps[i]
                rp, erp = self.rps[i], self.e_rps[i]
                self.Ps_MC[i], self.rps_MC[i] = sample_planets(P, eP,
                                                               rp, erp,
                                                               Ntrials)
            else:
                self.Ps_MC[i], self.rps_MC[i] = np.repeat(np.nan, Ntrials), \
                                                np.repeat(np.nan, Ntrials)

        # compute Ndet map over P and rp for each star
        self.logPgrid = np.logspace(np.log10(self.Plims[0]),
                                    np.log10(self.Plims[1]), self._xlen+1)
        self.logrpgrid = np.logspace(np.log10(self.rplims[0]),
                                     np.log10(self.rplims[1]), self._ylen+1)
        self.Ndet_i = np.zeros((self.Nstars, self._xlen, self._ylen))
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
                    self.Ndet_i[i,j,k] = g2.sum() / float(Ntrials)



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
        Nmax = 20
        self.Nplanets_inj = np.zeros(self.Nstars_wdet, 0)
        self.Nplanets_rec = np.zeros(self.Nstars_wdet, 0)
        self.Ps_inj = np.zeros((self.Nstars_wdet, 0, Nmax))
        self.rps_inj = np.zeros((self.Nstars_wdet, 0, Nmax))
        self.is_rec = np.zeros((self.Nstars_wdet, 0, Nmax))
        self.Ps_rec = np.zeros((self.Nstars_wdet, 0, Nmax))
        self.rps_rec = np.zeros((self.Nstars_wdet, 0, Nmax))
        self.is_FP = np.zeros((self.Nstars_wdet, 0, Nmax))
        
        for i in range(self.Nstars_wdet):

            epicnum = self.epicnums_wdet[i]
            fs = np.array(glob.glob('PipelineResults/EPIC_%i/K2LC*'%epicnum))

	    # remove planet search result (i.e. with index -99)
            g = np.in1d(self.fs, '%s/EPIC_%i/K2LC_-00099'%(self.folder,epicnum))
	    if np.any(g):
	        fs = np.delete(fs, np.where(g)[0][0])
            if self.fs.size == 0:
	        return None

            # get params of injected and recovered planets for this star
            for j in range(fs.size):

                print float(j) / fs.size
                d = loadpickle(fs[j])
                if d.DONE:
                    self.Nsims[i] += 1
                    self.Nplanets_inj[i] = np.append(self.Nplanets_inj[i],
                                                     d.Ptrue.size)
                    filler = np.repeat(np.nan, Nmax-self.Nplanets_inj[i,-1])
	            Pin = np.append(d.Ptrue, filler)
            	    rpin = np.append(d.rptrue, filler)
            	    isrecin = np.append(d.is_detected, filler)
            	    self.Ps_inj[i] = np.append(self.Ps_inj[i],
                                               Pin.reshape(1,Nmax), axis=0)
            	    self.rps_inj[i] = np.append(self.rps_inj[i],
                                                rpin.reshape(1,Nmax), axis=0)
            	    self.is_rec[i] = np.append(self.is_rec[i],
                                               isrecin.reshape(1,Nmax), axis=0)

                    # get false positives 
            	    params = d.params_guess
		    self.Nplanets_rec[i] = np.append(self.Nplanets_rec[i],
                                                     params.shape[0])
		    filler2 = np.repeat(np.nan, Nmax-self.Nplanets_rec[i,-1])
                    Pinrec = np.append(params[:,0], filler2)
                    rp = rvs.m2Rearth(rvs.Rsun2m(np.sqrt(params[:,2])*d.Rs))
                    rpinrec = np.append(rp, filler2)
                    isFPin = np.append(d.is_FP, filler2)
                    self.Ps_rec[i] = np.append(self.Ps_rec[i],
                                               Pinrec.reshape(1,Nmax), axis=0)
                    self.rps_rec[i] = np.append(self.rps_rec[i],
                                                rpinrec.reshape(1,Nmax), axis=0)
            	    self.is_FP[i] = np.append(self.is_FP[i],
                                              isFPin.reshape(1,Nmax), axis=0)

        # trim excess planets
        end = np.where(np.all(np.isnan(self.Ps_inj[:,:,])))[0][0] #TEMP
        self.Ps_inj = self.Ps_inj[:,:,:end]
        self.rps_inj = self.rps_inj[:,:,:end]
        self.is_rec = self.is_rec[:,:,:end]
        self.Ps_rec = self.Psfound[:,:,:end]
	self.rps_rec = self.rpsfound[:,:,:end]
        self.is_FP = self.isFP[:,:,:end]

        # compute sensitivity and transit probability maps
        self.compute_sens_maps()
        self.compute_transitprob_maps()

        
        
    def compute_sens_maps(self):
        '''Get all the simulations for all stars and compute the
        sensitivity and the number of FPs as functions of P and rp.'''
        self.logPgrid = np.logspace(np.log10(self.Plims[0]),
                                    np.log10(self.Plims[1]), self._xlen+1)
        self.logrpgrid = np.logspace(np.log10(self.rplims[0]),
                                     np.log10(self.rplims[1]), self._ylen+1)
        self.Nrec_i = np.zeros((self.Nstars_wdet, self._xlen, self._ylen))
        self.Ninj_i = np.zeros((self.Nstars_wdet, self._xlen, self._ylen))
        self.NFP_i  = np.zeros((self.Nstars_wdet, self._xlen, self._ylen))

        for i in range(self.Nstars_wdet):
            for j in range(self._xlen):
                for k in range(self._ylen):

                    g = (self.Ps_inj[i] >= self.logPgrid[j]) & \
                        (self.Ps_inj[i] <= self.logPgrid[j+1]) & \
                        (self.rps_inj[i] >= self.logrpgrid[k]) & \
                        (self.rps_inj[i] <= self.logrpgrid[k+1])
                    self.Nrec_i[i,j,k] = self.is_rec[g].sum()
                    self.Ninj_i[i,j,k] = self.is_rec[g].size

		    g = (self.Ps_rec[i] >= self.logPgrid[j]) & \
                        (self.Ps_rec[i] <= self.logPgrid[j+1]) & \
                        (self.rps_rec[i] >= self.logrpgrid[k]) & \
                        (self.rps_rec[i] <= self.logrpgrid[k+1])
                    self.NFP_i[i,j,k] = self.is_FP[g].sum()

        # compute sensitivity
        self.sens_i = self.Nrec_i / self.Ninj_i.astype(float)
        self.e_sens_i = np.sqrt(self.Nrec_i) / self.Ninj_i.astype(float)

        # compute yield correction to multiply Ndet by
        self.yield_corr_i = 1 - self.NFP_i / \
                            (self.Nrec_i + self.NFP_i.astype(float))
        self.e_yield_corr_i = np.sqrt(self.NFP_i) / \
                            (self.Nrec_i + self.NFP_i.astype(float))

        

    def compute_transitprob_maps(self):
        '''Compute the transiting probability maps for each star with a
        detected planet candidate.'''
        self.transit_prob_i = np.zeros_like(self.sens_i)
        self.e_transit_prob_i = np.zeros_like(self.sens_i)
        for i in range(self.Nstars_wdet):
            for j in range(self._xlen):
                for k in range(self._ylen):
                    
                    Pmid = 10**(np.log10(self.logPgrid[j]) + \
                                np.diff(np.log10(self.logPgrid[:2])/2))
                    rpmid = 10**(np.log10(self.logrpgrid[k]) + \
                                 np.diff(np.log10(self.logrpgrid[:2])/2))
                    
                    epicnum = self.epicnums_wdet[i]
                    g = np.where(self.epicnums_planetsearch == epicnum)[0][0]
                    Ms = unp.uarray(self.Mss[g], self.e_Ms[g])
                    sma = rvs.AU2m(rvs.semimajoraxis(Pmid, Ms, 0))
                    
                    Rs = unp.uarray(self.Rss[g], self.e_Rss[g])
                    prob = (rvs.Rsun2m(Rs) + rvs.Rearth2m(rpmid)) / sma
                    self.transit_prob_i[i,j,k] = unp.nominal_values(prob)
                    self.e_transit_prob_i[i,j,k] = unp.std_devs(prob)

	# correction from beta distribution fit (Kipping 2013)
	self.transit_prob_factor = 1.08
	self.transit_prob_i *= self.transit_prob_factor
	self.e_transit_prob_i *= self.transit_prob_factor


        
    def compute_occurrence_rate():
        '''Use the maps of Ndet, yield_corr, sensitivity, and transit 
        probability to compute the occurrence rate over P and rp.'''
        assert self.Ndet_i.shape[0] == self.Nstars
        assert self.sens_i.shape[0] == self.Nstars_wdet
        assert self.Ndet_i.shape[1:] == self.sens_i.shape[1:]
        
        self.occurrence_rate_i = np.zeros_like(self.sens_i)
        for i in range(self.Nstars_wdet):

            g = self.epicnums_planetsearch[self.unique_inds] == \
                self.epicnums_wdet[i] 
            Ndet_i = self.Ndet_i[g]
            S_i = self.sens_i[i]
            C_i = self.yield_corr_i[i]
            t_i = self.transit_prob_i[i]
            
            self.occurrence_rate_i[i] = Ndet_i * C_i / (S_i * t_i)

        # compute occurrence rate over P,rp
        self.occurrence_rate = np.nansum(self.occurrence_rate_i, axis=0)
            


    def _pickleobject(self, fname):
        fObj = open(fname, 'wb')
        pickle.dump(self, fObj)
        fObj.close()

                    

def get_1sigma(results):
    '''results = med, plus_1sig, minus_1sig'''
    assert results.size == 3
    return np.mean(results[1:])


def rpRs2rp(rpRs, Rs):
    rp = rvs.m2Rearth(rvs.Rsun2m(rpRs*Rs))
    return unp.nominal_values(rp), unp.std_devs(rp)


def sample_planets(P, eP, rp, erp, N):
    Psout = np.random.normal(P, eP, int(N))
    rpsout = np.random.normal(rp, erp, int(N))
    return Psout, rpsout


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
