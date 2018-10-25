from K2LCclass import *
from uncertainties import unumpy as unp
import rvs

class K2results:

    def __init__(self, folder, xlen=120, ylen=60):
	self.folder = folder
	self.fname_out = '%s/EPIC_K2results'%self.folder
	self._xlen, self._ylen = int(xlen), int(ylen)	

	self.get_results()
        self.compute_detections()
	self.save_stars_with_detections()
	self._pickleobject()


    def get_results(self):
        # planet search files
        fs = np.array(glob.glob('%s/EPIC_*/K2LC_-00099'%self.folder))
	self.fs, self.epicnums = [], np.zeros(0)

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

	            self.fs.append(fs[i])
		    self.epicnums = np.append(self.epicnums, d.epicnum)
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
                    
                    
  	_, self.unique_inds = np.unique(self.epicnums, return_index=True)
	self.Nstar = self.unique_inds.size
        assert self.unique_inds.size == self.Nstar
	self.Nplanets = self.Ndetected[self.unique_inds].sum()
	self.fs = np.array(self.fs)
	#assert self.Nplanets == self.params_guess.shape[0]

        

    def compute_detections(self, Ntrials=1e3, Plims=(.5,80), rplims=(.5,10)):
        '''compute detections over P and rp using MC simulations over the 
        radius uncertainties.'''
        Nsims, Ntrials = self.fs.size, int(Ntrials)
        self.Ps_MC = np.zeros((Nsims, Ntrials))
        self.rps_MC = np.zeros((Nsims, Ntrials))
        
        # for each detected planet around each star
        for i in range(Nsims):

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

        # compute detection map over P and rp
        self.Pgrid = np.logspace(np.log10(Plims[0]), np.log10(Plims[1]),
                                 self._xlen+1)
        self.rpgrid = np.logspace(np.log10(rplims[0]), np.log10(rplims[1]),
                                  self._ylen+1)
        self.Ndet = np.zeros((self._xlen, self._ylen))
        for i in range(self._xlen):
            for j in range(self._ylen):
                g = (self.Ps_MC >= self.Pgrid[i]) & \
                    (self.Ps_MC <= self.Pgrid[i+1]) & \
                    (self.rps_MC >= self.rpgrid[j]) & \
                    (self.rps_MC <= self.rpgrid[j+1])
                self.Ndet[i,j] = g.sum() / float(Ntrials)


    def save_stars_with_detections(self):
        epic_tmp = self.epicnums[self.unique_inds]
        epicnum2save = epic_tmp[self.Ndetected[self.unique_inds] > 0]
        np.savetxt('input_data/K2targets/K2Mdwarfs_withdetections.csv',
                   epicnum2save, delimiter=',', fmt='%i')

                
    def _pickleobject(self):
        fObj = open(self.fname_out, 'wb')
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
    
        

if __name__ == '__main__':
    folder = sys.argv[1]  # PipelineResults
    self = K2results(folder)
