from K2LCclass import *


class K2results:

    def __init__(self, folder):
	self.folder = folder
	self.fname_out = '%s/EPIC_K2results'%self.folder
	
	self.get_results()
	self._pickleobject()


    def get_results(self):
	fs = np.array(glob.glob('%s/EPIC_*/K2LC_-00099'%self.folder))
	self.fs, self.epicnames = [], np.zeros(0)
	self.Ndetected, self.params_guess, self.cond_vals = np.zeros(0), np.zeros((0,4)), \
                                                            np.zeros((0,5))
        self.Kepmags, self.efs = np.zeros(0), np.zeros(0)
	self.Mss, self.e_Mss = np.zeros(0), np.zeros(0)
        self.Rss, self.e_Rss = np.zeros(0), np.zeros(0)
        self.Teffs, self.e_Teffs = np.zeros(0), np.zeros(0)
        self.loggs, self.e_loggs = np.zeros(0), np.zeros(0)
	for i in range(fs.size):
	    print float(i)/fs.size
	    d = loadpickle(fs[i])
	    if d.DONE:
		Ndet = int(d.params_guess.shape[0])
		for j in range(Ndet+1):
		    self.fs.append(fs[i])
		    self.epicnames = np.append(self.epicnames, int(fs[i].split('_')[1].split('/')[0]))
                    self.Kepmags = np.append(self.Kepmags, d.Kepmag)
                    self.efs = np.append(self.efs, d.ef.mean())
                    self.Mss = np.append(self.Mss, d.Ms)
		    self.e_Mss = np.append(self.e_Mss, d.e_Ms)
                    self.Rss = np.append(self.Rss, d.Rs)
                    self.e_Rss = np.append(self.e_Rss, d.e_Rs)
                    self.Teffs = np.append(self.Teffs, d.Teff)
                    self.e_Teffs = np.append(self.e_Teffs, d.e_Teff)
		    self.Ndetected = np.append(self.Ndetected, Ndet)
                    
                    params = d.params_guess[j-1] if j > 0 else np.repeat(np.nan,4)
		    self.params_guess = np.append(self.params_guess, params.reshape(1,4), axis=0)
                    P = params[0]

                    #g = np.isclose(d.params_guess_priorto_confirm[:,0], P, rtol=.05)
                    Pss = d.params_guess_priorto_confirm[:,0]
                    g = abs(Pss-P) == np.min(abs(Pss-P))
		    print self.epicnames[-1], g.sum()
                    cond_vals = [d.transit_condition_scatterin_val[g], \
                                 d.transit_condition_depth_val[g], \
                                 d.transit_condition_no_bimodal_val[g], \
				 d.transit_condition_timesym_val[g], \
				 d.transit_condition_indiv_transit_frac_val[g]] if j > 0 else [np.nan]*5
                    self.cond_vals = np.append(self.cond_vals, np.array(cond_vals).reshape(1,5),
                                               axis=0)
                    
  	_, self.unique_inds = np.unique(self.epicnames, return_index=True)
	self.Nstar = self.unique_inds.size
	self.Nplanets = self.Ndetected[self.unique_inds].sum()
	self.fs = np.array(self.fs)
	#assert self.Nplanets == self.params_guess.shape[0]


    def _pickleobject(self):
        fObj = open(self.fname_out, 'wb')
        pickle.dump(self, fObj)
        fObj.close()


if __name__ == '__main__':
    folder = sys.argv[1]  # PipelineResults
    self = K2results(folder)
