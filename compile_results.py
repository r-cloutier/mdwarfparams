from K2LCclass import *


class K2results:

    def __init__(self, folder, prefix):
	self.folder, self.prefix = folder, prefix
	self.fname_out = '%s/%sK2results'%(self.folder, self.prefix)
	
	self.get_results()
	self._pickleobject()


    def get_results(self):
	fs = np.array(glob.glob('%s/%s*/*'%(self.folder, self.prefix)))
	self.fs, self.epicnames = [], np.zeros(0)
	self.Ndetected, self.params_guess, self.cond_vals = np.zeros(0), np.zeros((0,4)), \
                                                            np.zeros((0,3))
        self.Kepmags, self.Mss, self.Rss, self.Teffs = np.zeros(0), np.zeros(0), \
                                                       np.zeros(0), np.zeros(0)
	self.efs = np.zeros(0)
	for i in range(fs.size):
	    print float(i)/fs.size
	    d = loadpickle(fs[i])
	    if d.DONE:
		Ndet = int(d.params_guess.shape[0])
		for j in range(Ndet+1):
		    self.fs.append(fs[i])
		    self.epicnames = np.append(self.epicnames, int(fs[i].split('_')[1].split('/')[0]))
                    self.Kepmags = np.append(self.Kepmags, d.Kepmag)
                    self.Mss = np.append(self.Mss, d.Ms)
                    self.Rss = np.append(self.Rss, d.Rs)
                    self.Teffs = np.append(self.Teffs, d.Teff)
		    self.efs = np.append(self.efs, d.ef.mean())
		    self.Ndetected = np.append(self.Ndetected, Ndet)
                    
                    params = d.params_guess[j-1] if j > 0 else np.repeat(np.nan,4)
		    self.params_guess = np.append(self.params_guess, params.reshape(1,4), axis=0)
                    P = params[0]

                    g = np.isclose(d.params_guess_priorto_confirm[:,0], P, rtol=.05)
		    print self.epicnames[-1], g.sum()
                    cond_vals = [d.transit_condition_scatterin_val[g], \
                                 d.transit_condition_depth_val[g], \
                                 d.transit_condition_no_bimodal_val[g]] if j > 0 else [np.nan]*3
                    self.cond_vals = np.append(self.cond_vals, np.array(cond_vals).reshape(1,3),
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
    folder, prefix = sys.argv[1], sys.argv[2]  # PipelineResults EPIC_
    self = K2results(folder, prefix)
