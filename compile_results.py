from K2LCclass import *
from depth2rp import *


class K2results:

    def __init__(self, folder):
	self.folder = folder
	self.fname_out = '%s/EPIC_K2results'%self.folder
	
	self.get_results()
        self.compute_detections()
	self._pickleobject()


    def get_results(self):
	fs = np.array(glob.glob('%s/EPIC_*/K2LC_-00099'%self.folder))
	self.fs, self.epicnames = [], np.zeros(0)
	self.Ndetected, self.params_guess, self.cond_vals = np.zeros(0), \
                                                            np.zeros((0,4)), \
                                                            np.zeros((0,5))
        self.Kepmags, self.efs = np.zeros(0), np.zeros(0)
	self.Mss, self.e_Mss = np.zeros(0), np.zeros(0)
        self.Rss, self.e_Rss = np.zeros(0), np.zeros(0)
        self.Teffs, self.e_Teffs = np.zeros(0), np.zeros(0)
        self.loggs, self.e_loggs = np.zeros(0), np.zeros(0)
	for i in range(fs.size):
	    print float(i)/fs.size, fs[i]
	    d = loadpickle(fs[i])
	    if d.DONE:
		Ndet = int(d.params_guess.shape[0])
		for j in range(Ndet+1):
		    self.fs.append(fs[i])
		    self.epicnames = np.append(self.epicnames,
                                               int(fs[i].split('_')[1].split('/')[0]))
                    self.Kepmags = np.append(self.Kepmags, d.Kepmag)
                    self.efs = np.append(self.efs, d.ef.mean())
                    self.Mss = np.append(self.Mss, d.Ms)
		    self.e_Mss = np.append(self.e_Mss, d.e_Ms)
                    self.Rss = np.append(self.Rss, d.Rs)
                    self.e_Rss = np.append(self.e_Rss, d.e_Rs)
                    self.Teffs = np.append(self.Teffs, d.Teff)
                    self.e_Teffs = np.append(self.e_Teffs, d.e_Teff)
		    self.Ndetected = np.append(self.Ndetected, Ndet)
                    
                    params = d.params_guess[j-1] if j > 0 \
                             else np.repeat(np.nan,4)
		    self.params_guess = np.append(self.params_guess,
                                                  params.reshape(1,4), axis=0)
                    P = params[0]

                    Pss = d.params_guess_priorto_confirm[:,0]
		    if Pss.size > 0:    
                      	g = abs(Pss-P) == np.min(abs(Pss-P))
		    	assert g.sum() in range(2)
                    	cond_vals = [d.transit_condition_scatterin_val[g], \
                                     d.transit_condition_depth_val[g], \
                                     d.transit_condition_no_bimodal_val[g], \
				     d.transit_condition_timesym_val[g], \
				     d.transit_condition_indiv_transit_frac_val[g]] if j > 0 else [np.nan]*5
                    else:
			cond_vals = np.repeat(np.nan, 5)
		    self.cond_vals = np.append(self.cond_vals,
                                               np.array(cond_vals).reshape(1,5),
                                               axis=0)
                    
  	_, self.unique_inds = np.unique(self.epicnames, return_index=True)
	self.Nstar = self.unique_inds.size
	self.Nplanets = self.Ndetected[self.unique_inds].sum()
	self.fs = np.array(self.fs)
	#assert self.Nplanets == self.params_guess.shape[0]


    def compute_detections(self, Ntrials=1e3, Plims=(.5,30), rplims=(.5,4)):
        '''compute detections over P and rp using MC simulations over the 
        radius uncertainties.'''
        Nsim, Ntrials = self.fs.size, int(Ntrials)
        self.rps_MC = np.zeros((Nsim, Ntrials))
        self.Ndet_MC = np.zeros((xlen, ylen))
        assert self.unique_inds.size == self.Nstar
        
        # for each detected planet around each star
        for i in range(self.fs.size):
            
            # do N MC realizations of rp if planet is detected
            if self.Ndetected[i] > 0:
                for j in range(Ntrials):

                    sample_rp(self.paraself.Rss[i], self.)

                    

        xlen, ylen = 11, 7
        self.Pgrid = np.logspace(np.log10(Plims[0]), np.log10(Plims[1]), xlen+1)
        self.rpgrid = np.linspace(rplims[0], rplims[1], ylen+1)

        

        
    def _pickleobject(self):
        fObj = open(self.fname_out, 'wb')
        pickle.dump(self, fObj)
        fObj.close()



def sample_rp(depth, e_depth, Rs, e_Rs):
    depth2rp(P_days, depth, duration_days, Ms, Rs)
    
    
        

if __name__ == '__main__':
    folder = sys.argv[1]  # PipelineResults
    self = K2results(folder)
