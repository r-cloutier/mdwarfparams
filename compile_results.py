from K2LCclass import *


class K2results:

    def __init__(self, folder, prefix):
	self.folder, self.prefix = folder, prefix
	self.fname_out = '%s/%sK2results'%(self.folder, self.prefix)
	
	self.get_results()
	self._pickleobject()


    def get_results(self):
	fs = np.array(glob.glob('%s/%s*/*'%(self.folder, self.prefix)))
	self.fs = []
	self.Ndetected, self.params_guess = np.zeros(0), np.zeros((0,4))
        self.Kepmags, self.Mss, self.Rss, self.Teffs = np.zeros(0), np.zeros(0), \
                                                       np.zeros(0), np.zeros(0)
	for i in range(fs.size):
	    print float(i)/fs.size
	    d = loadpickle(fs[i])
	    if d.DONE:
		Ndet = int(d.params_guess.shape[0])
		for j in range(Ndet):
		    self.fs.append(fs[i])
                    self.Kepmags = np.append(self.Kepmags, d.Kepmag)
                    self.Mss = np.append(self.Mss, d.Ms)
                    self.Rss = np.append(self.Rss, d.Rs)
                    self.Teffs = np.append(self.Teffs, d.Teff)
		    self.Ndetected = np.append(self.Ndetected, Ndet)
		    self.params_guess = np.append(self.params_guess, d.params_guess, axis=0)
  	
  	_, self.unique_inds = np.unique(self.fs, return_index=True)
	self.Nstar = self.unique_inds.size
	self.Nplanets = self.Ndetected[self.unique_inds].sum()
	assert self.Nplanets == self.params_guess.shape[0]


    def _pickleobject(self):
        fObj = open(self.fname_out, 'wb')
        pickle.dump(self, fObj)
        fObj.close()


if __name__ == '__main__':
    folder, prefix = sys.argv[1], sys.argv[2]  # PipelineResults EPIC_
    self = K2results(folder, prefix)
