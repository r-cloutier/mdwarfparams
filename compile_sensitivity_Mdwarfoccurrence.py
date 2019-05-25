from LCclass import *


class MdwarfOccurrence:

    def __init__(self, folder):
        self.folder = folder

        # get the simulated systems and save all parameters
        self.get_injected_systems()
        self.retrieve_data()
        self._pickleobject()
        

    def get_injected_systems(self):
        '''Get all of the simulated systems with injected planets that are
        complete.'''
        fs = np.array(glob.glob('%s/*/*'%self.folder))
        self.fs, self.Nfs = [], 0
        for f in fs:
            if loadpickle(f).DONE:
                self.fs.append(f)
                self.Nfs += 1


    def retrieve_data(self):
        '''Save all the stellar, planetary, observational, and detection
        parameters.'''
        # stellar parameters
        N = self.Nfs
        self.ras, self.decs = np.zeros(N), np.zeros(N)
        self.GBPmags, self.GRPmags = np.zeros((N,2)), np.zeros((N,2))
        self.Kepmags, self.Jmags = np.zeros(N), np.zeros((N,2))
        self.Hmags, self.Kmags = np.zeros((N,2)), np.zeros((N,2))
        self.pars, self.dists = np.zeros((N,2)), np.zeros((N,3))
        self.mus, self.AKs = np.zeros((N,3)), np.zeros((N,3))
        self.MKs, self.Rss = np.zeros((N,3)), np.zeros((N,3))
        self.Teffs, self.Mss = np.zeros((N,3)), np.zeros((N,3))
        self.loggs = np.zeros((N,3))

        # planet parameters
        Nmax = 20
        self.Ps_inj, self.rps_inj = np.zeros((N,Nmax)), np.zeros((N,Nmax))
        self.multiplicities_inj = np.zeros(N)

        # light curve parameters
        self.CDPPs_inj, self.CDPPs_rec = np.zeros((N,Nmax)), np.zeros((N,Nmax))
        self.Zs_inj, self.Zs_rec = np.zeros((N,Nmax)), np.zeros((N,Nmax))
        self.Ntransits_inj, self.Ntransits_rec = np.zeros((N,Nmax)), \
                                                 np.zeros((N,Nmax))
        self.SNRs_inj, self.SNRs_rec = np.zeros((N,Nmax)), np.zeros((N,Nmax))
        self.is_det = np.zeros((N,Nmax), dtype=bool)
        
        for i in range(N):
            # save stellar parameters
            d = loadpickle(self.fs[i])
            assert d.DONE
            self.ras[i], self.decs[i] = d.ra, d.dec
            self.GBPmags[i] = d.GBPmag, d.e_GBPmag
            self.GRPmags[i] = d.GRPmag, d.e_GRPmag
            self.Kepmags[i] = d.Kepmag
            self.Jmags[i] = d.Jmag, d.e_Jmag
            self.Hmags[i] = d.Hmag, d.e_Hmag
            self.Kmags[i] = d.Kmag, d.e_Kmag
            self.pars[i] = d.par, d.e_par
            self.dists[i] = d.dist, d.ehi_dist, d.elo_dist
            self.mus[i] = d.mu, d.ehi_mu, d.elo_mu
            self.AKs[i] = d.AK, d.ehi_AK, d.elo_AK
            self.MKs[i] = d.MK, d.ehi_MK, d.elo_MK
            self.Rss[i] = d.Rs, d.ehi_Rs, d.elo_Rs
            self.Teffs[i] = d.Teff, d.ehi_Teff, d.elo_Teff
            self.Mss[i] = d.Ms, d.ehi_Ms, d.elo_Ms
            self.loggs[i] = d.logg, d.ehi_logg, d.elo_logg

            # planet parameters
            filler = np.repeat(np.nan, Nmax-d.Ptrue.size)
            self.Ps_inj[i] = np.append(d.Ptrue, filler)
            self.rps_inj[i] = np.append(d.rptrue, filler)
            self.multiplicities[i] = d.Ptrue.size

            # light cure parameters
            self.CDPPs_inj[i], self.CDPPs_rec[i] = d.CDPPS_inj, d.CDPPS_rec
            self.Zs_inj[i], self.Zs_rec[i] = d.depths_inj, d.depths_rec
            self.Ntransits_inj[i] = d.Ntransits_inj
            self.Ntransits_rec[i] = d.Ntransits_rec
            self.SNRs_inj[i] = d.SNRtransits_inj
            self.SNRs_rec[i] = d.SNRtransits_rec
            self.is_det[i] = d.is_detected

            
    def _pickleobject(self):
        fObj = open('%s/sens_Mdwarf'%self.folder, 'wb')
        pickle.dump(self, fObj)
        fObj.close()


if __name__ == '__main__':
    folder = sys.argv[1]
    self = MdwarfOccurrence(folder)
