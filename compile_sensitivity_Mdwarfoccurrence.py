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
            if loadpickle(f).DONE and '-' not in f:
                self.fs.append(f)
                self.Nfs += 1


    def retrieve_data(self):
        '''Save all the stellar, planetary, observational, and detection
        parameters.'''
        # stellar parameters
        N = self.Nfs
        self.epicnums = np.zeros(N)
        self.ras, self.decs = np.zeros(N), np.zeros(N)
        self.GBPmags, self.GRPmags = np.zeros((N,2)), np.zeros((N,2))
        self.Kepmags, self.Jmags = np.zeros(N), np.zeros((N,2))
        self.Hmags, self.Kmags = np.zeros((N,2)), np.zeros((N,2))
        self.pars, self.dists = np.zeros((N,2)), np.zeros((N,3))
        self.mus, self.AKs = np.zeros((N,3)), np.zeros((N,2))
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
        self.cond_free_params = np.zeros((N,11))
        self.cond_vals, self.cond_bools = np.zeros((N,Nmax,9)), np.zeros((N,Nmax,9))
        self.is_det = np.zeros((N,Nmax), dtype=bool)
        
        for i in range(N):
            # save stellar parameters
            d = loadpickle(self.fs[i])
            assert d.DONE
            self.epicnums[i] = d.epicnum
            self.ras[i], self.decs[i] = d.ra, d.dec
            self.GBPmags[i] = d.GBPmag, d.e_GBPmag
            self.GRPmags[i] = d.GRPmag, d.e_GRPmag
            self.Kepmags[i] = d.Kepmag
            self.Jmags[i] = d.Jmag, d.e_Jmag
            self.Hmags[i] = d.Hmag, d.e_Hmag
            self.Kmags[i] = d.Kmag, d.e_Kmag
            self.pars[i] = d.par, d.e_par
            self.dists[i] = d.dist, d.ehi_dist, d.elo_dist
            #self.mus[i] = d.mu, d.ehi_mu, d.elo_mu
            self.AKs[i] = d.AK, d.e_AK
            self.MKs[i] = d.MK, d.ehi_MK, d.elo_MK
            self.Rss[i] = d.Rs, d.ehi_Rs, d.elo_Rs
            self.Teffs[i] = d.Teff, d.ehi_Teff, d.elo_Teff
            self.Mss[i] = d.Ms, d.ehi_Ms, d.elo_Ms
            self.loggs[i] = d.logg, d.ehi_logg, d.elo_logg

            # planet parameters
            filler = np.repeat(np.nan, Nmax-d.Ptrue.size)
            self.Ps_inj[i] = np.append(d.Ptrue, filler)
            self.rps_inj[i] = np.append(d.rptrue, filler)
            self.multiplicities_inj[i] = d.Ptrue.size

            # light cure parameters
            self.CDPPs_inj[i] = np.append(d.CDPPs_inj, filler)
            #self.CDPPs_rec[i] = np.append(d.CDPPs_rec, filler)
            self.Zs_inj[i] = np.append(d.depths_inj, filler)
            #self.Zs_rec[i] = np.append(d.depths_rec, filler)
            self.Ntransits_inj[i] = np.append(d.Ntransits_inj, filler)
            #self.Ntransits_rec[i] = np.append(d.Ntransits_rec, filler)
            self.SNRs_inj[i] = np.append(d.SNRtransits_inj, filler)
            #self.SNRs_rec[i] = np.append(d.SNRtransits_rec, filler)
            self.is_det[i] = np.append(d.is_detected, np.repeat(False, Nmax-d.Ptrue.size))

            # transit parameters for human analysis
            Noi = d.transit_condition_values.shape[0]
            self.cond_vals[i] = np.append(d.transit_condition_values, np.repeat(np.nan,9*(Nmax-Noi)).reshape(Nmax-Noi,9), 0)
            self.cond_bools[i] = np.append(d.transit_condition_bool, np.repeat(np.nan,9*(Nmax-Noi)).reshape(Nmax-Noi,9), 0)
            self.cond_free_params[i] = d.transit_condition_free_params
            
            
    def _pickleobject(self):
        fObj = open('%s/sens_Mdwarf'%self.folder, 'wb')
        pickle.dump(self, fObj)
        fObj.close()


if __name__ == '__main__':
    folder = sys.argv[1]
    self = MdwarfOccurrence(folder)
