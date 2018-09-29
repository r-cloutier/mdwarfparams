from K2LCclass import *

class K2sensitivity:

    def __init__(self, epicnum):
        self.epicnum = epicnum
        self.fname_full = 'PipelineResults/EPIC_%isens'%self.epicnum
        self.get_data()
        self.compute_sensitivity()
        self._pickleobject()


    def get_data(self):
        self.fs = np.array(glob.glob('PipelineResults/EPIC_%i_*/*'%self.epicnum))
        self.Nsim = self.fs.size
        assert self.Nsim > 0
        d = loadpickle(self.fs[0])
        self.Kepmag, self.logg, self.Ms, self.Rs, self.Teff = d.Kepmag, \
                                                              d.logg, \
                                                              d.Ms, d.Rs, \
                                                              d.Teff
        Nmax = 20
        Nplanets = np.zeros(self.Nsim)
        self.Ps, self.rps = np.zeros((0,Nmax)), np.zeros((0,Nmax))
        self.isdet = np.zeros((0,Nmax))
        for i in range(self.Nsim):

            print float(i) / self.Nsim
            d = loadpickle(self.fs[i])

            Nplanets[i] = d.Ptrue.size
            filler = np.repeat(np.nan, Nmax-Nplanets[i])
            Pin = np.append(d.Ptrue, filler)
            rpin = np.append(d.rptrue, filler)
            isdetin = np.append(d.is_detected, filler)
            
            self.Ps  = np.append(self.Ps, Pin.reshape(1,Nmax), axis=0)
            self.rps = np.append(self.rps, rpin.reshape(1,Nmax), axis=0)
            self.isdet = np.append(self.isdet, isdetin.reshape(1,Nmax), axis=0)

        # trim excess planets
        
        

    def compute_sensitivity(self):
        '''Get all the simulations for this star and compute the
        sensitivity.'''
        xlen, ylen = 11, 9
        self.Pgrid = np.logspace(-1, np.log10(30), xlen+1)
        self.rpgrid = np.linspace(.5, 4, ylen+1)
        self.Ndet, self.Ntrue = np.zeros((xlen, ylen)), \
                                np.zeros((xlen, ylen))
        for i in range(xlen):
            for j in range(ylen):
                g = (self.Ps >= self.Pgrid[i]) & \
                    (self.Ps <= self.Pgrid[i+1]) & \
                    (self.rps >= self.rpgrid[j]) & \
                    (self.rpgrid <= self.rpgrid[j+1])
                self.Ndet[i,j] = self.isdet[g].sum()
                self.Ntrue[i,j] = self.isdet[g].size

        # compute sensitivity
        self.sens = self.Ndet / self.Ntrue.astype(float)
        self.esens = np.sqrt(self.Ndet) / self.Ntrue.astype(float)

                      
    def _pickleobject(self):
        fObj = open(self.fname_full, 'wb')
        pickle.dump(self, fObj)
        fObj.close()


if __name__ == '__main__':
    self = K2sensitivity(215096532)
