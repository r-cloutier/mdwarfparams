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
        self.isdet, self.isFP = np.zeros((0,Nmax)), np.zeros((0,Nmax))
        for i in range(self.Nsim):

            print float(i) / self.Nsim
            d = loadpickle(self.fs[i])

            Nplanets[i] = d.Ptrue.size
            filler = np.repeat(np.nan, Nmax-Nplanets[i])
            Pin = np.append(d.Ptrue, filler)
            rpin = np.append(d.rptrue, filler)
            isdetin = np.append(d.is_detected, filler)
            
            self.Ps = np.append(self.Ps, Pin.reshape(1,Nmax), axis=0)
            self.rps = np.append(self.rps, rpin.reshape(1,Nmax), axis=0)
            self.isdet = np.append(self.isdet, isdetin.reshape(1,Nmax), axis=0)

            # get false positives TEMP because d.is_FP doesnt exist yet until
            # the simulation is rerun
            params = d.params_guess
            is_FP =  np.array([int(np.invert(np.any(np.isclose(d.Ptrue,
                                                               params[i,0],
                                                               rtol=.02))))
                               for i in range(params.shape[0])]).astype(bool)
            filler2 = np.repeat(np.nan, Nmax-is_FP.size)
            self.isFP = np.append(self.isFP,
                                  np.append(is_FP,filler2).reshape(1,Nmax),
                                  axis=0)
            
        # trim excess planets
        end = np.where(np.all(np.isnan(self.Ps), axis=0))[0][0]
        self.Ps = self.Ps[:,:end]
        self.rps = self.rps[:,:end]
        self.isdet = self.isdet[:,:end]
        end2 = np.where(np.all(np.isnan(self.isFP), axis=0))[0][0]
        self.isFP = self.isFP[:,:end2]
        

    def compute_sensitivity(self):
        '''Get all the simulations for this star and compute the
        sensitivity and the number of FPs as functions of P and rp.'''
        xlen, ylen = 11, 9
        self.Pgrid = np.logspace(-1, np.log10(30), xlen+1)
        self.rpgrid = np.linspace(.5, 4, ylen+1)
        self.Ndet, self.Ntrue, self.NFP = np.zeros((xlen, ylen)), \
                                          np.zeros((xlen, ylen)), \
                                          np.zeros((xlen, ylen))
        for i in range(xlen):
            for j in range(ylen):
                g = (self.Ps >= self.Pgrid[i]) & \
                    (self.Ps <= self.Pgrid[i+1]) & \
                    (self.rps >= self.rpgrid[j]) & \
                    (self.rps <= self.rpgrid[j+1])
                self.Ndet[i,j] = self.isdet[g].sum()
                self.Ntrue[i,j] = self.isdet[g].size
                self.NFP[i,j] = self.isFP[g].sum()
                
        # compute sensitivity
        self.sens = self.Ndet / self.Ntrue.astype(float)
        self.esens = np.sqrt(self.Ndet) / self.Ntrue.astype(float)

        # compute yield correction to multiply the yield by
        self.yield_corr = 1 - self.NFP / (self.Ndet + self.NFP.astype(float))
        
                              
    def _pickleobject(self):
        fObj = open(self.fname_full, 'wb')
        pickle.dump(self, fObj)
        fObj.close()


if __name__ == '__main__':
    self = K2sensitivity(215096532)
