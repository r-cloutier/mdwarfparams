from K2LCclass import *
import rvs

class K2sensitivity:

    def __init__(self, epicnum):
        self.epicnum = epicnum
        self.fname_full = 'PipelineResults/EPIC_%i/EPIC_%i_sens'%(self.epicnum, self.epicnum)
        self.get_data()
        #self.compute_sensitivity()
	#self.compute_transit_prob()
        self._pickleobject()


    def get_data(self):
        self.fs = np.array(glob.glob('PipelineResults/EPIC_%i/K2LC*'%self.epicnum))
        assert self.fs.size > 0
        self.Nsim = 0
        d = loadpickle(self.fs[0])
        self.Kepmag, self.logg, self.Ms, self.Rs, self.Teff = d.Kepmag, \
                                                              d.logg, \
                                                              d.Ms, d.Rs, \
                                                              d.Teff
        Nmax = 20
        self.Nplanets = np.zeros(0)
        self.Ps, self.rps, self.isdet = np.zeros((0,Nmax)), np.zeros((0,Nmax)), np.zeros((0,Nmax))
        self.PsFP, self.rpsFP, self.isFP = np.zeros((0,Nmax)), np.zeros((0,Nmax)), np.zeros((0,Nmax))
        for i in range(self.fs.size):

            print float(i) / self.fs.size
            d = loadpickle(self.fs[i])

	    if d.DONE:
		self.Nsim += 1
            	self.Nplanets = np.append(self.Nplanets, d.Ptrue.size)
            	filler = np.repeat(np.nan, Nmax-self.Nplanets[-1])
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
		PinFP = np.append(params[is_FP,0], filler2)
		rpinFP = np.append(rvs.m2Rearth(rvs.Rsun2m(params[is_FP,2]*d.Rs)), filler2)
		self.PsFP = np.append(self.PsFP, PinFP.reshape(1,Nmax), axis=0)
		self.rpsFP = np.append(self.rpsFP, rpinFP.reshape(1,Nmax), axis=0)
            	self.isFP = np.append(self.isFP,
                                      np.append(is_FP,filler2).reshape(1,Nmax),
                                      axis=0)
            
        # trim excess planets
        end = np.where(np.all(np.isnan(self.Ps), axis=0))[0][0]
        self.Ps = self.Ps[:,:end]
        self.rps = self.rps[:,:end]
        self.isdet = self.isdet[:,:end]
        end2 = np.where(np.all(np.isnan(self.isFP), axis=0))[0][0]
	self.PsFP = self.PsFP[:,:end2]
	self.rpsFP = self.rpsFP[:,:end2]
        self.isFP = self.isFP[:,:end2]
        

    def compute_sensitivity(self, Plims=(.5,30), rplims=(.5,4)):
        '''Get all the simulations for this star and compute the
        sensitivity and the number of FPs as functions of P and rp.'''
        xlen, ylen = 11, 7
        self.Pgrid = np.logspace(np.log10(Plims[0]), np.log10(Plims[1]), xlen+1)
        self.rpgrid = np.linspace(rplims[0], rplims[1], ylen+1)
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
		
		g = (self.PsFP >= self.Pgrid[i]) & \
                    (self.PsFP <= self.Pgrid[i+1]) & \
                    (self.rpsFP >= self.rpgrid[j]) & \
                    (self.rpsFP <= self.rpgrid[j+1])
                self.NFP[i,j] = self.isFP[g].sum()
                
        # compute sensitivity
        self.sens = self.Ndet / self.Ntrue.astype(float)
        self.esens = np.sqrt(self.Ndet) / self.Ntrue.astype(float)

        # compute yield correction to multiply the yield by
        self.yield_corr = 1 - self.NFP / (self.Ndet + self.NFP.astype(float))
        
    
    def compute_transit_prob(self):
	self.transit_prob = np.zeros((self.Pgrid.size-1, self.rpgrid.size-1))
	for i in range(self.Pgrid.size-1):
	    for j in range(self.rpgrid.size-1):
		Pmid = 10**(np.log10(self.Pgrid[i]) + np.diff(np.log10(self.Pgrid[i:i+2])/2))
		rpmid = np.mean(self.rpgrid[j:j+2])
		sma = rvs.AU2m(rvs.semimajoraxis(Pmid, self.Ms, 0))
		self.transit_prob[i,j] = (rvs.Rsun2m(self.Rs) + rvs.Rearth2m(rpmid)) / sma

	# correction from beta distribution fit (Kipping 2013)
	self.transit_prob *= 1.08

                          
    def _pickleobject(self):
        fObj = open(self.fname_full, 'wb')
        pickle.dump(self, fObj)
        fObj.close()


if __name__ == '__main__':
    #fs = np.array(glob.glob('PipelineResults/EPIC_*'))
    #print fs.size
    #for i in range(fs.size):
#	print fs[i]
#	epicnum = int(fs[i].split('_')[-1])
    self = K2sensitivity(201134718)
