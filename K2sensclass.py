from K2LCclass import *
import rvs

class K2sensitivity:

    def __init__(self, epicnum):
        self.epicnum = epicnum
        self.fname_full = 'PipelineResults/EPIC_%i/EPIC_%i_sens'%(self.epicnum, self.epicnum)
        self.get_data()
        #self.get_probable_detections()
        self.compute_sensitivity()
	self.compute_transit_prob()
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
        self.Nplanets_true, self.Nplanets_found = np.zeros(0), np.zeros(0)
        self.Ps, self.rps, self.isdet = np.zeros((0,Nmax)), np.zeros((0,Nmax)), np.zeros((0,Nmax))
        self.Psfound, self.rpsfound, self.isFP = np.zeros((0,Nmax)), np.zeros((0,Nmax)), np.zeros((0,Nmax))
        for i in range(self.fs.size):

            print float(i) / self.fs.size
            d = loadpickle(self.fs[i])

	    if d.DONE:
		self.Nsim += 1
            	self.Nplanets_true = np.append(self.Nplanets_true, d.Ptrue.size)
            	filler = np.repeat(np.nan, Nmax-self.Nplanets_true[-1])
	        Pin = np.append(d.Ptrue, filler)
            	rpin = np.append(d.rptrue, filler)
            	isdetin = np.append(d.is_detected, filler)
            
            	self.Ps = np.append(self.Ps, Pin.reshape(1,Nmax), axis=0)
            	self.rps = np.append(self.rps, rpin.reshape(1,Nmax), axis=0)
            	self.isdet = np.append(self.isdet, isdetin.reshape(1,Nmax), axis=0)

            	# get false positives TEMP because d.is_FP doesnt exist yet until
            	# the simulation is rerun
            	params = d.params_guess
		self.Nplanets_found = np.append(self.Nplanets_found, params.shape[0])
		filler2 = np.repeat(np.nan, Nmax-self.Nplanets_found[-1])
                Pinfound = np.append(params[:,0], filler2)
                rpinfound = np.append(rvs.m2Rearth(rvs.Rsun2m(np.sqrt(params[:,2])*d.Rs)), filler2)
                self.Psfound = np.append(self.Psfound, Pinfound.reshape(1,Nmax), axis=0)
                self.rpsfound = np.append(self.rpsfound, rpinfound.reshape(1,Nmax), axis=0)

            	is_FPin =  np.array([int(np.invert(np.any(np.isclose(d.Ptrue,
                                                               	     params[i,0],
                                                                     rtol=.02))))
                                     for i in range(params.shape[0])]).astype(bool)
            	self.isFP = np.append(self.isFP,
                                      np.append(is_FPin,filler2).reshape(1,Nmax),
                                      axis=0)
            
        # trim excess planets
        end = np.where(np.all(np.isnan(self.Ps), axis=0))[0][0]
        self.Ps = self.Ps[:,:end]
        self.rps = self.rps[:,:end]
        self.isdet = self.isdet[:,:end]
        #end2 = np.where(np.all(np.isnan(self.PsFP), axis=0))[0][0]
	self.Psfound = self.Psfound[:,:end]
	self.rpsfound = self.rpsfound[:,:end]
        self.isFP = self.isFP[:,:end]


    def get_probable_detections(maxfrac=.5, rtol=.02):
        '''planets frequently flagged as an FP are probable actual planet 
        detections in the light curve but are flagged as a FP because they 
        are not injected by the algorithm. They should be flagged as a planet 
        detection and removed from the list of FPs to improve the accuracy of
        the ensuing FP correction.'''
        # probably shouldnt do this because multiple FPs could also be due
        # to a poor K2SFF systematics correction which means that it must be
        # included in the FP correction
        assert False
        assert 0 < maxfrac < 1
        assert rtol > 0
        P_FP = self.Psfound[self.isFP == 1]
        Pdet = np.zeros(0)
        for i in range(P_FP.size):
            
            N_at_this_P = float(np.isclose(P_FP[i], P_FP, rtol=rtol).sum())

            # save this planet if it is 'falsely' detected in more than
            # maxfrac of the simulations and correct the number of FPs
            if N_at_this_P / self.Nsim > maxfrac:
                Pdet = np.append(Pdet, P_FP[i])
                g = np.isclose(P_FP[i], self.Psfound, rtol=rtol)
                self.isFP[g] = 0


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
		
		g = (self.Psfound >= self.Pgrid[i]) & \
                    (self.Psfound <= self.Pgrid[i+1]) & \
                    (self.rpsfound >= self.rpgrid[j]) & \
                    (self.rpsfound <= self.rpgrid[j+1])
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


    def plot_map(self, zmap, zlabel='', xarr=np.zeros(0), yarr=np.zeros(0), avgtitle=False, sumtitle=False, issens=False, pltt=True, label=False):
	fig = plt.figure()
	ax = fig.add_subplot(111)
	xarr = self.Pgrid if xarr.size == 0 else xarr
        yarr = self.rpgrid if yarr.size == 0 else yarr
	img = ax.pcolormesh(xarr, yarr, zmap.T, cmap=plt.get_cmap('hot_r'))
	cbar_axes = fig.add_axes([.1,.1,.87,.04])
	cbar = fig.colorbar(img, cax=cbar_axes, orientation='horizontal')
        cbar.set_label(zlabel)
	ax.set_xscale('log')
	ax.set_xlabel('Period [days]'), ax.set_ylabel('Planet Radius [R$_{\oplus}$]')
	if sumtitle:
	    ax.set_title('Total = %i'%np.nansum(zmap), fontsize=12)
	if avgtitle:
	    ax.set_title('Average = %.3f'%np.nanmean(zmap), fontsize=12)

	# fill nans
	g = np.where(np.isnan(zmap))
	for i in range(g[0].size):
	    x1, x2 = xarr[g[0][i]], xarr[g[0][i]+1]
	    y1, y2 = yarr[g[1][i]], yarr[g[1][i]+1]
	    ax.fill([x1,x2,x2,x1], [y1,y1,y2,y2], fill=False, hatch='\\')

	# fill low sens if plotting a sensitivity map
	if issens:
	    g = np.where(zmap<.15)
            for i in range(g[0].size):
            	x1, x2 = xarr[g[0][i]], xarr[g[0][i]+1]
            	y1, y2 = yarr[g[1][i]], yarr[g[1][i]+1]
            	ax.fill([x1,x2,x2,x1], [y1,y1,y2,y2], fill=False, hatch='//')

	fig.subplots_adjust(bottom=.24, top=.95)
	if pltt:
	    plt.show()
	plt.close('all')

                          
    def _pickleobject(self):
        fObj = open(self.fname_full, 'wb')
        pickle.dump(self, fObj)
        fObj.close()


if __name__ == '__main__':
    fs = np.array(glob.glob('PipelineResults/EPIC_*'))
    print fs.size
    for i in range(fs.size):
	print fs[i]
	epicnum = int(fs[i].split('_')[-1])
    	self = K2sensitivity(epicnum)
