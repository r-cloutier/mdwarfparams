from LCclass import *
import rvs
from uncertainties import unumpy as unp

class Sensitivityclass:

    def __init__(self, epicnum, prefix):
        self.epicnum, self.prefix = epicnum, prefix
	self.prefix2 = 'EPIC' if prefix == 'K2' else 'KIC'
        self.fname_full = 'PipelineResults/%s_%i/%s_%i_sens'%(self.prefix2,
							      self.epicnum,
							      self.prefix2,
                                                              self.epicnum)
        self.get_data()
	if self.fs.size > 0:
            #self.get_probable_detections()
            self.compute_sensitivity()
	    self.compute_transit_prob()
        self._pickleobject()


    def get_data(self):
        self.fs = np.array(glob.glob('PipelineResults/%s_%i/*LC*'%(self.prefix2,self.epicnum)))
	# remove planet search result (i.e. with index -99)
	if np.any(np.in1d(self.fs, 'PipelineResults/%s_%i/%sLC_-00099'%(self.prefix2,
									self.epicnum,self.prefix))):
	    g = np.where(np.in1d(self.fs, 'PipelineResults/%s_%i/%sLC_-00099'%(self.prefix2,self.epicnum,
									       self.prefix))[0][0]
	    self.fs = np.delete(self.fs, g)
        if self.fs.size == 0:
	    return None
        self.Nsim = 0
        d = loadpickle(self.fs[0])
        self.Kepmag, self.logg, self.Ms, self.Rs, self.Teff = d.Kepmag, \
                                                              d.logg, \
                                                              d.Ms, d.Rs, \
                                                              d.Teff
        self.e_logg, self.e_Ms, self.e_Rs, self.e_Teff = d.e_logg, \
                                                         d.e_Ms, d.e_Rs, \
                                                         d.e_Teff
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


    def compute_sensitivity(self, xlen=120, ylen=60, Plims=(.5,80),
                            rplims=(.5,10)):
        '''Get all the simulations for this star and compute the
        sensitivity and the number of FPs as functions of P and rp.'''
        self._xlen, self._ylen = int(xlen), int(ylen)
        self.logPgrid = np.logspace(np.log10(Plims[0]), np.log10(Plims[1]),
                                    self._xlen+1)
        self.rpgrid = np.linspace(rplims[0], rplims[1], self._ylen+1)
        self.logrpgrid = np.logspace(np.log10(rplims[0]), np.log10(rplims[1]),
                                     self._ylen+1)
        self.Ndet, self.Ntrue, self.NFP = np.zeros((self._xlen, self._ylen)), \
                                          np.zeros((self._xlen, self._ylen)), \
                                          np.zeros((self._xlen, self._ylen))
        self.logNdet, self.logNtrue, self.logNFP = np.zeros_like(self.Ndet), \
                                                   np.zeros_like(self.Ntrue), \
                                                   np.zeros_like(self.NFP)
        for i in range(self._xlen):
            for j in range(self._ylen):
                # get detections in log-linear space
                g = (self.Ps >= self.logPgrid[i]) & \
                    (self.Ps <= self.logPgrid[i+1]) & \
                    (self.rps >= self.rpgrid[j]) & \
                    (self.rps <= self.rpgrid[j+1])
                self.Ndet[i,j] = self.isdet[g].sum()
                self.Ntrue[i,j] = self.isdet[g].size

		# get FPs in log-linear space
		g = (self.Psfound >= self.logPgrid[i]) & \
                    (self.Psfound <= self.logPgrid[i+1]) & \
                    (self.rpsfound >= self.rpgrid[j]) & \
                    (self.rpsfound <= self.rpgrid[j+1])
                self.NFP[i,j] = self.isFP[g].sum()
                
                # get detections in log-log space
                g = (self.Ps >= self.logPgrid[i]) & \
                    (self.Ps <= self.logPgrid[i+1]) & \
                    (self.rps >= self.logrpgrid[j]) & \
                    (self.rps <= self.logrpgrid[j+1])
                self.logNdet[i,j] = self.isdet[g].sum()
                self.logNtrue[i,j] = self.isdet[g].size

		# get FPs in log-log space
		g = (self.Psfound >= self.logPgrid[i]) & \
                    (self.Psfound <= self.logPgrid[i+1]) & \
                    (self.rpsfound >= self.logrpgrid[j]) & \
                    (self.rpsfound <= self.logrpgrid[j+1])
                self.logNFP[i,j] = self.isFP[g].sum()

        # compute sensitivity
        self.sens = self.Ndet / self.Ntrue.astype(float)
        self.esens = np.sqrt(self.Ndet) / self.Ntrue.astype(float)
        self.logsens = self.logNdet / self.logNtrue.astype(float)
        self.logesens = np.sqrt(self.logNdet) / self.logNtrue.astype(float)

        # compute yield correction to multiply the yield by
        self.yield_corr = 1 - self.NFP / (self.Ndet + self.NFP.astype(float))
        self.logyield_corr = 1 - self.logNFP / \
                             (self.logNdet + self.logNFP.astype(float))
        
    
    def compute_transit_prob(self):
	self.transit_prob, self.etransit_prob = np.zeros_like(self.sens), \
                                                np.zeros_like(self.sens)
	self.logtransit_prob,self.elogtransit_prob = np.zeros_like(self.sens), \
                                                     np.zeros_like(self.sens)
	for i in range(self._xlen):
	    for j in range(self._ylen):

                # compute transit probability in log-linear space 
                Pmid = 10**(np.log10(self.logPgrid[i]) + \
                            np.diff(np.log10(self.logPgrid[:2])/2))
		rpmid = self.rpgrid[j] + np.diff(self.rpgrid)[0]/2.
		sma = rvs.AU2m(rvs.semimajoraxis(Pmid, unp.uarray(self.Ms,
                                                                  self.e_Ms),
                                                 0))
                transit_prob =  (rvs.Rsun2m(unp.uarray(self.Rs, self.e_Rs)) + \
                                 rvs.Rearth2m(rpmid)) / sma
                self.transit_prob[i,j] = unp.nominal_values(transit_prob)
                self.etransit_prob[i,j] = unp.std_devs(transit_prob)
                
                # compute transit probability in log-linear space 
                rpmid = 10**(np.log10(self.logrpgrid[j]) + \
                             np.diff(np.log10(self.logrpgrid[:2])/2))
		transit_prob =  (rvs.Rsun2m(unp.uarray(self.Rs, self.e_Rs)) + \
                                 rvs.Rearth2m(rpmid)) / sma
                self.logtransit_prob[i,j] = unp.nominal_values(transit_prob)
                self.elogtransit_prob[i,j] = unp.std_devs(transit_prob)

	# correction from beta distribution fit (Kipping 2013)
	factor = 1.08
	self.transit_prob *= factor
	self.etransit_prob *= factor
        self.logtransit_prob *= factor
	self.elogtransit_prob *= factor


    def plot_map(self, xarr, yarr, zmap, zlabel='', avgtitle=False,
                 sumtitle=False, issens=False, pltt=True, label=False):
	fig = plt.figure()
	ax = fig.add_subplot(111)
	img = ax.pcolormesh(xarr, yarr, zmap.T, cmap=plt.get_cmap('hot_r'))
	cbar_axes = fig.add_axes([.1,.1,.87,.04])
	cbar = fig.colorbar(img, cax=cbar_axes, orientation='horizontal')
        cbar.set_label(zlabel)
	if not np.all(np.isclose(np.diff(xarr), np.diff(xarr)[0])):
            ax.set_xscale('log')
        if not np.all(np.isclose(np.diff(yarr), np.diff(yarr)[0])):
            ax.set_yscale('log')

        ax.set_xlabel('Period [days]')
        ax.set_ylabel('Planet Radius [R$_{\oplus}$]')
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


        
class SensitivityFULL:

    def __init__(self, folder, prefix, xlen=120, ylen=60):
        self.folder, self.prefix = folder, prefix
	self.prefix2 = 'EPIC' if prefix == 'K2' else 'KIC'
        self._xlen, self._ylen = int(xlen), int(ylen)
	self.fs = np.array(glob.glob('%s/%s_*/%s_*_sens'%(self.folder,
							  self.prefix2,
							  self.prefix2)))
	self.Nstars = self.fs.size
        self.fname_full = '%s/%s_%ssens'%(self.folder,self.prefix2,self.prefix)
        self.get_full_maps()
        self._pickleobject()


    def get_full_maps(self):
        '''combine the maps from the injection/recovery simulations on 
        individual systems to get the master maps.'''
        # firstly, get the number of simulations run on each systems
        # for weighting
        self.Nsims = np.array([loadpickle(self.fs[i]).Nsim
                               for i in range(self.Nstars)])
        norm = float(self.Nsims.sum())
        
        # combine maps
        d = loadpickle(self.fs[0])
        self.logPgrid, self.logrpgrid = d.logPgrid, d.logrpgrid
        self.sensFULL         = np.zeros((self._xlen, self._ylen))
        self.yield_corrFULL   = np.zeros_like(self.sensFULL)
        self.transit_probFULL = np.zeros_like(self.sensFULL)

        for i in range(self.Nstars):
            print float(i) / self.Nstars
            d = loadpickle(self.fs[i])

            # compute bin the maps if necessary
            if (d.sens.shape[0] != self._xlen) or \
               (d.sens.shape[1] != self._ylen):
                d.compute_sensitivity(xlen=self._xlen, ylen=self._ylen)
                d.compute_transit_prob()
                self.logPgrid, self.logrpgrid = d.logPgrid, d.logrpgrid

            weight = d.Nsim / norm
            self.sensFULL         += weight * zero_nans(d.logsens)
            self.yield_corrFULL   += weight * zero_nans(d.logyield_corr)
            self.transit_probFULL += weight * zero_nans(d.logtransit_prob)


    def _pickleobject(self):
        fObj = open(self.fname_full, 'wb')
        pickle.dump(self, fObj)
        fObj.close()

        

def zero_nans(arr):
    arrv2 = np.copy(arr)
    g = (np.isnan(arr)) | (np.isinf(arr))
    arrv2[g] = 0.
    return arrv2



if __name__ == '__main__':
    epicnums = np.loadtxt('input_data/K2targets/K2Mdwarfs_withdetections.csv',
			  delimiter=',')
    for i in range(epicnums.size):
	print epicnums[i]
    	self = Sensitivityclass(epicnums[i])

    xlen, ylen = 14, 8
    self = SensitivityFULLclass('PipelineResults', xlen=xlen,
                             	ylen=ylen)
