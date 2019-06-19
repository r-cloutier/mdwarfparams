from imports import *
from linear_lnlike import boxcar
from occurrencerateclass import *


def get_tics_with_candidates(self):
    '''Get TICs with a detected planet candidate for use in copy_detections.py 
    on bubbles.'''
    N = self.fs_planetsearch.size
    ticsout=np.array([int(self.fs_planetsearch[i].split('/')[1].split('_')[1])
                      for i in range(N)])
    g = self.Ps > 0
    return np.unique(ticsout[g])


def download_LCs(tics, folder='PipelineResults_TIC_sector34'):
    '''After running copy_detections.py on bubbles which copies light curves to
    where they can be accessed, download and save the light curves.'''
    url = 'https://www.cita.utoronto.ca/~cloutier'
    for i in range(tics.size):
	
	# make directory
	try:
	    os.mkdir('%s/EPIC_%i'%(folder, tics[i]))
	except OSError:
	    pass

	# download light curve
	os.system('wget %s/LC%i_-00099'%(url, tics[i]))
	
	# save light curve in the proper directory
	os.system('mv LC%i_-00099 %s/EPIC_%i/LC_-00099'%(tics[i],folder,tics[i]))


def get_tics(self):
    '''Create the array of tics which is often lacking in the detection 
    object.'''
    if not hasattr(self, 'tics'):
        N = self.fs_planetsearch.size
        self.tics = np.zeros(N, dtype=int)
        for i in range(N):
            self.tics[i] = int(self.fs_planetsearch[i].split('/')[1].split('_')[1])
    self._pickleobject() 

        
def diagnostic_plots(self, tic, sctr, folder='PipelineResults_TIC_sector34',
                     pltts=[0,1,2]):
    '''Load a light curve with a detected planet candidate and flag each planet
    candidate according to its user-defined disposition based on the plots 
    created in this function.'''
    # get data
    d = loadpickle('%s/TIC_%i/LC_-00099'%(folder, tic))
    assert d.Ndet >= 1
    assert d.tic == tic
    print '\nTIC %i'%tic

    # get data folder
    tic_str = '%.16d'%int(tic)
    tid1 = tic_str[:4]
    tid2 = tic_str[4:8]
    tid3 = tic_str[8:12]
    tid4 = tic_str[12:]
    url = 'https://archive.stsci.edu/missions/tess/tid/'
    folder = 's%.4d/%s/%s/%s/%s/'%(sctr,tid1,tid2,tid3,tid4)
    print url+folder   
 
    # plot full light curve
    t0 = 2457000
    if 0 in pltts:
        plt.plot(d.bjd-t0, d.f, 'k.', alpha=.6, label='raw flux')
        plt.plot(d.tbin-t0, d.fbin, 'g.', label='binned raw flux')
        plt.plot(d.bjd-t0, d.mu, 'b-', label='mean GP')
        plt.xlabel('Time [BJD - 2,457,000]')
        plt.ylabel('Normalized flux')
        plt.title('Raw TESS light curve')
        plt.legend(loc='upper left')
        plt.show()

    # plot detrended flux
    cols = ['b','g','r','c','y','o','v']
    if 1 in pltts:
        plt.plot(d.bjd-t0, d.fcorr, 'k.', alpha=.5,
                 label='detrended flux')
        plt.xlim((d.bjd.min()-t0-.3, d.bjd.max()-t0+.3))
        plt.xlabel('Time [BJD - 2,457,000]')
        plt.ylabel('Normalized flux')
        plt.title('Detrended TESS light curve')
        # mark transits
        for i in range(d.Ndet):
            P,T0 = d.params_optimized[i,:2]
            for j in range(-20,20):
                if j == -20:
                    plt.axvline(T0+P*j-t0, ymax=.05, lw=2, color=cols[i],
                                label='P%i = %.3f days'%(i+1,P))
                else:
                    plt.axvline(T0+P*j-t0, ymax=.05, lw=2, color=cols[i])
                
        plt.legend(loc='upper left', fontsize=11)
        plt.show()

    # plot phase-folded light curves
    if 2 in pltts:
        for i in range(d.Ndet):
            P,T0 = d.params_optimized[i,:2]
            phase = foldAt(d.bjd, P, T0)
            phase[phase>.5] -= 1
            s = np.argsort(phase)
            # plot transit
            plt.figure(figsize=(12,5))
            plt.subplot(121)
            plt.plot(phase, d.fcorr, 'k.', alpha=.5)
            plt.plot(phase[s], d.fmodels[i,s], '%s-'%cols[i])
            plt.xlabel('Phase (P%i = %.3f days)'%(i+1,P))
            plt.ylabel('Normalized flux')

            # zoom-in
            plt.subplot(122)
            plt.plot(phase, d.fcorr, 'k.', alpha=.7)
            plt.plot(phase[s], d.fmodels[i,s], '%s-'%cols[i])
            tbin, fbin, efbin = boxcar(phase[s], d.fcorr[s], d.ef[s], dt=.001)
            plt.plot(tbin, fbin, 'y.')
            plt.xlabel('Phase (P%i = %.3f days)'%(i+1,P))
            D = d.params_guess[i,3]
            plt.xlim((-5*D/P,5*D/P))
            
            plt.show()

    # ask for user-defined dispositions
    disposition_labels = {-1: 'FP',
                          -2: 'FP due to bad pointing',
                          -3: 'EB',
                          -4: 'edge effect',
                          -6: 'ruled out as BEB',
                          -7: 'ruled out as BEB with double period',
                          -8: 'undetermined FP',
                          0: 'maybe a planet',
                          1: 'PC',
                          2: 'possible single transit',
                          2.5: 'putative ST',
                          5: 'wtf'}
    for key, value in sorted(disposition_labels.iteritems(),
                             key=lambda (k,v): (v,k)):
        print '%s: %s'%(key, value)

    # create arrays if not already
    if not hasattr(self, 'tics'):
        get_tics(self)
    if not hasattr(self, 'disposition_human'):
        self.disposition_human = np.zeros(self.tics.size) + np.nan

    # get human dispositions based on plots
    for i in range(d.Ndet):
        P = d.params_optimized[i,0]
        g = (self.tics == tic) & (self.Ps == P)
        assert g.sum() == 1
        prompt = 'Disposition on planet %i (%.3f days) '%(i+1, P)
        self.disposition_human[g] = int(raw_input(prompt))
