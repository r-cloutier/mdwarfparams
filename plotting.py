from GAIAMdwarfs import *
import matplotlib.gridspec as gridspec
import matplotlib
import string
from occurrencerateclass import *
import linear_lnlike as llnl
from TESS_search import *
from matplotlib.ticker import NullFormatter
import matplotlib.colors as colors
import matplotlib.patches as patches
from followup import *
from query_gaia_2MASS import query_nearby_gaia


matplotlib.rcParams.update({'ytick.labelsize': 9})
matplotlib.rcParams.update({'xtick.labelsize': 9})


def truncate_cmap(cmap, minval=0, maxval=1, n=100):
    return colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f},)'.format(n=cmap.name,a=minval,b=maxval),
        cmap(np.linspace(minval, maxval, n)))


def plot_distance_scatter(self, pltt=True, label=False):
    '''plot the NASA archive distances vs the GAIA distances along with the 
    uncertainties as a function of distances.'''
    g = (self.loggN>3) & (self.RsN<1) & (self.TeffN<4600) & \
        np.isfinite(self.parallaxG) & np.isfinite(self.distN) & \
        (self.parallax_reliable==1)
    _,indstmp = np.unique(self.names[g], return_index=True)
    g = np.arange(self.names.size)[g][indstmp]
    assert g.size == 163

    fig = plt.figure(figsize=(7,5))
    gs = gridspec.GridSpec(10,5)
    ax1 = plt.subplot(gs[2:8,:4])
    ax2 = plt.subplot(gs[:2,:4])
    ax3 = plt.subplot(gs[2:8,4:])
    ax4 = plt.subplot(gs[8:,:4])
    
    x, ux, lx, y, uy, ly = self.distN[g], self.udistN[g], self.ldistN[g], \
                           self.distG[g], self.edistG[g], self.edistG[g]
    ax1.errorbar(x, y, xerr=[lx,ux], yerr=[ly,uy], fmt='ko', elinewidth=.5,
                 ms=1)
    ax1.plot([.7,2e3],[.7,2e3], 'b--')
    ax1.set_xscale('log'), ax1.set_yscale('log')
    ax1.set_xlim((.7,2e3)), ax1.set_ylim((.7,2e3))
    ax1.set_xticklabels('')
    ax1.set_ylabel('GAIA Distance\n[pc]', fontsize=12)
    
    ax2.plot(x, np.mean([ux,lx],axis=0), 'ko', ms=2)
    ax2.set_xscale('log'), ax2.set_yscale('log')
    ax2.set_xlim((.7,2e3)), ax2.set_ylim((1e-4,1e3))
    ax2.set_xticklabels('')
    ax2.set_ylabel('Archive Distance\nUncertainty [pc]', fontsize=11)

    ax3.plot(np.mean([uy,ly],axis=0), y, 'ko', ms=2)
    ax3.set_xscale('log'), ax3.set_yscale('log')
    ax3.set_xlim((1e-4,1e3)), ax3.set_ylim((.7,2e3))
    ax3.set_yticklabels('')
    ax3.set_xlabel('GAIA\nDistance\nUncertainty\n[pc]', fontsize=11)

    uratio = unp.uarray(y, uy) / unp.uarray(x, ux)
    lratio = unp.uarray(y, ly) / unp.uarray(x, lx)
    ax4.errorbar(x, unp.nominal_values(uratio), xerr=[lx,ux],
                 yerr=[unp.std_devs(lratio), unp.std_devs(uratio)], fmt='k.',
                 elinewidth=.8)
    ax4.plot([.7,2e3],np.ones(2), 'b--')
    ax4.set_xscale('log'), ax4.set_yscale('log'), ax4.set_xlim((.7,2e3))
    ax4.set_ylim((.25,4))
    ax4.set_yticks([1./3,1,3])
    ax4.set_yticklabels(['1/3','1','3'], fontsize=11)
    ax4.set_xlabel('Archive Distance [pc]', fontsize=12)
    ax4.set_ylabel('GAIA Distance /\nArchive Distance', fontsize=11)

    fig.subplots_adjust(bottom=.09, top=.94, wspace=0, hspace=0)
    if label:
        plt.savefig('plots/distance_scatter.png')
    if pltt:
        plt.show()
    plt.close('all')


    
def plot_stellar_radius_uncertainty(pltt=True, label=False):
    # get initial sample
    f1 = fits.open('input_data/Keplertargets/kepler_dr2_1arcsec.fits')[1]
    f4 = fits.open('input_data/Keplertargets/kepler_dr2_4arcsec.fits')[1]
    kepids = np.append(f1.data['kepid'], f4.data['kepid'])
    _,g = np.unique(kepids, return_index=True)
    Teff = np.append(f1.data['teff'], f4.data['teff'])[g]
    eTeff = np.append(f1.data['teff_err2'], f4.data['teff_err2'])[g]
    logg = np.append(f1.data['logg'], f4.data['logg'])[g]
    elogg = np.append(f1.data['logg_err1'], f4.data['logg_err1'])[g]
    Rs = np.append(f1.data['radius'], f4.data['radius'])[g]
    eRs = np.append(f1.data['radius_err1'], f4.data['radius_err1'])[g]
    
    g = (Teff-eTeff <= 4e3) & (logg+elogg > 3.5) & (Rs-eRs < .75)
    Rs, eRs = Rs[g], eRs[g]
    ratio = eRs / Rs
    print np.nanmean(ratio)

    fig = plt.figure(figsize=(8,6))
    ax = fig.add_subplot(111)
    #ax.hist(ratio, bins=np.logspace(np.log10(.01), 0, 30))
    ax.plot(Rs, ratio, 'ko', ms=2, alpha=.3,
            label='Initial stellar sample (%i stars)'%ratio.size)
    ax.set_yscale('log'), ax.set_xlim((0,.8))
    ax.set_yticks([.03,.06,.1,.3,.6,1])
    ax.set_yticklabels(['0.03','0.06','0.1','0.3','0.6','1'])
    ax.set_xlabel('Stellar Radius, R$_s$')
    ax.set_ylabel('Fractional Stellar Radius\nUncertainty, $\sigma_{R_s}$/R$_s$')

    ax.legend()
    if label:
        plt.savefig('plots/sigRs.png')
    if pltt:
        plt.show()
    plt.close('all')


def plot_detected_population_Prp(occurrencerateclass, pltt=True, label=False):
    '''Plot the empirical distribution of detected planets in the P,rp plane.'''
    self = occurrencerateclass
    assert 'det' in self.fname_out

    sigrp = self.rps / np.mean([self.elo_rps,self.ehi_rps],0)
    g = np.isfinite(self.rps) & np.isfinite(self.Ps) & (sigrp > 3)
    
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.errorbar(self.Ps[g], self.rps[g], xerr=self.e_Ps[g],
                yerr=[self.elo_rps[g],self.ehi_rps[g]], fmt='k.',
                elinewidth=.9)
    ax.set_xscale('log'), ax.set_yscale('log')
    ax.set_xlabel('Period [days]')
    ax.set_ylabel('Planet Radius [R$_{\oplus}$]')
    ax.set_title('%i planet candidates'%self.Ps[g].size, fontsize=12)
    
    if label:
        plt.savefig('plots/Ndet_Prp.png')    
    if pltt:
        plt.show()
    plt.close('all')


def plot_detected_population_rphist(occurrencerateclass, pltt=True,
                                    label=False):
    '''Plot the rp histogram of detected planets.'''
    self = occurrencerateclass
    assert 'det' in self.fname_out

    sigrp = self.rps / np.mean([self.elo_rps,self.ehi_rps],0)
    g = np.isfinite(self.rps) & np.isfinite(self.Ps) & (sigrp > 3)
    g = (g) & (self.Mss > .6)

    rpbins = np.logspace(np.log10(.5), np.log10(10), 31)
    N, rp_edges = np.histogram(self.rps[g], bins=rpbins)
    rp = 10**(np.log10(rp_edges[1:]) - np.diff(np.log10(rp_edges))[0]/2)

    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.errorbar(rp, N, yerr=np.sqrt(N), color='k', elinewidth=.9,
                drawstyle='steps-mid', capsize=2, lw=2.5)
    ax.set_xscale('log')
    ax.set_xlabel('Planet Radius [R$_{\oplus}$]')
    ax.set_ylabel('Number of planet canadidates')
    ax.set_title('%i planet candidates'%N.sum(), fontsize=12)
    title = 'all M$_s$'# <$ 0.6 M$_{\odot}$'
    ax.text(.8, .85, title, fontsize=11, transform=ax.transAxes)
    
    if label:
        plt.savefig('plots/Ndet_rphist_all.png')
    if pltt:
        plt.show()
    plt.close('all')



def plot_1_planet_candidate_summary(self, rp, e_rp, F, e_F, SNRoffset=10,
                                    phaseylims=(.99,1.01), fluxlim=(.995,1.005),
                                    phasexlims=(-.1,.1), pltt=True,
                                    label=False):
    assert self.Ndet == 1

    # setup stuff
    cols = ['k','#ff4700','#fe9500','#a80000','#ffff5f','#fe9500']
    cols_linear = ['k','#a80000','#ff4700']

    Z = self.params_guess[0,2]
    P, T0, aRs, rpRs, inc = self.params_optimized[0]
    rp = rvs.m2Rearth(rvs.Rsun2m(rpRs*self.Rs))
    b = rvs.impactparam_inc(P, self.Ms, self.Rs, inc)
    D = rvs.transit_width(P, self.Ms, self.Rs, rp, b)
    title_nums = self.Ndet, np.round(self.TESSmag,3), \
                 np.round(self.Teff), np.round(self.ehi_Teff), \
                 np.round(self.logg,3), np.round(self.ehi_logg,3), \
                 np.round(self.Rs,3), np.round(self.ehi_Rs,3)
    title = r'$\bf{'+'OI'+'}$'+' '+r'$\bf{'+str(self.tic)+'.01}$'+' (%i planet candidate)\nT$_{mag}$ = %.3f,\tT$_{eff}$ = %i $\pm$ %i K,\t$\log{g}$ = %.3f $\pm$ %.3f,\tR$_s$ = %.3f $\pm$ %.3f R$_{\odot}$'%title_nums
    
    # plot the extracted 2 minute pdcsap light curve
    fig = plt.figure(figsize=(6.5,8.5))
    #gs = gridspec.GridSpec(10,5)
    ax1 = fig.add_axes([.13,.83,.83,.12])
    t0 = 2457000
    bjdlim = self.bjd.min()-t0-.3, self.bjd.max()-t0+.3
    ax1.plot(self.bjd-t0, self.f, 'o', c=cols[0], ms=1, alpha=.5)
    ax1.fill_between(self.bjd-t0, self.mu-self.sig, self.mu+self.sig, alpha=.3,
                     color=cols[1])
    ax1.plot(self.bjd-t0, self.mu, '-', c=cols[1], lw=.9)
    for i in range(-50,50):
        ax1.axvline(T0-t0+i*P, ymax=.08, lw=2, color=cols[2])
        ax1.axvline(T0-t0+i*P, ymin=.92, lw=2, color=cols[2])
    ax1.set_xlim(bjdlim)
    ax1.set_ylim(fluxlim)
    ax1.set_xticklabels('')
    ax1.set_ylabel('Normalized PDCSAP\nextracted flux', fontsize=9)
    ax1.set_title(title, fontsize=9, y=1)

    # plot the detrended light curve
    ax2 = fig.add_axes([.13,.7,.83,.12])
    ax2.plot(self.bjd-t0, self.fcorr, 'o', c=cols[0], ms=1, alpha=.5)
    for i in range(-50,50):
        ax2.axvline(T0-t0+i*P, ymax=.08, lw=2, color=cols[2])
        ax2.axvline(T0-t0+i*P, ymin=.92, lw=2, color=cols[2])
    ax2.set_xlim(bjdlim)
    ax2.set_ylim(fluxlim)
    ax2.set_xticklabels('')
    ax2.set_ylabel('Normalized\nde-trended flux', fontsize=9)

    # plot linear transit search
    ax3 = fig.add_axes([.13,.45,.83,.24])
    for i in range(3):
        ax3.plot(self.transit_times-t0, self.SNRs_linearsearch[:,i]-i*SNRoffset,
                 '-', lw=.9, c=cols_linear[i])
        ax3.axhline(-i*SNRoffset+5, lw=.9, ls='--', c=cols_linear[i])
        ax3.text(self.bjd.min()-t0+.2, -i*SNRoffset+3.1,
                 'D = %.1f hrs'%(self.durations[i]*24), fontsize=7,
                 color=cols_linear[i], weight='semibold')
    for i in range(-50,50):
        ax3.axvline(T0-t0+i*P, ymax=.04, lw=2, color=cols[2])
        ax3.axvline(T0-t0+i*P, ymin=.96, lw=2, color=cols[2])
    ax3.set_xlim(bjdlim)
    ax3.set_xlabel('Time [BJD - 2,457,000]', fontsize=9)
    ax3.set_ylabel('Linear search S/N =\n($\ln{\mathcal{L}}$ - median($\ln{\mathcal{L}}$)) / MAD($\ln{\mathcal{L}}$)', fontsize=9)
    
    # plot phase-folded light curve
    phase = foldAt(self.bjd, P, T0)
    phase[phase > .5] -= 1
    s = np.argsort(phase)
    tb, fb, efb = llnl.boxcar(phase[s], self.fcorr[s], self.ef[s], dt=.2*D/P)
    ax4 = fig.add_axes([.13,.24,.45,.15])
    ax4.plot(phase, self.fcorr, 'o', c=cols[0], ms=1, alpha=.5)
    ax4.plot(tb, fb, '-', c=cols[2], lw=1)
    ax4.set_xlim((-.5,.5))
    ax4.set_ylim(phaseylims)
    ax4.set_ylabel('Normalized\nde-trended flux', fontsize=9)

    # zoom-in on transit
    ax5 = fig.add_axes([.13,.05,.45,.15])
    ax5.plot(phase, self.fcorr, 'o', c=cols[0], ms=1, alpha=.3)
    ax5.errorbar(tb, fb, efb, fmt='o', markerfacecolor=cols[2], ms=5,
                 markeredgecolor='k', markeredgewidth=.5, elinewidth=1)
    ax5.plot(phase[s], self.fmodels[0][s], '-', c=cols[0], lw=2)
    ax5.set_xlim(phasexlims)
    ax5.set_ylim(phaseylims)
    ax5.set_xlabel('Orbital Phase', fontsize=9)
    ax5.set_ylabel('Normalized\nde-trended flux', fontsize=9)

    # data info
    ax4.text(1.48, .96, 'Reduced light curve parameters', fontsize=8,
             weight='semibold', transform=ax4.transAxes,
             horizontalalignment='center')
    rms = np.std(self.fcorr)*1e6
    ax4.text(1.05, .84, 'Total LC rms = %i [ppm]'%rms, fontsize=8,
             transform=ax4.transAxes)
    cdpp = np.round(self.CDPPs[0]*1e6)
    D = np.round(self.params_guess[0,3]*24,1)
    ax4.text(1.05, .72, 'CDPP$_{transit}$ = %i [ppm] (D = %.1f hrs)'%(cdpp,D),
             fontsize=8, transform=ax4.transAxes)
    ax4.text(1.05, .6, 'S/N$_{transit}$ = %.1f'%(np.round(self.SNRtransits[0],
                                                          1)), fontsize=8,
             transform=ax4.transAxes)

    # GP parameters
    ax4.text(1.48, .46, 'Systematic GP model parameters', fontsize=8,
             weight='semibold', transform=ax4.transAxes,
             horizontalalignment='center')
    ax4.text(1.05, .34, '$\ln$a$_{GP}$ = %.3f'%(np.round(self.thetaGPout[0,0],3)),
             fontsize=8, transform=ax4.transAxes)
    ax4.text(1.05, .22, '$\ln{\lambda}$/days = %.3f'%(np.round(self.thetaGPout[0,1],3)), fontsize=8, transform=ax4.transAxes)
    ax4.text(1.05, .1, '$\ln{\Gamma}$ = %.3f'%(np.round(self.thetaGPout[0,2],3)), fontsize=8, transform=ax4.transAxes)
    ax4.text(1.05, -.02, '$\ln$P$_{GP}$/days = %.3f'%(np.round(self.thetaGPout[0,3],3)), fontsize=8, transform=ax4.transAxes)
             
    # planet info
    e_P,e_T0,e_aRs,e_rpRs,e_inc = np.mean(self.params_results[0,1:],0)
    ax5.text(1.48, 1.1, 'Measured planet parameters', fontsize=8,
             weight='semibold', transform=ax5.transAxes,
             horizontalalignment='center')
    ax5.text(1.05, .98, 'P = %.5f $\pm$ %.5f [days]'%(P,e_P), fontsize=8,
             transform=ax5.transAxes)
    ax5.text(1.05, .86, 'T$_0$ = %.5f $\pm$ %.5f [BJD-2,457,000]'%(T0-t0,e_T0),
             fontsize=8, transform=ax5.transAxes)
    ax5.text(1.05, .74, 'a/R$_s$ = %.2f $\pm$ %.2f'%(np.round(aRs,2),
                                                     np.round(e_aRs,2)),
             fontsize=8, transform=ax5.transAxes)
    ax5.text(1.05, .62, 'rp/R$_s$ = %.4f $\pm$ %.4f'%(np.round(rpRs,4),
                                                      np.round(e_rpRs,4)),
             fontsize=8, transform=ax5.transAxes)
    ax5.text(1.05, .5, 'i = %.2f $\pm$ %.2f [deg]'%(np.round(inc,2),
                                                    np.round(e_inc,2)),
             fontsize=8, transform=ax5.transAxes)
    ax5.text(1.05, .38, 'a$_{LDC}$ = %.3f (fixed)'%np.round(self.u1,3),
             fontsize=8, transform=ax5.transAxes)
    ax5.text(1.05, .26, 'b$_{LDC}$ = %.3f (fixed)'%np.round(self.u2,3),
             fontsize=8, transform=ax5.transAxes)
    #e_rp = unp.std_devs(rvs.m2Rearth(rvs.Rsun2m(unp.uarray(rpRs,e_rpRs) * \
    #                                        unp.uarray(self.Rs,self.ehi_Rs))))
    ax5.text(1.48, .12, 'Derived parameters', fontsize=8, weight='semibold',
             transform=ax5.transAxes, horizontalalignment='center')
    ax5.text(1.1, 0,'r$_p$ = %.2f $\pm$ %.2f [R$_{\oplus}$]'%(np.round(rp,2),
                                                        np.round(e_rp,2)),
             fontsize=8, transform=ax5.transAxes)
    #sma =rvs.semimajoraxis(unp.uarray(P,e_P),unp.uarray(self.Ms,self.ehi_Ms),0)
    #Teq = unp.uarray(self.Teff,self.ehi_Teff) * \
    #      unp.sqrt(rvs.Rsun2m(unp.uarray(self.Rs,self.ehi_Rs)) / \
    #               rvs.AU2m((2*sma))) * (1-.3)**(.25)
    #Teq, e_Teq = unp.nominal_values(Teq), unp.std_devs(Teq)
    ax5.text(1.1,-.12,'S = %.1f $\pm$ %.1f [S$_{\oplus}]$'%(np.round(F),
                                                            np.round(e_F)),
             fontsize=8, transform=ax5.transAxes)

    #fig.subplots_adjust(top=.93, right=.96, hspace=.2)
    if label:
        plt.savefig('plots/planetcandidatesummary_1_%i.png'%self.tic)
    if pltt:
        plt.show()
    plt.close('all')



def plot_linearsearch(self, dur_ind=0, P=0, T0=0, pltt=True, label=''):

    # black, red, yellow, orange
    cols = ['k','#a80000','#ffff5f','#ff4700']
    t0 = 2457000
    dur_ind = int(dur_ind)
    bjdlim = self.bjd.min()-t0-.3, self.bjd.max()-t0+.3
    gapts = 1338., 1339.66
    g=(self.transit_times-t0 <= gapts[0]) | (self.transit_times-t0 >= gapts[1])
    
    fig = plt.figure(figsize=(3.7,2.7))
    ax = fig.add_subplot(111)    
    ax.plot(self.transit_times[g]-t0, self.SNRs_linearsearch[g,dur_ind], '-',
            c=cols[dur_ind], lw=.4)
    ax.axhline(5, ls='--', color=cols[dur_ind], lw=1.5)
    ax.set_xlim(bjdlim)
    ax.set_ylim((-6,10))
    ax.set_xlabel('T$_0$ [BJD - 2,457,000]', fontsize=9)
    ax.set_ylabel('Linear search S/N =\n($\ln{\mathcal{L}}$ - median($\ln{\mathcal{L}}$)) / MAD($\ln{\mathcal{L}}$)', fontsize=9)
    #ax.text(1346, 5.5, 'S/N$_{LS}$ = 5', fontsize=7)
    ax.set_title('D = %.1f hrs'%(self.durations[dur_ind]*24), fontsize=9, y=.98)

    # add periodic markers
    P,T0 = self.params_optimized[0,:2]
    #fint = interp1d(self.transit_times, self.SNRs_linearsearch[:,dur_ind])
    t0s, snrs = np.zeros(0), np.zeros(0)
    for i in range(-30,30):
        t = T0+i*P
        if self.bjd.min() <= t <= self.bjd.max():
            t0s = np.append(t0s, t)
            g = abs(self.transit_times-t) < .05
            snr = self.SNRs_linearsearch[g,dur_ind].max() if g.sum() > 0 \
                  else np.nan
            snrs = np.append(snrs, snr)

    g1 = snrs < 5
    g2 = snrs >= 5
    ax.plot(t0s[g1]-t0, snrs[g1], 'o', markerfacecolor=cols[1], ms=4,
            markeredgecolor='k', markeredgewidth=.8)
    ax.plot(t0s[g2]-t0, snrs[g2], 'o', markerfacecolor=cols[2], ms=6,
            markeredgecolor='k', markeredgewidth=.9)

    # plot transit ticks
    if P != 0 and T0 != 0:
        for i in range(-27,27):
            if not (gapts[0] <= T0-t0+P*i <= gapts[1]):
                ax.axvline(T0-t0+P*i, ls='-', lw=.9, c=cols[3], ymin=.96)
                ax.axvline(T0-t0+P*i, ls='-', lw=.9, c=cols[3], ymax=.04)
    
    fig.subplots_adjust(bottom=.16, top=.92, right=.98, left=.16)
    if label:
        plt.savefig('plots/linearsearch_%i.png'%self.tic)
    if pltt:
        plt.show()
    plt.close('all')


# try: d=loadpickle('PipelineResults_TIC/TIC_234994474/LC_-00099')
def plot_periodicsearch(self, P, dur_ind=0, pltt=True, label=False):
    assert self.tic == 234994474
    cols = ['k','#e31a1c','#fed976','#fc4e2a']
    t0 = 2457000
    dur_ind = int(dur_ind)
    bjdlim = self.bjd.min()-t0-.3, self.bjd.max()-t0+.3
    gapts = 1338., 1339.66
    g=(self.transit_times-t0 <= gapts[0]) | (self.transit_times-t0 >= gapts[1])

    # get periodic markers
    dur_ind = int(dur_ind)
    transit_times = self.transit_times[self.SNRs_linearsearch[:,dur_ind] >= 5] 
    g = np.append(np.where(np.diff(transit_times)>P)[0], -1)
    transit_times = transit_times[g]
    N = transit_times.size

    # get lnLOIs of TIC 234994474 POIs
    POIs, lnLOIs = np.load('input_data/POIs_lnLOIs_tic234994474.npy').T
    POIs = np.append(POIs, [16.8,19.6])
    lnLOIs = np.append(lnLOIs, np.repeat(np.nanmedian(lnLOIs),2))
    
    # create period matrix from high S/N transit times and the lnL matrix
    indarr = np.arange(1,N+2)
    Pmat = np.zeros((N,N))
    lnL = np.zeros((N,N))
    for i in range(N):
        for j in range(N):
            Pmat[i,j] = abs(transit_times[i] - transit_times[j])
            g = np.isclose(POIs, abs(Pmat[i,j]), rtol=.1)
            if np.any(g):
                lnL[i,j] = lnLOIs[g][lnLOIs[g] == lnLOIs[g].max()]
            else:
                lnL[i,j] = np.nan
    lnL -= np.nanmedian(lnL)
                
    fig = plt.figure(figsize=(6.5,3.3))
    ax = fig.add_subplot(121)
    cax = ax.pcolormesh(Pmat.T,
                        cmap=truncate_cmap(plt.get_cmap('hot_r'),.01,.85))
    cbar_axes = fig.add_axes([.05,.1,.43,.04])
    cbar = fig.colorbar(cax, cax=cbar_axes, orientation='horizontal')
    cbar.set_label('P$_{i,j}$ = T$_{0,i}$-T$_{0,j}$ [days]', fontsize=8,
                   labelpad=.1)
    
    # plot grid
    for i in range(1,N):
        ax.axvline(i, ls='-', lw=.6, color='k')
        ax.axhline(i, ls='-', lw=.6, color='k')

    # highlight detected period
    for i in range(N):
        for j in range(N):
            col = 'k' if Pmat.T[i,j] < 14 else 'w'
            if i > j:
                rect = patches.Rectangle((i, j), 1, 1, color='k')
                ax.add_patch(rect)
            else:
                ax.text(i+.95, j+.68, '%.3f\ndays'%Pmat.T[i,j], fontsize=5,
                        color=col, weight='semibold',
                        verticalalignment='bottom',
                        horizontalalignment='right')
                ax.text(i+.95, j+.89, '(%.3f)'%(Pmat.T[i,j]/P), fontsize=5,
                        color=col, weight='semibold',
                        verticalalignment='bottom',
                        horizontalalignment='right')
            
    # plot lnL matrix
    ax1 = fig.add_subplot(122)
    cax1 = ax1.pcolormesh(lnL,cmap=truncate_cmap(plt.get_cmap('hot_r'),.01,.85))
    cbar_axes1 = fig.add_axes([.55,.1,.43,.04])
    cbar1 = fig.colorbar(cax1, cax=cbar_axes1, orientation='horizontal')
    cbar1.set_label('$\ln{\mathcal{L}}$ - median($\ln{\mathcal{L}}$)',
                    fontsize=8, labelpad=.1)

    # plot grid
    for i in range(1,N):
        ax1.axvline(i, ls='-', lw=.6, color='k')
        ax1.axhline(i, ls='-', lw=.6, color='k')

    # highlight detected period ratios
    for i in range(N):
        for j in range(N):
            col = 'k' if (lnL[i,j] < 40) | np.isnan(lnL[i,j]) else 'w'
            if i > j:
                rect = patches.Rectangle((i, j), 1, 1, color='k')
                ax1.add_patch(rect)
            else:
                entry = '%.3f'%lnL[i,j] if np.isfinite(lnL[i,j]) else 'nan' 
                ax1.text(i+.95, j+.85, entry, fontsize=5.5, color=col,
                         weight='semibold', verticalalignment='bottom',
                         horizontalalignment='right')
                
                
    ax.set_ylim((N,0))
    ax.set_xlim((0,N))
    ax.set_xlabel('index $i$', fontsize=8, labelpad=.9)
    ax.set_ylabel('index $j$', fontsize=8, labelpad=1)
    ax.set_xticks(np.arange(N)+.5)
    ax.set_xticklabels(['%i'%(i+.5) for i in np.arange(N)+.5], fontsize=7)
    ax.set_yticks(np.arange(N)+.5)
    ax.set_yticklabels(['%i'%(i+.5) for i in np.arange(N)+.5], fontsize=7)
    
    ax1.set_ylim((N,0))
    ax1.set_xlim((0,N))
    ax1.set_xlabel('index $i$', fontsize=8, labelpad=.9)
    ax1.set_ylabel('index $j$', fontsize=8, labelpad=1)
    ax1.set_xticks(np.arange(N)+.5)
    ax1.set_xticklabels(['%i'%(i+.5) for i in np.arange(N)+.5], fontsize=7)
    ax1.set_yticks(np.arange(N)+.5)
    ax1.set_yticklabels(['%i'%(i+.5) for i in np.arange(N)+.5], fontsize=7)

    fig.subplots_adjust(bottom=.24, top=.98, right=.98, left=.05, wspace=.15)
    if label:
        plt.savefig('plots/periodicsearch_%i.png'%self.tic)
    if pltt:
        plt.show()
    plt.close('all')
    


# try: d=loadpickle('PipelineResults_TIC/TIC_234994474/LC_-00099')
def plot_periodicsearch(self, P, dur_ind=0, pltt=True, label=False):
    assert self.tic == 234994474
    cols = ['k','#e31a1c','#fed976','#fc4e2a']
    t0 = 2457000
    dur_ind = int(dur_ind)
    bjdlim = self.bjd.min()-t0-.3, self.bjd.max()-t0+.3
    gapts = 1338., 1339.66
    g=(self.transit_times-t0 <= gapts[0]) | (self.transit_times-t0 >= gapts[1])

    # get periodic markers
    dur_ind = int(dur_ind)
    transit_times = self.transit_times[self.SNRs_linearsearch[:,dur_ind] >= 5] 
    g = np.append(np.where(np.diff(transit_times)>P)[0], -1)
    transit_times = transit_times[g]
    N = transit_times.size

    # get lnLOIs of TIC 234994474 POIs
    POIs, lnLOIs = np.load('input_data/POIs_lnLOIs_tic234994474.npy').T
    POIs = np.append(POIs, [16.8,19.6])
    lnLOIs = np.append(lnLOIs, np.repeat(np.nanmedian(lnLOIs),2))
    
    # create period matrix from high S/N transit times and the lnL matrix
    indarr = np.arange(1,N+2)
    Pmat = np.zeros((N,N))
    lnL = np.zeros((N,N))
    for i in range(N):
        for j in range(N):
            Pmat[i,j] = abs(transit_times[i] - transit_times[j])
            g = np.isclose(POIs, abs(Pmat[i,j]), rtol=.1)
            if np.any(g):
                lnL[i,j] = lnLOIs[g][lnLOIs[g] == lnLOIs[g].max()]
            else:
                lnL[i,j] = np.nan
    lnL -= np.nanmedian(lnL)
                
    fig = plt.figure(figsize=(6.5,3.3))
    ax = fig.add_subplot(121)
    cax = ax.pcolormesh(Pmat.T,
                        cmap=truncate_cmap(plt.get_cmap('hot_r'),.01,.85))
    cbar_axes = fig.add_axes([.05,.1,.43,.04])
    cbar = fig.colorbar(cax, cax=cbar_axes, orientation='horizontal')
    cbar.set_label('P$_{i,j}$ = T$_{0,i}$-T$_{0,j}$ [days]', fontsize=8,
                   labelpad=.1)
    
    # plot grid
    for i in range(1,N):
        ax.axvline(i, ls='-', lw=.6, color='k')
        ax.axhline(i, ls='-', lw=.6, color='k')

    # highlight detected period
    for i in range(N):
        for j in range(N):
            col = 'k' if Pmat.T[i,j] < 14 else 'w'
            if i > j:
                rect = patches.Rectangle((i, j), 1, 1, color='k')
                ax.add_patch(rect)
            else:
                ax.text(i+.95, j+.68, '%.3f\ndays'%Pmat.T[i,j], fontsize=5,
                        color=col, weight='semibold',
                        verticalalignment='bottom',
                        horizontalalignment='right')
                ax.text(i+.95, j+.89, '(%.3f)'%(Pmat.T[i,j]/P), fontsize=5,
                        color=col, weight='semibold',
                        verticalalignment='bottom',
                        horizontalalignment='right')
            
    # plot lnL matrix
    ax1 = fig.add_subplot(122)
    cax1 = ax1.pcolormesh(lnL,cmap=truncate_cmap(plt.get_cmap('hot_r'),.01,.85))
    cbar_axes1 = fig.add_axes([.55,.1,.43,.04])
    cbar1 = fig.colorbar(cax1, cax=cbar_axes1, orientation='horizontal')
    cbar1.set_label('$\ln{\mathcal{L}}$ - median($\ln{\mathcal{L}}$)',
                    fontsize=8, labelpad=.1)

    # plot grid
    for i in range(1,N):
        ax1.axvline(i, ls='-', lw=.6, color='k')
        ax1.axhline(i, ls='-', lw=.6, color='k')

    # highlight detected period ratios
    for i in range(N):
        for j in range(N):
            col = 'k' if (lnL[i,j] < 40) | np.isnan(lnL[i,j]) else 'w'
            if i > j:
                rect = patches.Rectangle((i, j), 1, 1, color='k')
                ax1.add_patch(rect)
            else:
                entry = '%.3f'%lnL[i,j] if np.isfinite(lnL[i,j]) else 'nan' 
                ax1.text(i+.95, j+.85, entry, fontsize=5.5, color=col,
                         weight='semibold', verticalalignment='bottom',
                         horizontalalignment='right')
                
                
    ax.set_ylim((N,0))
    ax.set_xlim((0,N))
    ax.set_xlabel('index $i$', fontsize=8, labelpad=.9)
    ax.set_ylabel('index $j$', fontsize=8, labelpad=1)
    ax.set_xticks(np.arange(N)+.5)
    ax.set_xticklabels(['%i'%(i+.5) for i in np.arange(N)+.5], fontsize=7)
    ax.set_yticks(np.arange(N)+.5)
    ax.set_yticklabels(['%i'%(i+.5) for i in np.arange(N)+.5], fontsize=7)
    
    ax1.set_ylim((N,0))
    ax1.set_xlim((0,N))
    ax1.set_xlabel('index $i$', fontsize=8, labelpad=.9)
    ax1.set_ylabel('index $j$', fontsize=8, labelpad=1)
    ax1.set_xticks(np.arange(N)+.5)
    ax1.set_xticklabels(['%i'%(i+.5) for i in np.arange(N)+.5], fontsize=7)
    ax1.set_yticks(np.arange(N)+.5)
    ax1.set_yticklabels(['%i'%(i+.5) for i in np.arange(N)+.5], fontsize=7)

    fig.subplots_adjust(bottom=.24, top=.98, right=.98, left=.05, wspace=.15)
    if label:
        plt.savefig('plots/periodicsearch_%i.png'%self.tic)
    if pltt:
        plt.show()
    plt.close('all')


    
# try plotting TIC 235037759
def plot_GPdetrend(self, compute_thetaGP=True, pltt=True, label=False,
                   ylim=(.75,1.32), ylimres=(-5e3,5e3), P=0, T0=0,
                   Npntsmin=5e2, Npntsmax=1e3, Nsig=3, medkernel=9,
                   sector=0, hi=2):
    #assert np.all(self.quarters==0)
    thetaGPin_final, thetaGPout_final = self.thetaGPin[0], self.thetaGPout[0]

    # black, red, yellow, orange
    cols = ['k','#a80000','#ffff5f','#','#ff4700']
    t0 = 2457000
    Nsects = np.unique(self.quarters).size
    #Nsects = 1 # TEMP
    
    # get thetaGPs
    if compute_thetaGP:
        #g = (self.bjd-self.bjd.min()<=21.88)|(self.bjd-self.bjd.min()>=24.119)
        thetain_comp2final, thetaout_comp2final, thetaGPsin, thetaGPsout = \
            _do_optimize_0_custom(self.bjd, self.f, self.ef, self.quarters)
        np.save('PipelineResults_TIC/thetain_comp2final_%i'%self.tic,
                thetain_comp2final)
        np.save('PipelineResults_TIC/thetaout_comp2final_%i'%self.tic,
                thetaout_comp2final)
        np.save('PipelineResults_TIC/thetaGPsin_%i'%self.tic, thetaGPsin)
        np.save('PipelineResults_TIC/thetaGPsout_%i'%self.tic, thetaGPsout)
    else:
        thetain_comp2final = \
            np.load('PipelineResults_TIC/thetain_comp2final_%i.npy'%self.tic)
        thetaout_comp2final = \
            np.load('PipelineResults_TIC/thetaout_comp2final_%i.npy'%self.tic)
        #thetaout_comp2finalv2 = \
        #    np.load('PipelineResults_TIC/thetaout_comp2final_%iv2.npy'%self.tic)
        thetaGPsin = np.load('PipelineResults_TIC/thetaGPsin_%i.npy'%self.tic)
        thetaGPsout = np.load('PipelineResults_TIC/thetaGPsout_%i.npy'%self.tic)
        #thetaGPsoutv2 = np.load('PipelineResults_TIC/thetaGPsout_%iv2.npy'%self.tic)

        
    #thetaGPout_finalv2 = np.load('PipelineResults_TIC/thetaout_comp2final_%iv2.npy'%self.tic)[0]

    fig = plt.figure(figsize=(7,4.2))
    ax = fig.add_axes([.1,.3,.88,.66])
    ax2 = fig.add_axes([.1,.1,.88,.2])
    dts = np.zeros(Nsects)
    for j in range(Nsects):
        # get time binning from best-fit P_GP
        Npnts_per_timescale = 8.
        timescale_to_resolve = np.exp(thetaGPsout[j,0,3]) / \
                               Npnts_per_timescale
        # bin the light curve
        s = self.quarters == j
        Ttot = self.bjd[s].max() - self.bjd[s].min()
        if Ttot/timescale_to_resolve < Npntsmin:
            dt = Ttot / Npntsmin
        elif Ttot/timescale_to_resolve > Npntsmax:
            dt = Ttot / Npntsmax
        else:
            dt = timescale_to_resolve

        # trim outliers and median filter to avoid fitting deep transits
        dts[j] = dt*24
        g = abs(self.f[s]-np.median(self.f[s])) <= Nsig*np.std(self.f[s])
        tbin, fbin, efbin = llnl.boxcar(self.bjd[s][g],
                                        medfilt(self.f[s][g],medkernel),
                                        self.ef[s][g], dt=dt)

        # plot best GP models
        ax.plot(self.bjd[s]-t0, self.f[s], 'o', c=cols[0], ms=1, alpha=.15)

        # cover up bad point
        g3 = np.in1d(np.arange(tbin.size), np.where(np.diff(tbin) > .5)[0][-1],
                     invert=True)
        ax.plot(tbin[g3]-t0, fbin[g3], 'o', ms=3, markerfacecolor=cols[2],
                markeredgecolor='k', markeredgewidth=.4)

        # plot gp model
        tmodel, mu, sig = _get_GP_model(thetaGPout_final, tbin, fbin, efbin,
                                        np.linspace(self.bjd[s].min(),
                                                    self.bjd[s].max(), 1000))
        ax.fill_between(tmodel-t0, mu-sig, mu+sig, color=cols[1], alpha=.3)
        ax.plot(tmodel-t0, mu, '-', c=cols[1])
        ax.set_xticklabels('')
        ax.set_title('TIC %i (sector %i)'%(self.tic,sector), fontsize=8,
                     weight='semibold', y=.98)
        
        # plot GP models from other iterations
        '''
        NGP = thetaGPsin.shape[1]
        assert NGP >= 9
        for i in range(9):
            tmodel, mu, sig = _get_GP_model(thetaGPsout[j,i], tbin, fbin,
                                            efbin, tmodel)
            if i == 0:
                ax.plot(tmodel-t0, mu, '-', c=cols[3], lw=.4,
                        label='remaining GP models ')
            else:
                ax.plot(tmodel-t0, mu, '-', c=cols[3], lw=.4)
        '''
        # plot residuals
        trbin, rbin, erbin = llnl.boxcar(self.bjd[s][g],
                            medfilt((self.f[s]-self.mu[s])[g]*1e6,medkernel),
                                        self.ef[s][g], dt=dt)
        ax2.plot(self.bjd[s]-t0, (self.f[s]-self.mu[s])*1e6, 'o', c=cols[0],
                 ms=1, alpha=.15) 
        ax2.plot(trbin-t0, rbin, 'o', ms=3, markerfacecolor=cols[2],
                 markeredgecolor='k', markeredgewidth=.4)
        # cover up pnt in second sector
        if j == 0:
            g3 = np.where(np.diff(trbin) > .5)[0][-1]
            ax2.plot(trbin[g3]-t0, rbin[g3], 'wo', ms=4)
        
    if Nsects > 1:
        x = self.bjd.mean()-t0+.2
        ax.axvline(x, ls='--', color='k', lw=2)
        ax2.axvline(x, ls='--', color='k', lw=2)
        ax.text(x-.7, 1.0052, 'sector 1', horizontalalignment='right',
                fontsize=6, weight='semibold')
        ax.text(x+.7, 1.0052, 'sector 2', horizontalalignment='left',
                fontsize=6, weight='semibold')
        
    ax2.set_xlabel('Time [BJD - 2,457,000]', fontsize=9)
    ax.set_ylabel('Normalized flux', fontsize=9)
    ax2.set_ylabel('De-trended Flux\n[ppm]', fontsize=9, labelpad=-3)
    ax.set_ylim(ylim)
    bjdlim = self.bjd.min()-t0-.3, self.bjd.max()-t0+.3
    ax.set_xlim(bjdlim)
    ax2.set_xlim(bjdlim)
    ax2.set_ylim(ylimres)
    print ylimres
    ax2.set_yticklabels(['','-2$\cdot$10$^{3}$','1','2$\cdot$10$^{3}$'])

    # plot transit ticks
    if P != 0 and T0 != 0:
        for i in range(-27,27):
            ax.axvline(T0-t0+P*i, ls='-', lw=1.5, c=cols[4], ymax=.05)
            ax2.axvline(T0-t0+P*i, ls='-', lw=1.5, c=cols[4], ymax=.1)
    
    # custom legend
    ax.plot(.05, .94, 'o', c=cols[0], ms=2, alpha=.8, transform=ax.transAxes)
    ax.text(.07, .94, 'raw photometry', fontsize=7, transform=ax.transAxes,
            verticalalignment='center', weight='semibold')
    ax.plot(.05, .88, 'o', markerfacecolor=cols[2], markeredgecolor='k',
            markeredgewidth=.5, ms=5, transform=ax.transAxes)
    endings = [')',')']
    dt_labels = ''.join(['dt$_%i$ = %.1f min%s'%(i+1,dts[i]*60,endings[i])
                         for i in range(Nsects)])
    ax.text(.07, .88, 'binned photometry (%s'%dt_labels, fontsize=7,
            transform=ax.transAxes, verticalalignment='center',
            weight='semibold')
    
    ax.fill_between([.63,.67], [.92,.92], [.94,.94], color=cols[1], alpha=.3,
                    transform=ax.transAxes)
    ax.plot([.63,.67], [.93,.93], '-', color=cols[1], lw=2,
            transform=ax.transAxes)
    ax.text(.69, .93,
            'max($\ln{\mathcal{L}}$) GP models mean and\nstandard deviation',
            fontsize=7, transform=ax.transAxes, verticalalignment='center',
            weight='semibold')
    if P != 0 and T0 != 0:
        ax.plot([.655,.655], [.84,.88], '-', color=cols[4], lw=1.5,
                transform=ax.transAxes)
        ax.text(.69,.86, 'transit markers (P = %.3f days)'%P, fontsize=7,
                transform=ax.transAxes, verticalalignment='center',
                weight='semibold')
    
    fig.subplots_adjust(bottom=.1, top=.96, right=.98, left=.05)
    if label:
        plt.savefig('plots/GPdetrend_%i.png'%self.tic)
    if pltt:
        plt.show()
    plt.close('all')



def _do_optimize_0_custom(bjd, f, ef, quarters, N=10,
                          Npntsmin=5e2, Npntsmax=1e3, Nsig=3, medkernel=9):
    '''First fit the PDC LC with GP and a no planet model using an optimization
    routine and tested with N different initializations.'''
    # fit GP separately to each quarter (or K2 campaign)
    assert bjd.size == f.size
    assert bjd.size == ef.size
    assert bjd.size == quarters.size
    NGP = np.unique(quarters).size

    # test various hyperparameter initializations and keep the one resulting
    # in the most gaussian-like residuals
    N = int(N)
    thetaGPs_in_tmp, thetaGPs_out_tmp = np.zeros((NGP,N,4)), \
                                        np.zeros((NGP,N,4))
    thetaGPs_in, thetaGPs_out = np.zeros((NGP,4)), np.zeros((NGP,4))
    for i in range(NGP):
        g1 = quarters == i

        #pvalues = np.zeros(N)
        lnLs = np.zeros(N)
        for j in range(N):
            # get initial GP parameters (P, and l from periodogram)
            thetaGPs_in_tmp[i,j] = initialize_GP_hyperparameters(bjd[g1], f[g1],
                                                                 ef[g1],
                                                                 Pindex=j)
            Npnts_per_timescale = 8.
            inds = np.array([1,3])
            timescale_to_resolve = np.exp(thetaGPs_in_tmp[i,j,inds]).min() / \
                                   Npnts_per_timescale
            # bin the light curve
            Ttot = bjd[g1].max() - bjd[g1].min()
            if Ttot/timescale_to_resolve < Npntsmin:
                dt = Ttot / Npntsmin
            elif Ttot/timescale_to_resolve > Npntsmax:
                dt = Ttot / Npntsmax
            else:
                dt = timescale_to_resolve

            # trim outliers and median filter to avoid fitting deep transits
            g = abs(f[g1]-np.median(f[g1])) <= Nsig*np.std(f[g1])
            tbin, fbin, efbin = llnl.boxcar(bjd[g1][g], medfilt(f[g1][g],medkernel),
                                            ef[g1][g], dt=dt)
            gp,mu,sig,thetaGPs_out_tmp[i,j] = fit_GP_0(thetaGPs_in_tmp[i,j],
                                                       tbin, fbin, efbin)
            # compute residuals and the normality test p-value
            #_,pvalues[j] = normaltest(fbin-mu)
            lnLs[j] = llnl.lnlike(tbin, fbin, efbin, mu)
            
        # select the most gaussian-like residuals
        g = np.argsort(lnLs)[-1]
        thetaGPs_in[i]  = thetaGPs_in_tmp[i,g]
        thetaGPs_out[i] = thetaGPs_out_tmp[i,g]

    return thetaGPs_in, thetaGPs_out, thetaGPs_in_tmp, thetaGPs_out_tmp


def _get_GP_model(thetaGP, tbin, fbin, efbin, tmodel):
    assert len(thetaGP) == 4
    a, l, G, Pgp = np.exp(thetaGP)
    k1 = george.kernels.ExpSquaredKernel(l)
    k2 = george.kernels.ExpSine2Kernel(G,Pgp)
    gp = george.GP(a*(k1+k2))
    try:
        gp.compute(tbin, efbin)
    except (ValueError, np.linalg.LinAlgError):
        return np.repeat(None,4)
    mu, cov = gp.predict(fbin, tmodel)
    sig = np.sqrt(np.diag(cov))
    return tmodel, mu, sig    







def _do_optimize_0_custom(bjd, f, ef, quarters, N=10,
                          Npntsmin=5e2, Npntsmax=1e3, Nsig=3, medkernel=9):
    '''First fit the PDC LC with GP and a no planet model using an optimization
    routine and tested with N different initializations.'''
    # fit GP separately to each quarter (or K2 campaign)
    assert bjd.size == f.size
    assert bjd.size == ef.size
    assert bjd.size == quarters.size
    NGP = np.unique(quarters).size

    # test various hyperparameter initializations and keep the one resulting
    # in the most gaussian-like residuals
    N = int(N)
    thetaGPs_in_tmp, thetaGPs_out_tmp = np.zeros((NGP,N,4)), \
                                        np.zeros((NGP,N,4))
    thetaGPs_in, thetaGPs_out = np.zeros((NGP,4)), np.zeros((NGP,4))
    for i in range(NGP):
        g1 = quarters == i

        #pvalues = np.zeros(N)
        lnLs = np.zeros(N)
        for j in range(N):
            # get initial GP parameters (P, and l from periodogram)
            thetaGPs_in_tmp[i,j] = initialize_GP_hyperparameters(bjd[g1], f[g1],
                                                                 ef[g1],
                                                                 Pindex=j)
            Npnts_per_timescale = 8.
            inds = np.array([1,3])
            timescale_to_resolve = np.exp(thetaGPs_in_tmp[i,j,inds]).min() / \
                                   Npnts_per_timescale
            # bin the light curve
            Ttot = bjd[g1].max() - bjd[g1].min()
            if Ttot/timescale_to_resolve < Npntsmin:
                dt = Ttot / Npntsmin
            elif Ttot/timescale_to_resolve > Npntsmax:
                dt = Ttot / Npntsmax
            else:
                dt = timescale_to_resolve

            # trim outliers and median filter to avoid fitting deep transits
            g = abs(f[g1]-np.median(f[g1])) <= Nsig*np.std(f[g1])
            tbin, fbin, efbin = llnl.boxcar(bjd[g1][g], medfilt(f[g1][g],medkernel),
                                            ef[g1][g], dt=dt)
            gp,mu,sig,thetaGPs_out_tmp[i,j] = fit_GP_0(thetaGPs_in_tmp[i,j],
                                                       tbin, fbin, efbin)
            # compute residuals and the normality test p-value
            #_,pvalues[j] = normaltest(fbin-mu)
            lnLs[j] = llnl.lnlike(tbin, fbin, efbin, mu)
            
        # select the most gaussian-like residuals
        g = np.argsort(lnLs)[-1]
        thetaGPs_in[i]  = thetaGPs_in_tmp[i,g]
        thetaGPs_out[i] = thetaGPs_out_tmp[i,g]

    return thetaGPs_in, thetaGPs_out, thetaGPs_in_tmp, thetaGPs_out_tmp


def _get_GP_model(thetaGP, tbin, fbin, efbin, tmodel):
    assert len(thetaGP) == 4
    a, l, G, Pgp = np.exp(thetaGP)
    k1 = george.kernels.ExpSquaredKernel(l)
    k2 = george.kernels.ExpSine2Kernel(G,Pgp)
    gp = george.GP(a*(k1+k2))
    try:
        gp.compute(tbin, efbin)
    except (ValueError, np.linalg.LinAlgError):
        return np.repeat(None,4)
    mu, cov = gp.predict(fbin, tmodel)
    sig = np.sqrt(np.diag(cov))
    return tmodel, mu, sig


def plot_stellar_corner(dTIC, dTIC_sect1, dTIC_sect2, pltt=True, label=False):
    assert dTIC.shape == (93090,89)  #np.genfromtxt('input_data/TESStargets/TICv7_Mdwarfsv1.csv', delimiter=',', skip_header=5)
    assert dTIC_sect1.shape == (905, 39) #np.loadtxt('input_data/TESStargets/TESSMdwarfs_sector1_v2.csv', delimiter=',')
    assert dTIC_sect2.shape == (1062,39) #np.loadtxt('input_data/TESStargets/TESSMdwarfs_sector2_v2.csv', delimiter=',')

    # get only TICs in sects 1 or 2
    tic1 = np.loadtxt('input_data/TESStargets/all_targets_S001_v1.csv',
                      delimiter=',', skiprows=6)[:,0]
    tic2 = np.loadtxt('input_data/TESStargets/all_targets_S002_v1.csv',
                      delimiter=',', skiprows=6)[:,0]
    g = np.in1d(dTIC[:,0], tic1) | np.in1d(dTIC[:,0], tic2)
    dTIC = dTIC[g]
    
    # black, red, yellow, orange
    cols = ['#ff4700','#a80000','#ffff5f']
    #cols = ['#fc4e2a','#800026','#e31a1c','#bd0026','#fed976']
    
    # setup arrays for pre and post GAIA stellar parameters
    inds1 = np.array([60,64,70,72])
    TESSmag1, Teff1, Rs1, Ms1 = dTIC[:,inds1].T
    list1 = [TESSmag1, Teff1, Rs1, Ms1]

    inds2 = np.array([0,7,27,30,33])
    TICs21, TESSmag21, Rs21, Teff21, Ms21 = dTIC_sect1[:,inds2].T
    TICs22, TESSmag22, Rs22, Teff22, Ms22 = dTIC_sect2[:,inds2].T
    _,inds = np.unique(np.append(TICs21,TICs22), return_index=True)
    TESSmag2 = np.append(TESSmag21,TESSmag22)[inds]
    Teff2 = np.append(Teff21,Teff22)[inds]
    Rs2 = np.append(Rs21,Rs22)[inds]
    Ms2 = np.append(Ms21,Ms22)[inds]
    list2 = [TESSmag2, Teff2, Rs2, Ms2]

    labels = ['$T$','$T_{eff}$ [K]','$R_s$ [R$_{\odot}$]','$M_s$ [M$_{\odot}$]']
    lims = [(5,16),(2700,4100),(.1,.65),(.1,.65)]
    titles1 = ['%.1f,','%i,','%.2f,','%.2f,']
    titles2 = ['%.1f','%i K','%.2f R$_{\odot}$','%.2f M$_{\odot}$']
    tickvals = [range(6,16,2),np.arange(28,41,4)*1e2,np.arange(.2,.7,.2),
                np.arange(.2,.7,.2)]
    ticklabels = [['6','8','10','12','14'], ['2800','3200','3600','4000'],
                  ['0.2','0.4','0.6'], ['0.2','0.4','0.6']]
    bins = [np.linspace(7,16,20), np.linspace(2700,4050,20),
            np.linspace(.1,.65,20), np.linspace(.1,.65,20)]
    
    fig = plt.figure(figsize=(3.7,4))
    subplot = 1
    for i in range(4):
        for j in range(4):

            ax = fig.add_subplot(4,4,subplot)
            subplot += 1
            
            # plot 1d histogram
            if i == j:
                y1,_,_=ax.hist(list1[i], bins=bins[i], color=cols[0], normed=1,
                              histtype='step',lw=1)
                ax.hist(list1[i], bins=bins[i], color=cols[0], normed=1,
                        alpha=.5)
                y2,_,_=ax.hist(list2[i], bins=bins[i], color=cols[1], normed=1,
                               histtype='step',lw=1.5)
                ax.hist(list2[i], bins=bins[i], color=cols[1], normed=1,
                        alpha=.6)
                ax.set_xlim(lims[i])
                if j < 3: ax.set_xticklabels('')
                if j == 3:
                    ax.set_xlabel(labels[j], fontsize=7, labelpad=1)
                    ax.set_xticks(tickvals[j])
                    ax.set_xticklabels(ticklabels[j], fontsize=5)
                ax.set_yticklabels('')
                ax.plot(np.repeat(np.nanmedian(list1[i]),2),[0,np.max([y1,y2])],
                        '--', lw=.8,c=cols[0])
                ax.plot(np.repeat(np.nanmedian(list2[i]),2),[0,np.max([y1,y2])],
                        '--', lw=.8, c=cols[1])
                ax.text(.49, 1.05, titles1[i]%np.nanmedian(list1[i]),
                        fontsize=7,
                        transform=ax.transAxes, horizontalalignment='right',
                        weight='semibold', color=cols[0])
                ax.text(.51, 1.05, titles2[i]%np.nanmedian(list2[i]),
                        fontsize=7,
                        transform=ax.transAxes, horizontalalignment='left',
                        weight='semibold', color=cols[1])
                
            # 2d histogram
            elif i > j:
                print len(list1[j]), len(list2[j])
                ax.plot(list1[j], list1[i], '.', c=cols[0], ms=1, alpha=.8)
                ax.plot(list2[j], list2[i], '.', c=cols[1], ms=1.5, alpha=.5)
                ax.set_xlim(lims[j]), ax.set_ylim(lims[i])
                if i < 3: ax.set_xticklabels('')
                if j > 0: ax.set_yticklabels('')
                if j == 0:
                    ax.set_ylabel(labels[i], fontsize=7, labelpad=1.5)
                    ax.set_yticks(tickvals[i])
                    ax.set_yticklabels(ticklabels[i], fontsize=5)
                if i == 3:
                    ax.set_xlabel(labels[j], fontsize=7, labelpad=1.5)
                    ax.set_xticks(tickvals[j])
                    ax.set_xticklabels(ticklabels[j], fontsize=5, rotation=45)
                    
            # no plot
            else:
                ax.axis('off')

    fig.subplots_adjust(top=.96, right=.98, left=.1, bottom=.09, hspace=.05, wspace=.05)
    if label:
        plt.savefig('plots/stellar_corner.png')
    if pltt:
        plt.show()
    plt.close('all')

    


def plot_planet_population(self, pltt=True, label=False):
    assert hasattr(self, 'isTESSalert')
    ta = (self.isTESSalert==1) & np.in1d(self.disposition_human, range(2)) 
    g0 = self.disposition_human == 0
    g1 = self.disposition_human == 1
    g2 = self.disposition_human == 2
    g25 = self.disposition_human == 2.5
    
    # light orange, red, yellow, orange
    cols = ['#fe9500','#a80000','#ffff5f','#ff4700','#fe9500']
    #cols = ['b','k','g','r','c']
    marker = ['o','^','p','d']
    rplim = (.7,4.5)

    # undetected TOIs
    tois_UD = [203.01,221.01,175.02,175.03]
    assert np.any(np.in1d(tois_UD, self.tois, invert=True)) 
    tics = [259962054,316937670,307210830,307210830]
    gt = np.in1d(self.tics, tics)
    Ls_UD = compute_Ls(self.Rss[gt], self.Teffs[gt])
    Ps_UD = [52,.624, 7.45,2.25]
    Fs_UD = compute_F(Ls_UD, rvs.semimajoraxis(np.array(Ps_UD),self.Mss[gt],0))
    rps_UD = [1.22,1.67, 1.43,.794]
    
    # plot planet candidates that passed human vetting (and maybe vespa?)
    fig = plt.figure(figsize=(6.5,2.9))

    # plot rp vs P
    ax1 = fig.add_subplot(121)
    ax1.errorbar(self.Ps[g0], self.rps[g0], xerr=self.e_Ps[g0],
                 yerr=self.ehi_rps[g0], fmt=marker[1], ms=5,
                 markerfacecolor=cols[1],
                 markeredgecolor='k', markeredgewidth=.5, ecolor='k',
                 elinewidth=.9, capsize=0, label='pPC')
    ax1.errorbar(self.Ps[g1], self.rps[g1], xerr=self.e_Ps[g1],
                 yerr=self.ehi_rps[g1], fmt=marker[0], ms=6,
                 markerfacecolor=cols[1],
                 markeredgecolor='k', markeredgewidth=.5, ecolor='k',
                 elinewidth=.9, capsize=0, label='PC')
    
    ax1.errorbar(self.Ps_singletransit[g2], self.rps_singletransit[g2],
                 xerr=[self.elo_Ps_singletransit[g2],
                       self.ehi_Ps_singletransit[g2]],
		 yerr=[self.elo_rps_singletransit[g2],
                       self.ehi_rps_singletransit[g2]], fmt=marker[2], ms=5,
                 markerfacecolor=cols[2],
                 markeredgecolor='k', markeredgewidth=.5, ecolor='k',
                 elinewidth=.9, capsize=0, label='single transit')
    ax1.errorbar(self.Ps_singletransit[g25], self.rps_singletransit[g25],
                 xerr=[self.elo_Ps_singletransit[g25],
                       self.ehi_Ps_singletransit[g25]],
		 yerr=[self.elo_rps_singletransit[g25],
                       self.ehi_rps_singletransit[g25]], fmt=marker[1], ms=5,
                 markerfacecolor=cols[2],
                 markeredgecolor='k', markeredgewidth=.5, ecolor='k',
                 elinewidth=.9, capsize=0)

    
    ax1.plot(self.Ps[ta], self.rps[ta], marker[3], markerfacecolor=cols[3],
             markeredgecolor='k', markeredgewidth=.3, ms=10,
             label='TESS alert')
    ax1.plot(Ps_UD, rps_UD, marker[3], markerfacecolor=cols[4],
             markeredgecolor='k', markeredgewidth=.3, ms=4)

    # TOI labels
    dx = [2,-2.3,-.7,-1.5,-5,-3,-1.85,-2.5]
    dy = [.2,.17,-.15,-.1,-.23,.2,-.1,.2]
    #for i in range(ta.sum()):
        #ax1.text(self.Ps[ta][i]+dx[i], self.rps[ta][i]+dy[i],
        #         '%.2f'%self.tois[ta][i], fontsize=6, weight='semibold')

    dx = [5,-1.4,-.07,.6,.25,.1,.1,.1]
    dy = [.02,-.17,.1,-.1,.01,.1,.1,.1]
    #for i in range(len(Ps_UD)):
    #    print Ps_UD[i], tois_UD[i]
        #ax1.text(Ps_UD[i]+dx[i], rps_UD[i]+dy[i], '%.2f'%tois_UD[i],
        #         fontsize=5, weight='semibold')
        
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xlabel('Orbital Period [days]', fontsize=10)
    ax1.set_ylabel('Planet Radius [R$_{\oplus}$]', fontsize=10)
    ax1.set_xlim((.5,2e2))
    ax1.set_ylim(rplim)
    ax1.yaxis.set_major_formatter(NullFormatter())
    ax1.yaxis.set_minor_formatter(NullFormatter())
    ax1.set_xticks([1,10,100])
    ax1.set_xticklabels(['1','10','100'])
    ax1.set_yticks(range(1,5))
    ax1.set_yticklabels(['%i'%i for i in range(1,5)])
    
    # plot rp vs F
    ax2 = fig.add_subplot(122)
    ax2.errorbar(self.Fs[g0], self.rps[g0], xerr=[self.elo_Fs[g0],
                                                  self.ehi_Fs[g0]],
                 yerr=self.ehi_rps[g0], fmt=marker[1], ms=5,
                 markerfacecolor=cols[1],
                 markeredgecolor='k', markeredgewidth=.5, ecolor='k',
                 elinewidth=.9, capsize=0, label='pFC')
    ax2.errorbar(self.Fs[g1], self.rps[g1], xerr=[self.elo_Fs[g1],
                                                  self.ehi_Fs[g1]],
                 yerr=self.ehi_rps[g1], fmt=marker[0], ms=6,
                 markerfacecolor=cols[1],
                 markeredgecolor='k', markeredgewidth=.5, ecolor='k',
                 elinewidth=.9, capsize=0, label='FC')
    
    ax2.errorbar(self.Fs_singletransit[g2], self.rps_singletransit[g2],
                 xerr=[self.elo_Fs_singletransit[g2],
                       self.ehi_Fs_singletransit[g2]],
		 yerr=[self.elo_rps_singletransit[g2],
                       self.ehi_rps_singletransit[g2]], fmt=marker[2], ms=5,
                 markerfacecolor=cols[2],
                 markeredgecolor='k', markeredgewidth=.5, ecolor='k',
                 elinewidth=.9, capsize=0, label='single transit')
    ax2.errorbar(self.Fs_singletransit[g25], self.rps_singletransit[g25],
                 xerr=[self.elo_Fs_singletransit[g25],
                       self.ehi_Fs_singletransit[g25]],
		 yerr=[self.elo_rps_singletransit[g25],
                       self.ehi_rps_singletransit[g25]], fmt=marker[1], ms=5,
                 markerfacecolor=cols[2],
                 markeredgecolor='k', markeredgewidth=.5, ecolor='k',
                 elinewidth=.9, capsize=0)
    
    ax2.plot(self.Fs[ta], self.rps[ta], marker[3], markerfacecolor=cols[3],
             markeredgecolor='k', markeredgewidth=.3, ms=10,
             label='TESS alert')
    ax2.plot(Fs_UD, rps_UD, marker[3], markerfacecolor=cols[4],
             markeredgecolor='k', markeredgewidth=.3, ms=4)

    # plot the HZ
    ax2.fill_between([.2,1.5], np.repeat(np.min(rplim),2),
                     np.repeat(np.max(rplim),2), color='k', alpha=.15)
    ax2.fill_between([.22,.9], np.repeat(np.min(rplim),2),
                     np.repeat(np.max(rplim),2), color='k', alpha=.15)
    ax2.axvline(1.5, ls='--', lw=.6, color='k')
    ax2.axvline(.2, ls='--', lw=.6, color='k')
    ax2.axvline(.22, ls='-.', lw=.6, color='k')
    ax2.axvline(.9, ls='-.', lw=.6, color='k')
    
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    ax2.set_xlabel('Insolation [S$_{\oplus}$]', fontsize=10)
    ax2.set_xlim((2e2,.1))
    ax2.set_ylim(rplim)
    ax2.yaxis.set_major_formatter(NullFormatter())
    ax2.yaxis.set_minor_formatter(NullFormatter())
    ax2.set_xticks([.1,1,10,100])
    ax2.set_xticklabels(['0.1','1','10','100'])
    ax2.set_yticks(range(1,5))
    ax2.set_yticklabels('')

    # custom legend
    ax2.plot(.07, .93, marker[0], markerfacecolor=cols[1], ms=6,
             transform=ax2.transAxes, markeredgecolor='k', markeredgewidth=.5)
    ax2.text(.11, .93, 'PC', fontsize=6, transform=ax2.transAxes,
             verticalalignment='center', weight='semibold')
    ax2.plot(.07, .88, marker[1], markerfacecolor=cols[1], ms=5,
             transform=ax2.transAxes, markeredgecolor='k', markeredgewidth=.5)
    ax2.text(.11, .88, 'pPC', fontsize=6, transform=ax2.transAxes,
             verticalalignment='center', weight='semibold')
    ax2.plot(.07, .83, marker[2], markerfacecolor=cols[2], ms=5,
             transform=ax2.transAxes, markeredgecolor='k', markeredgewidth=.5)
    ax2.text(.11, .83, 'ST', fontsize=6, transform=ax2.transAxes,
             verticalalignment='center', weight='semibold')
    ax2.plot(.07, .78, marker[1], markerfacecolor=cols[2], ms=5,
             transform=ax2.transAxes, markeredgecolor='k', markeredgewidth=.5)
    ax2.text(.11, .78, 'pST', fontsize=6, transform=ax2.transAxes,
             verticalalignment='center', weight='semibold')
    ax2.plot(.07, .73, marker[3], markerfacecolor=cols[3], ms=6,
             transform=ax2.transAxes, markeredgecolor='k', markeredgewidth=.3)
    ax2.text(.11, .73, 'detected TOI', fontsize=6, transform=ax2.transAxes,
             verticalalignment='center', weight='semibold')
    ax2.plot(.07, .68, marker[3], markerfacecolor=cols[4], ms=4,
             transform=ax2.transAxes, markeredgecolor='k', markeredgewidth=.3)
    ax2.text(.11, .68, 'undetected TOI', fontsize=6, transform=ax2.transAxes,
             verticalalignment='center', weight='semibold')
    
    fig.subplots_adjust(bottom=.15, top=.97, right=.98, left=.06, wspace=.05)
    if label:
        plt.savefig('plots/planetsample.png')
    if pltt:
        plt.show()
    plt.close('all')



def plot_planet_population_star(self, pltt=True, label=False):
    assert hasattr(self, 'isTESSalert')
    ta = (self.isTESSalert==1) & np.in1d(self.disposition_human, range(2)) 
    g0 = self.disposition_human == 0
    g1 = self.disposition_human == 1
    g2 = self.disposition_human == 2
    g25 = self.disposition_human == 2.5

    
    # light orange, red, yellow, orange
    cols = ['#fe9500','#a80000','#ffff5f','#ff4700','#fe9500']
    #cols = ['b','k','g','r','c']
    marker = ['o','^','p','d']
    rplim = (.7,5)

    # undetected TOIs
    tois_UD = [203.01,221.01,175.02,175.03]
    assert np.any(np.in1d(tois_UD, self.tois, invert=True)) 
    tics = [259962054,316937670,307210830,307210830]
    udta = np.where(np.in1d(self.tics, tics))[0]
    #FPtoi = np.where(np.array(tics) == 307210830)[0]
    #if udta.size != len(tics):
    #    udta = np.append(udta, udta[-1])
    gt = np.in1d(self.tics, tics)
    rps_UD = np.array([1.22,1.67, 1.43,.794])

    # get known transiting planets (NASA archive Dec 15)
    fname = 'input_data/TESStargets/NASAarchive_transitingMdwarfplanets.csv'
    distT,TeffT,rpT,JmagT = np.loadtxt(fname, delimiter=',',
                                       usecols=(12,16,33,38), skiprows=57).T
    
    # plot planet candidates that passed human vetting (and maybe vespa?)
    fig = plt.figure(figsize=(6.5,2.9))

    # plot rp vs Teff
    ax1 = fig.add_subplot(131)
    ax1.errorbar(self.Teffs[g0], self.rps[g0], xerr=self.ehi_Teffs[g0],
                 yerr=self.ehi_rps[g0], fmt=marker[1], ms=5,
                 markerfacecolor=cols[1],
                 markeredgecolor='k', markeredgewidth=.5, ecolor='k',
                 elinewidth=.9, capsize=0, label='pFC')
    ax1.errorbar(self.Teffs[g1], self.rps[g1], xerr=self.elo_Teffs[g1],
                 yerr=self.ehi_rps[g1], fmt=marker[0], ms=6,
                 markerfacecolor=cols[1],
                 markeredgecolor='k', markeredgewidth=.5, ecolor='k',
                 elinewidth=.9, capsize=0, label='FC')
    ax1.errorbar(self.Teffs[g2], self.rps_singletransit[g2],
                 xerr=self.ehi_Teffs[g2],
		 yerr=[self.elo_rps_singletransit[g2],
                       self.ehi_rps_singletransit[g2]], fmt=marker[2], ms=5,
                 markerfacecolor=cols[2],
                 markeredgecolor='k', markeredgewidth=.5, ecolor='k',
                 elinewidth=.9, capsize=0, label='single transit')
    ax1.errorbar(self.Teffs[g25], self.rps_singletransit[g25],
                 xerr=self.elo_Teffs[g25],
		 yerr=[self.elo_rps_singletransit[g25],
                       self.ehi_rps_singletransit[g25]], fmt=marker[1], ms=5,
                 markerfacecolor=cols[2],
                 markeredgecolor='k', markeredgewidth=.5, ecolor='k',
                 elinewidth=.9, capsize=0)

    ax1.plot(self.Teffs[ta], self.rps[ta], marker[3], markerfacecolor=cols[3],
             markeredgecolor='k', markeredgewidth=.3, ms=10,
             label='TESS alert')
    ax1.plot(self.Teffs[udta], rps_UD, marker[3], markerfacecolor=cols[4],
             markeredgecolor='k', markeredgewidth=.3, ms=4)
    ax1.plot(self.Teffs[udta], rps_UD, marker[3], markerfacecolor=cols[4],
             markeredgecolor='k', markeredgewidth=.3, ms=4)
    #ax1.plot(self.Teffs[udta][FPtoi], rps_UD[FPtoi], marker[3],
    #         markerfacecolor='k', markeredgecolor='k', markeredgewidth=.3,
    #         ms=4)
    ax1.plot(TeffT, rpT, 'o', ms=1, c='k', alpha=.5)

    # TOI labels
    dx = [2,-2.3,-.7,-1.5,-5,-3,-1.85,-2.5]
    dy = [.2,.17,-.15,-.1,-.23,.2,-.1,.2]
    #for i in range(ta.sum()):
    #    ax1.text(self.Ps[ta][i]+dx[i], self.rps[ta][i]+dy[i],
    #             '%.2f'%self.tois[ta][i], fontsize=6, weight='semibold')

    dx = [5,-1.4,-.07,.6,.25,.1,.1,.1]
    dy = [.02,-.17,.1,-.1,.01,.1,.1,.1]
    
    ax1.set_yscale('log')
    ax1.set_xlabel('Stellar Effective Temperature [K]', fontsize=10)
    ax1.set_ylabel('Planet Radius [R$_{\oplus}$]', fontsize=10)
    #ax1.set_xlim((.5,2e2))
    ax1.set_ylim(rplim)
    ax1.yaxis.set_major_formatter(NullFormatter())
    ax1.yaxis.set_minor_formatter(NullFormatter())
    #ax1.set_xticks([1,10,100])
    #ax1.set_xticklabels(['1','10','100'])
    ax1.set_yticks(range(1,6))
    ax1.set_yticklabels(['%i'%i for i in range(1,6)])
    
    # plot rp vs Jmag
    ax2 = fig.add_subplot(132)
    ax2.errorbar(self.Jmags[g0], self.rps[g0], xerr=self.e_Jmags[g0],
                 yerr=self.ehi_rps[g0], fmt=marker[1], ms=5,
                 markerfacecolor=cols[1],
                 markeredgecolor='k', markeredgewidth=.5, ecolor='k',
                 elinewidth=.9, capsize=0, label='pFC')
    ax2.errorbar(self.Jmags[g1], self.rps[g1], xerr=self.e_Jmags[g1],
                 yerr=self.ehi_rps[g1], fmt=marker[0], ms=6,
                 markerfacecolor=cols[1],
                 markeredgecolor='k', markeredgewidth=.5, ecolor='k',
                 elinewidth=.9, capsize=0, label='FC')
    ax2.errorbar(self.Jmags[g2], self.rps_singletransit[g2],
                 xerr=self.e_Jmags[g2],
		 yerr=[self.elo_rps_singletransit[g2],
                       self.ehi_rps_singletransit[g2]], fmt=marker[2], ms=5,
                 markerfacecolor=cols[2],
                 markeredgecolor='k', markeredgewidth=.5, ecolor='k',
                 elinewidth=.9, capsize=0, label='single transit')
    ax2.errorbar(self.Jmags[g25], self.rps_singletransit[g25],
                 xerr=self.e_Jmags[g25],
		 yerr=[self.elo_rps_singletransit[g25],
                       self.ehi_rps_singletransit[g25]], fmt=marker[1], ms=5,
                 markerfacecolor=cols[2],
                 markeredgecolor='k', markeredgewidth=.5, ecolor='k',
                 elinewidth=.9, capsize=0)

    ax2.plot(self.Jmags[ta], self.rps[ta], marker[3], markerfacecolor=cols[3],
             markeredgecolor='k', markeredgewidth=.3, ms=10,
             label='TESS alert')
    ax2.plot(self.Jmags[udta], rps_UD, marker[3], markerfacecolor=cols[4],
             markeredgecolor='k', markeredgewidth=.3, ms=4)
    ##ax2.plot(self.Jmags[udta][FPtoi], rps_UD[FPtoi], marker[3],
    ##         markerfacecolor='k',
    ##         markeredgecolor='k', markeredgewidth=.3, ms=4)
    ax2.plot(JmagT, rpT, 'o', ms=1, c='k', alpha=.5)
    
    ax2.set_yscale('log')
    ax2.set_xlabel('$J$', fontsize=10)
    #ax2.set_xlim((2e2,.1))
    ax2.set_ylim(rplim)
    ax2.yaxis.set_major_formatter(NullFormatter())
    ax2.yaxis.set_minor_formatter(NullFormatter())
    #ax2.set_xticks([.1,1,10,100])
    #ax2.set_xticklabels(['0.1','1','10','100'])
    ax2.set_yticks(range(1,6))
    ax2.set_yticklabels('')

    # custom legend
    ax2.plot(.04, .95, marker[0], markerfacecolor='k', ms=2, alpha=.7,
             transform=ax2.transAxes)
    ax2.text(.08, .95, 'confirmed transiting planets', fontsize=5,
             transform=ax2.transAxes,verticalalignment='center',
             weight='semibold')
    ax2.plot(.04, .9, marker[0], markerfacecolor=cols[1], ms=6,
             transform=ax2.transAxes, markeredgecolor='k', markeredgewidth=.5)
    ax2.text(.08, .9, 'PC', fontsize=5, transform=ax2.transAxes,
             verticalalignment='center', weight='semibold')
    ax2.plot(.04, .85, marker[1], markerfacecolor=cols[1], ms=5,
             transform=ax2.transAxes, markeredgecolor='k', markeredgewidth=.5)
    ax2.text(.08, .85, 'pPC', fontsize=5, transform=ax2.transAxes,
             verticalalignment='center', weight='semibold')
    ax2.plot(.04, .8, marker[2], markerfacecolor=cols[2], ms=5,
             transform=ax2.transAxes, markeredgecolor='k', markeredgewidth=.5)
    ax2.text(.08, .8, 'ST', fontsize=5, transform=ax2.transAxes,
             verticalalignment='center', weight='semibold')
    ax2.plot(.04, .75, marker[1], markerfacecolor=cols[2], ms=5,
             transform=ax2.transAxes, markeredgecolor='k', markeredgewidth=.5)
    ax2.text(.08, .75, 'pST', fontsize=5, transform=ax2.transAxes,
             verticalalignment='center', weight='semibold')
    ax2.plot(.04, .7, marker[3], markerfacecolor=cols[3], ms=6,
             transform=ax2.transAxes, markeredgecolor='k', markeredgewidth=.3)
    ax2.text(.08, .7, 'detected TOI', fontsize=5, transform=ax2.transAxes,
             verticalalignment='center', weight='semibold')
    ax2.plot(.04, .65, marker[3], markerfacecolor=cols[4], ms=4,
             transform=ax2.transAxes, markeredgecolor='k', markeredgewidth=.3)
    ax2.text(.08, .65, 'undetected TOI', fontsize=5, transform=ax2.transAxes,
             verticalalignment='center', weight='semibold')
    #ax2.plot(.04, .6, marker[3], markerfacecolor='k', ms=4,
    #         transform=ax2.transAxes, markeredgecolor='k', markeredgewidth=.3)
    #ax2.text(.08, .6, 'AFP TOI', fontsize=5, transform=ax2.transAxes,
    #         verticalalignment='center', weight='semibold')

    # plot rp vs dist
    ax3 = fig.add_subplot(133)
    ax3.errorbar(self.dists[g0], self.rps[g0], xerr=self.ehi_dists[g0],
                 yerr=self.ehi_rps[g0], fmt=marker[1], ms=5,
                 markerfacecolor=cols[1],
                 markeredgecolor='k', markeredgewidth=.5, ecolor='k',
                 elinewidth=.9, capsize=0, label='pFC')
    ax3.errorbar(self.dists[g1], self.rps[g1], xerr=self.ehi_dists[g1],
                 yerr=self.ehi_rps[g1], fmt=marker[0], ms=6,
                 markerfacecolor=cols[1],
                 markeredgecolor='k', markeredgewidth=.5, ecolor='k',
                 elinewidth=.9, capsize=0, label='FC')
    ax3.errorbar(self.dists[g2], self.rps_singletransit[g2],
                 xerr=self.ehi_dists[g2],
		 yerr=[self.elo_rps_singletransit[g2],
                       self.ehi_rps_singletransit[g2]], fmt=marker[2], ms=5,
                 markerfacecolor=cols[2],
                 markeredgecolor='k', markeredgewidth=.5, ecolor='k',
                 elinewidth=.9, capsize=0, label='single transit')
    ax3.errorbar(self.dists[g25], self.rps_singletransit[g25],
                 xerr=self.ehi_dists[g25],
		 yerr=[self.elo_rps_singletransit[g25],
                       self.ehi_rps_singletransit[g25]], fmt=marker[1], ms=5,
                 markerfacecolor=cols[2],
                 markeredgecolor='k', markeredgewidth=.5, ecolor='k',
                 elinewidth=.9, capsize=0)

    ax3.plot(self.dists[ta], self.rps[ta], marker[3], markerfacecolor=cols[3],
             markeredgecolor='k', markeredgewidth=.3, ms=10,
             label='TESS alert')
    ax3.plot(self.dists[udta], rps_UD, marker[3], markerfacecolor=cols[4],
             markeredgecolor='k', markeredgewidth=.3, ms=4)
    #ax3.plot(self.dists[udta][FPtoi], rps_UD[FPtoi], marker[3],
    #         markerfacecolor='k',
    #         markeredgecolor='k', markeredgewidth=.3, ms=4)
    ax3.plot(distT, rpT, 'o', ms=1, c='k', alpha=.5)
    
    ax3.set_xscale('log')
    ax3.set_yscale('log')
    ax3.set_xlabel('Distance [pc]', fontsize=10)
    ax3.set_xlim((8,6e2))
    ax3.set_ylim(rplim)
    ax3.yaxis.set_major_formatter(NullFormatter())
    ax3.yaxis.set_minor_formatter(NullFormatter())
    ax3.set_xticks([10,30,100,300])
    ax3.set_xticklabels(['10','30','100','300'])
    ax3.set_yticks(range(1,6))
    ax3.set_yticklabels('')
    
    fig.subplots_adjust(bottom=.15, top=.97, right=.98, left=.06, wspace=.05)
    if label:
        plt.savefig('plots/planetsample_star_34.png')
    if pltt:
        plt.show()
    plt.close('all')


def plot_planet_population_followup(self, pltt=True, label=False):
    assert hasattr(self, 'isTESSalert')
    ta = (self.isTESSalert==1) & np.in1d(self.disposition_human, range(2)) 
    g0 = self.disposition_human == 0
    g1 = self.disposition_human == 1
    g2 = self.disposition_human == 2
    g25 = self.disposition_human == 2.5

    # light orange, red, yellow, orange
    cols = ['#fe9500','#a80000','#ffff5f','#ff4700','#fe9500']
    marker = ['o','^','p','d']
    Plim = (.5,2e2)
    Jlim = (7,15)
    
    # undetected TOIs
    tois_UD = [203.01,221.01,175.02,175.03]
    assert np.any(np.in1d(tois_UD, self.tois, invert=True)) 
    tics = [259962054,316937670,307210830,307210830]
    udta = np.in1d(self.tics, tics)
    Ps_UD = np.array([52,.624,7.45,2.25])
    rps_UD = np.array([1.,1.67, 1.43,.794])
    mpU,KU,TeqU,SFU,TSMU = estimate_transmission_metric(Ps_UD,rps_UD,
                                                        self.Jmags[udta],
                                                        self.Mss[udta],
                                                        self.Rss[udta],
                                                        self.Teffs[udta])
    _,_,_,BratioU,ESMU = estimate_eclipse_metric(Ps_UD,rps_UD,
                                                 self.Kmags[udta],
                                                 self.Mss[udta],
                                                 self.Rss[udta],
                                                 self.Teffs[udta])
    
    # get barclay planet parameters
    fname = 'input_data/TESStargets/apjsaae3e9t2_mrt.txt'
    p = np.loadtxt(fname, skiprows=54, usecols=(11,12,14,15,16,19,21,22))
    KmagB,JmagB,RsB,MsB,TeffB,detB,PB,rpB = p.T
    rpB += np.random.randn(rpB.size)*.001
    gB = (TeffB <= 4200) & (detB==1) & (rpB < 10)
    KmagB,JmagB,RsB,MsB,TeffB,detB,PB,rpB = p[gB].T
    
    # get known transiting planets (NASA archive Dec 15)
    fname = 'input_data/TESStargets/NASAarchive_transitingMdwarfplanets.csv'
    p = np.loadtxt(fname, delimiter=',', usecols=(2,16,20,24,33,38,44),
                   skiprows=57)
    PT,TeffT,MsT,RsT,rpT,JmagT,KmagT = p.T
    rpT += np.random.randn(rpT.size)*.001
    gT = (TeffT <= 4200) & (rpT < 10)
    PT,TeffT,MsT,RsT,rpT,JmagT,KmagT = p[gT].T

    # compute transmission metrics (Kempton+2018)
    _,Ktoi,_,_,TSMtoi = estimate_transmission_metric(self.Ps[ta],self.rps[ta],
                                                     self.Jmags[ta],
                                                     self.Mss[ta],self.Rss[ta],
                                                     self.Teffs[ta])
    mpB,KB,TeqB,SFB,TSMB = estimate_transmission_metric(PB,rpB,JmagB,
                                                        MsB,RsB,TeffB)
    mpT,KT,TeqT,SFT,TSMT = estimate_transmission_metric(PT,rpT,JmagT,
                                                        MsT,RsT,TeffT)
    mp0,K0,Teq0,SF0,TSM0 = estimate_transmission_metric(self.Ps[g0],
                                                        self.rps[g0],
                                                        self.Jmags[g0],
                                                        self.Mss[g0],
                                                        self.Rss[g0],
                                                        self.Teffs[g0])
    mp1,K1,Teq1,SF1,TSM1 = estimate_transmission_metric(self.Ps[g1],
                                                        self.rps[g1],
                                                        self.Jmags[g1],
                                                        self.Mss[g1],
                                                        self.Rss[g1],
                                                        self.Teffs[g1])
    mp2,K2,Teq2,SF2,TSM2 = estimate_transmission_metric(self.Ps_singletransit[g2],
                                                        self.rps_singletransit[g2],
                                                        self.Jmags[g2],
                                                        self.Mss[g2],
                                                        self.Rss[g2],
                                                        self.Teffs[g2])
    mp25,K25,Teq25,SF25,TSM25 = estimate_transmission_metric(
        self.Ps_singletransit[g25],self.rps_singletransit[g25],
        self.Jmags[g25],self.Mss[g25],self.Rss[g25],self.Teffs[g25])

    # compute emission metrics
    _,_,_,_,ESMtoi = estimate_eclipse_metric(self.Ps[ta],self.rps[ta],
                                             self.Kmags[ta], self.Mss[ta],
                                             self.Rss[ta], self.Teffs[ta])
    _,_,_,BratioB,ESMB = estimate_eclipse_metric(PB,rpB,KmagB,
                                                 MsB,RsB,TeffB)
    _,_,_,BratioT,ESMT = estimate_eclipse_metric(PT,rpT,KmagT,
                                                 MsT,RsT,TeffT)
    _,_,_,Bratio0,ESM0 = estimate_eclipse_metric(self.Ps[g0], self.rps[g0],
                                                 self.Kmags[g0], self.Mss[g0],
                                                 self.Rss[g0], self.Teffs[g0])
    _,_,_,Bratio1,ESM1 = estimate_eclipse_metric(self.Ps[g1], self.rps[g1],
                                                 self.Kmags[g1], self.Mss[g1],
                                                 self.Rss[g1], self.Teffs[g1])
    _,_,_,Bratio2,ESM2 = estimate_eclipse_metric(self.Ps_singletransit[g2],
                                                 self.rps_singletransit[g2],
                                                 self.Kmags[g2], self.Mss[g2],
                                                 self.Rss[g2], self.Teffs[g2])
    _,_,_,Bratio25,ESM25 = estimate_eclipse_metric(self.Ps_singletransit[g25],
                                                   self.rps_singletransit[g25],
                                                   self.Kmags[g25],self.Mss[g25],
                                                   self.Rss[g25], self.Teffs[g25])
    
    
    # plot planet candidates that passed human vetting (and maybe vespa?)
    fig = plt.figure(figsize=(3.3,5))

    # plot K vs P
    ax1 = fig.add_subplot(311)
    ax1.plot(self.Jmags[ta], Ktoi, marker[3], markerfacecolor=cols[3],
             markeredgecolor='k', markeredgewidth=.3, ms=9)
    ax1.plot(self.Jmags[g0], K0, marker[1], ms=5, markerfacecolor=cols[1],
             markeredgecolor='k', markeredgewidth=.5)
    ax1.plot(self.Jmags[g1], K1, marker[0], ms=5, markerfacecolor=cols[1],
             markeredgecolor='k', markeredgewidth=.5)
    ax1.plot(self.Jmags[g2], K2, marker[2], ms=5, markerfacecolor=cols[2],
             markeredgecolor='k', markeredgewidth=.5)
    ax1.plot(self.Jmags[g25], K25, marker[1], ms=5, markerfacecolor=cols[2],
             markeredgecolor='k', markeredgewidth=.5)
    ax1.plot(self.Jmags[udta], KU, marker[3], markerfacecolor=cols[4],
             markeredgecolor='k', markeredgewidth=.3, ms=4)
    ax1.plot(JmagT, KT, 'o', ms=2, c='k', alpha=.4)
    ax1.plot(JmagB, KB, 'v', ms=2, c='k', alpha=.4)

    # plot high Omega values (Eq 24 in Cloutier+2018)
    Omega0 = (self.Omega[g0] > .14*self.Jmags[g0]-.35) & (self.Jmags[g0]<11.7)
    Omega1 = (self.Omega[g1] > .14*self.Jmags[g1]-.35) & (self.Jmags[g1]<11.7)
    Omega2 = (self.Omega[g2] > .14*self.Jmags[g2]-.35) & (self.Jmags[g2]<11.7)
    Omega25 = (self.Omega[g25]>.14*self.Jmags[g25]-.35) & (self.Jmags[g25]<11.7)
    ax1.plot(self.Jmags[g0][Omega0], K0[Omega0], 'o', ms=8,
             markerfacecolor='none', markeredgecolor='k', markeredgewidth=1.2)
    ax1.plot(self.Jmags[g1][Omega1], K1[Omega1], 'o', ms=8,
             markerfacecolor='none', markeredgecolor='k', markeredgewidth=1.2)
    ax1.plot(self.Jmags[g2][Omega2], K2[Omega2], 'o', ms=8,
             markerfacecolor='none', markeredgecolor='k', markeredgewidth=1.2)
    ax1.plot(self.Jmags[g25][Omega25], K25[Omega25], 'o', ms=8,
             markerfacecolor='none', markeredgecolor='k', markeredgewidth=1.2)
    
    #ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_ylabel('Expected RV\nSemi-Amplitude [m/s]', fontsize=8, labelpad=1)
    ax1.set_ylim((.2,100))
    ax1.set_xlim(Jlim)
    ax1.set_xticklabels('')
    ax1.set_yticks([1,10,100])
    ax1.set_yticklabels(['1','10','100'], fontsize=8)
    
    # plot TSM vs P
    ax2 = fig.add_subplot(312)
    ax2.plot(self.Jmags[ta], TSMtoi, marker[3], markerfacecolor=cols[3],
             markeredgecolor='k', markeredgewidth=.3, ms=9)
    ax2.plot(self.Jmags[g0], TSM0, marker[1], ms=5, markerfacecolor=cols[1],
             markeredgecolor='k', markeredgewidth=.5)
    ax2.plot(self.Jmags[g1], TSM1, marker[0], ms=5, markerfacecolor=cols[1],
             markeredgecolor='k', markeredgewidth=.5)
    ax2.plot(self.Jmags[g2], TSM2, marker[2], ms=5, markerfacecolor=cols[2],
             markeredgecolor='k', markeredgewidth=.5)
    ax2.plot(self.Jmags[g25], TSM25, marker[1], ms=5, markerfacecolor=cols[2],
             markeredgecolor='k', markeredgewidth=.5)
    ax2.plot(self.Jmags[udta], TSMU, marker[3], markerfacecolor=cols[4],
             markeredgecolor='k', markeredgewidth=.3, ms=4)
    ax2.plot(JmagT, TSMT, 'o', ms=2, c='k', alpha=.4)
    ax2.plot(JmagB, TSMB, 'v', ms=2, c='k', alpha=.4)
    ax2.axvline(8.1, ls='--', lw=2, color='k')
    ax2.text(8.1, .7, '$\sim$ NIRISS SOSS\n     bright limit', fontsize=7,
             weight='semibold')

    # plot high TSMs values (Eq 24 in Cloutier+2018)
    gTSM0 = _does_exceed_TSM_cutoff(self.rps[g0], TSM0)
    gTSM1 = _does_exceed_TSM_cutoff(self.rps[g1], TSM1)
    gTSM2 = _does_exceed_TSM_cutoff(self.rps[g2], TSM2)
    gTSM25 = _does_exceed_TSM_cutoff(self.rps[g25], TSM25)
    ax2.plot(self.Jmags[g0][gTSM0], TSM0[gTSM0], 'o', ms=8,
             markerfacecolor='none', markeredgecolor='k', markeredgewidth=1.2)
    ax2.plot(self.Jmags[g1][gTSM1], TSM1[gTSM1], 'o', ms=8,
             markerfacecolor='none', markeredgecolor='k', markeredgewidth=1.2)
    ax2.plot(self.Jmags[g2][gTSM2], TSM2[gTSM2], 'o', ms=8,
             markerfacecolor='none', markeredgecolor='k', markeredgewidth=1.2)
    ax2.plot(self.Jmags[g25][gTSM25], TSM25[gTSM25], 'o', ms=8,
             markerfacecolor='none', markeredgecolor='k', markeredgewidth=1.2)

    #ax2.set_xscale('log')
    ax2.set_yscale('log')
    ax2.set_ylabel('Transmission\nSpectroscopy Metric', fontsize=8, labelpad=2)
    ax2.set_ylim((.5,2000))
    ax2.set_xlim(Jlim)
    ax2.set_xticklabels('')
    #ax2.set_yticks([.1,1,10,100])
    #ax2.set_yticklabels(['0.1','1','10','100'], fontsize=8)
    
    # plot ESM vs P
    ax3 = fig.add_subplot(313)
    ax3.plot(self.Jmags[ta], ESMtoi, marker[3], markerfacecolor=cols[3],
             markeredgecolor='k', markeredgewidth=.3, ms=9)
    ax3.plot(self.Jmags[g0], ESM0, marker[1], ms=5, markerfacecolor=cols[1],
             markeredgecolor='k', markeredgewidth=.5)
    ax3.plot(self.Jmags[g1], ESM1, marker[0], ms=5, markerfacecolor=cols[1],
             markeredgecolor='k', markeredgewidth=.5)
    ax3.plot(self.Jmags[g2], ESM2, marker[2], ms=5, markerfacecolor=cols[2],
             markeredgecolor='k', markeredgewidth=.5)
    ax3.plot(self.Jmags[g25], ESM25, marker[1], ms=5, markerfacecolor=cols[2],
             markeredgecolor='k', markeredgewidth=.5)
    ax3.plot(self.Jmags[udta], ESMU, marker[3], markerfacecolor=cols[4],
             markeredgecolor='k', markeredgewidth=.3, ms=4)
    ax3.plot(JmagT, ESMT, 'o', ms=2, c='k', alpha=.4)
    ax3.plot(JmagB, ESMB, 'v', ms=2, c='k', alpha=.4)

    # plot high ESMs values
    gESM0 = _does_exceed_ESM_cutoff(ESM0)
    gESM1 = _does_exceed_ESM_cutoff(ESM1)
    gESM2 = _does_exceed_ESM_cutoff(ESM2)
    gESM25 = _does_exceed_ESM_cutoff(ESM25)
    ax3.plot(self.Jmags[g0][gESM0], ESM0[gESM0], 'o', ms=8,
             markerfacecolor='none', markeredgecolor='k', markeredgewidth=1.2)
    ax3.plot(self.Jmags[g1][gESM1], ESM1[gESM1], 'o', ms=8,
             markerfacecolor='none', markeredgecolor='k', markeredgewidth=1.2)
    ax3.plot(self.Jmags[g2][gESM2], ESM2[gESM2], 'o', ms=8,
             markerfacecolor='none', markeredgecolor='k', markeredgewidth=1.2)
    ax3.plot(self.Jmags[g25][gESM25], ESM25[gESM25], 'o', ms=8,
             markerfacecolor='none', markeredgecolor='k', markeredgewidth=1.2)

    #ax3.set_xscale('log')
    ax3.set_yscale('log')
    ax3.set_xlabel('$J$', fontsize=8, labelpad=0)
    ax3.set_ylabel('Emission\nSpectroscopy Metric', fontsize=8, labelpad=-.5)
    ax3.set_ylim((1e-2,1e3))
    ax3.set_xlim(Jlim)
    #ax3.set_xticklabels('')
    
    fig.subplots_adjust(bottom=.06, top=.98, right=.98, left=.17, hspace=.05)
    if label:
        plt.savefig('plots/planetsample_followup.png')
    if pltt:
        plt.show()
    plt.close('all')


def _does_exceed_TSM_cutoff(rps, TSMs):
    assert rps.size == TSMs.size
    cutoffs = np.zeros(rps.size)
    for i in range(rps.size):
        if rps[i] < 1.5:
            cutoffs[i] = 12
        elif 1.5 < rps[i] < 2.75:
            cutoffs[i] = 92
        elif 2.75 < rps[i] < 4:
            cutoffs[i] = 84
        elif 4 < rps[i] < 10:
            cutoffs[i] = 96

    return TSMs >= cutoffs


def _does_exceed_ESM_cutoff(ESMs):
    ESM_gj1132b = 7.5
    return ESMs >= ESM_gj1132b
    
    
def plot_transit_LCs(self, simname='LC_-00099_8d4',
                     folder='PipelineResults_TIC', pltt=True, label=False):

    g = self.disposition_human >=0
    g = np.where(g)[0][np.argsort(self.tics_candidates[g])]
    Np = g.size
    assert Np == 15
    disp_dict = {0:'pPC', 1:'PC', 2:'ST', 2.5:'pST'}
    
    # black, red, yellow, orange, warm orange
    cols = ['k','#a80000','#ffff5f','#ff4700','#fe9500']
    limits = {12421862:{'xlim':[(-7,7)],'ylim':[(.997,1.003)],'dt':[.38],
                         'Ndec':[2]},
              33734143:{'xlim':[(-7,7)],'ylim':[(.997,1.003)],'dt':[.35],
                         'Ndec':[2]},
              47484268:{'xlim':[(-7,7)],'ylim':[(.983,1.015)],'dt':[.3],
                         'Ndec':[2]},
              49678165:{'xlim':[(-13,13)],'ylim':[(.975,1.025)],'dt':[.3],
                         'Ndec':[2]},
              55652063:{'xlim':[(-11,11)],'ylim':[(.995,1.005)],'dt':[.3],
                         'Ndec':[2]},
              92444219:{'xlim':[(-10,10)],'ylim':[(.985,1.014)],'dt':[.35],
                         'Ndec':[2]},
              100103200:{'xlim':[(-7,7)],'ylim':[(.997,1.003)],'dt':[.23],
                         'Ndec':[2]},
              100103201:{'xlim':[(-12,12),(-8,8)],
                         'ylim':[(.994,1.006),(.995,1.005)],'dt':[.3,.3],
                         'Ndec':[2,2]},
              141708335:{'xlim':[(-7,7)],'ylim':[(.99,1.01)],'dt':[.2],
                         'Ndec':[2]},
              206660104:{'xlim':[(-10,10)],'ylim':[(.997,1.003)],'dt':[.25],
                         'Ndec':[2]},
              231279823:{'xlim':[(-8,8)],'ylim':[(.997,1.0025)],'dt':[.2],
                         'Ndec':[2]},
              231702397:{'xlim':[(-6,6)],'ylim':[(.982,1.018)],'dt':[.2],
                         'Ndec':[2]},
              234994474:{'xlim':[(-4,4)],'ylim':[(.998,1.002)],'dt':[.2],
                         'Ndec':[2]},
              235037759:{'xlim':[(-6,6)],'ylim':[(.91,1.08)],'dt':[.3],
                         'Ndec':[2]},
              238027971:{'xlim':[(-9,9)],'ylim':[(.995,1.004)],'dt':[.55],
                         'Ndec':[2]},
              260004324:{'xlim':[(-4,4)],'ylim':[(.997,1.003)],'dt':[.4],
                         'Ndec':[2]},
              262530407:{'xlim':[(-3.5,3.5)],'ylim':[(.997,1.003)],'dt':[.2],
                         'Ndec':[2]},
              278661431:{'xlim':[(-12,12)],'ylim':[(.983,1.017)],'dt':[.4],
                         'Ndec':[2]},
              279574462:{'xlim':[(-12,12)],'ylim':[(.94,1.06)],'dt':[.5],
                         'Ndec':[2]},
              303586421:{'xlim':[(-10,10)],'ylim':[(.982,1.017)],'dt':[.4],
                         'Ndec':[2]},
              305048087:{'xlim':[(-7,7)],'ylim':[(.98,1.02)],'dt':[.3],
                         'Ndec':[2]},
              307210830:{'xlim':[(-6,6)],'ylim':[(.996,1.004)],'dt':[.18],
                         'Ndec':[2]},
              415969908:{'xlim':[(-8,8),(-12,12)],
                         'ylim':[(.993,1.007),(.993,1.007)],'dt':[.25,.5],
                         'Ndec':[2,2]},
              441026957:{'xlim':[(-2,2)],'ylim':[(.998,1.002)],'dt':[.2],
                         'Ndec':[2]},
              441056702:{'xlim':[(-7,7)],'ylim':[(.993,1.007)],'dt':[.3],
                         'Ndec':[2]}}
    
    fig = plt.figure(figsize=(6.5,6))
    ticsDONE, subplot_ind = np.zeros(0), 1
    for i in range(Np):

        # define planet parameters
        tic = self.tics[g][i]
	tic_candidate = self.tics_candidates[g][i]
        disposition = self.disposition_human[g][i]
        print i, tic
        P = self.Ps[g][i] if disposition < 2 else self.Ps_singletransit[g][i]
        aRs = self.aRss[g][i] if disposition<2 else self.aRs_singletransit[g][i]
        rpRs=self.rpRss[g][i] if disposition<2 else self.rpRs_singletransit[g][i]
        inc = self.incs[g][i] if disposition<2 else self.inc_singletransit[g][i]
        T0 = self.T0s[g][i]
        isalert = self.isTESSalert[g][i]

        f2add = .007 if tic == 279574462 else 0
            
        # compute transit model
        d = loadpickle('%s/TIC_%i/%s'%(folder, tic, simname))
        func = llnl.transit_model_func_curve_fit(d.u1, d.u2)
        fmodel = func(d.bjd, P, T0, aRs, rpRs, inc)

        # phase-fold the LC
        phase = foldAt(d.bjd, P, T0)
        phase[phase > .5] -= 1
        phase_hrs = phase*P*24  # phase in hours
        s = np.argsort(phase_hrs)
        PCind = np.in1d(ticsDONE, tic).sum()
        dt = limits[tic]['dt'][PCind]
        tb, fb, efb = llnl.boxcar(phase_hrs[s], d.fcorr[s], d.ef[s], dt=dt)

        ax = fig.add_subplot(4,4,subplot_ind)
        subplot_ind += 1
        ax.plot(phase_hrs, d.fcorr+f2add, 'o', ms=1, c=cols[0], alpha=.1)

        # plot color and marker depdending on the source
        if isalert:
            m, ms, c = 'd', 4, cols[4]
        if disposition == 0:
            m, ms, c = '^', 4, cols[3]
        elif disposition == 1:
            m, ms, c = 'o', 3, cols[1]
        elif disposition == 2.5:
            m, ms, c = '^', 4, cols[2]            
        else:
            m, ms, c = 'p', 4, cols[2]
        ax.plot(tb, fb+f2add, m, ms=ms, markerfacecolor=c, markeredgecolor='k',
                markeredgewidth=.2)
        ax.plot(phase_hrs[s], fmodel[s], '-', c=cols[0])
        ax.text(.95, .92, disp_dict[disposition], fontsize=7,
                horizontalalignment='right', verticalalignment='top',
                transform=ax.transAxes, weight='semibold')
        if isalert:
            toi = '%.2f'%self.tois[g][i] if self.tois[g][i] > 0 else ''
            ax.text(.05, .92, toi, fontsize=7,
                    horizontalalignment='left', verticalalignment='top',
                    transform=ax.transAxes, weight='semibold')
        ax.set_title('%.2f'%tic_candidate, fontsize=8,
                     weight='semibold', y=.95)
        ax.set_xlim(limits[tic]['xlim'][PCind])
        ax.set_ylim(limits[tic]['ylim'][PCind])
        yticks = ax.get_yticks()
        ax.set_yticks(yticks)
        ax.set_yticklabels(_normF2percent(yticks, limits[tic]['Ndec'][PCind]))
        ax.set_ylim(limits[tic]['ylim'][PCind])

        # this TIC is done now
        ticsDONE = np.append(ticsDONE, tic)
        
        
    fig.text(.01, .5, 'Transit depth [%]', rotation=90,
             verticalalignment='center')
    fig.text(.5, .01, 'Hours from mid-transit', horizontalalignment='center')
    fig.subplots_adjust(bottom=.07, left=.09, right=.98, top=.97,
                        wspace=.26, hspace=.28)
    if label:
        plt.savefig('plots/transitLC_PCs.png')
    if pltt:
        plt.show()
    plt.close('all')



def _normF2percent(Farr, Ndec):
    return np.round(1e2*(np.array(Farr)-1.), Ndec)


# try: d=loadpickle('PipelineResults_TIC/TIC_25200252/LC_-00099_12d5')
def plot_flare_LC(self, medkernel=9, Nsig_flare=8, flare_dur_days=30./60/24,
                  pltt=True, label=False):

    cols = ['k','#','#','#','#f36300']
    t0 = 2457000
    bjdlim = self.bjd.min()-t0-.3, self.bjd.max()-t0+.3
    tb, fb, ef = llnl.boxcar(self.bjd, medfilt(self.fcorr,medkernel),
                             self.ef, dt=30./24/60)
    inflare=self.fcorr>=np.median(self.fcorr)+Nsig_flare*llnl.MAD1d(self.fcorr)
    Nflares = (np.diff(np.where(np.diff(self.bjd[inflare]) > flare_dur_days)[0]) > 1).sum()
    
    fig = plt.figure(figsize=(3.7,2.7))
    ax = fig.add_subplot(111)    
    ax.plot(self.bjd-t0, self.fcorr, '.', c=cols[0], ms=1, alpha=.6)
    ax.plot(self.bjd[inflare]-t0, self.fcorr[inflare], '.', ms=6, alpha=1,
            markerfacecolor=cols[4], markeredgecolor='k', markeredgewidth=.4)

    ax1 = fig.add_axes([.54,.6,.33,.33])
    t0_flare = 1332.212
    dxi, dxf = -70, 120
    yi, yf = .95, 1.8
    ax1.plot((self.bjd-t0-t0_flare)*24*60, self.fcorr, '.', ms=3, alpha=.9,
             c=cols[0])
    ax1.plot((self.bjd-t0-t0_flare)*24*60, self.fcorr, '-', lw=.5, c=cols[0])
    ax1.plot((self.bjd[inflare]-t0-t0_flare)*24*60, self.fcorr[inflare], '.',
             ms=3, alpha=.9, markerfacecolor=cols[4], markeredgecolor='k',
             markeredgewidth=.2)
    ax.plot([t0_flare,1338.8],[yi,1.42], '--', c='k', lw=.5)
    ax.plot([t0_flare,1338.8],[yf,1.76], '--', c='k', lw=.5)
    ax.fill_between(np.array([dxi,dxf])/60./24+t0_flare, [yi,yi], [yf,yf],
                    color=cols[0], alpha=.2)
    ax1.set_facecolor('k')
    ax1.patch.set_alpha(.1)

    # flag flares
    bjd_flares = self.bjd[inflare][np.where(np.diff(self.bjd[inflare]) > flare_dur_days)[0]][np.where(np.diff(np.where(np.diff(self.bjd[inflare]) > flare_dur_days)[0])>1)[0]+1]
    for b in bjd_flares:
        ax.axvline(b-2457000, lw=1.5, color=cols[4], ymax=.03)
    
    ax1.set_xlim((dxi, dxf))
    ax1.set_ylim((yi,yf))
    ax.set_xlim(bjdlim)
    ax.set_ylim((yi,yf))
    ax1.set_ylabel('De-trended Flux', fontsize=7, labelpad=8, rotation=270)
    ax1.set_xlabel('Time from flare\npeak [minutes]', fontsize=7, labelpad=.1)
    ax1.set_xticks(range(-60,dxf+1,60))
    ax1.set_xticklabels(['%i'%i for i in range(-60,dxf+1,60)], fontsize=7)
    ax1.set_yticks(np.arange(1,1.9,.4))
    ax1.set_yticklabels(['%.1f'%i for i in np.arange(1,1.9,.4)], fontsize=7)
    ax1.yaxis.tick_right()
    ax1.yaxis.set_label_position("right")
    ax1.yaxis.set_ticks_position('both')

    ax.set_xlabel('Time [BJD - 2,4570,000]', fontsize=8)
    ax.set_ylabel('De-trended Flux', fontsize=8)

    fig.subplots_adjust(bottom=.14, top=.97, right=.98, left=.12)
    if label:
        plt.savefig('plots/flareLC_%i.png'%self.tic)
    if pltt:
        plt.show()
    plt.close('all')



def plot_WTF_Ls(pltt=True, label=False):

    cols = ['k','#fe9500']
    t0 = 2457000
    
    self6 = loadpickle('PipelineResults_TIC/TIC_63037741/LC_-00099_8d4')
    bjdlim6 = self6.bjd.min()-t0-.3, self6.bjd.max()-t0+.3
    bjdlim61 = 1340.6-.3, 1340.6+.3
    bjdlim62 = 1352.99-.14, 1352.99+.14
    ylim6 = .93, 1.04
    ylim61 = .985, 1.039
    ylim62 = .93, 1.02
    
    self4 = loadpickle('PipelineResults_TIC/TIC_434105091/LC_-00099_8d4')
    t04 = 1364.96
    bjdlim4 = self4.bjd.min()-t0-.3, self4.bjd.max()-t0+.3
    bjdlim41 = t04-.5, t04+.5
    ylim4 = .91, 1.02
    
    fig = plt.figure(figsize=(3.3,6))
    gs = gridspec.GridSpec(4,2)
    ax1 = plt.subplot(gs[0,:])
    ax2 = plt.subplot(gs[1,:1])
    ax3 = plt.subplot(gs[1,1:])
    ax4 = plt.subplot(gs[2,:])
    ax5 = plt.subplot(gs[3,:])

    # plot first full LC
    ax1.plot(self6.bjd-t0, self6.f, 'o', c=cols[0], ms=1, alpha=.3)    
    #ax1.plot(self6.bjd-t0, self6.f, '-', c=cols[0], lw=.6, alpha=1)    
    ax1.text(.06,.11,'TIC %i'%self6.tic, fontsize=7, weight='semibold',
             transform=ax1.transAxes)
    ax1.text(1339, 1.025, 'a', weight='semibold', fontsize=9)
    ax1.text(1351.5, 1.025, 'b', weight='semibold', fontsize=9)
    ax1.fill_between(bjdlim61, np.repeat(ylim6[0],2), np.repeat(ylim6[1],2),
                     color='k', alpha=.1)
    ax1.fill_between(bjdlim62, np.repeat(ylim6[0],2), np.repeat(ylim6[1],2),
                     color='k', alpha=.1)
    ax1.set_ylim(ylim6)
    ax1.set_xlim(bjdlim6)
    ax1.set_ylabel('Normalized\nde-trended flux', fontsize=8)
    ax1.set_xlabel('Time [BJD-2,457,000]', fontsize=8, labelpad=0)
    ax1.set_xticks(range(1325, 1355, 5))
    ax1.set_xticklabels(['%i'%i for i in range(1325, 1355, 5)], fontsize=8)
    ax1.set_yticks(np.arange(.94,1.05,.02))
    ax1.set_yticklabels(['%.2f'%i for i in np.arange(.94,1.05,.02)],
                        fontsize=8)
    
    # zoom-in first feature
    ax2.plot((self6.bjd-t0-1340.6)*24, self6.f, 'o',
             markerfacecolor=cols[1],
             markeredgecolor='k', markeredgewidth=.3, ms=1.5, alpha=1)    
    ax2.plot((self6.bjd-t0-1340.6)*24, self6.f, '-', c=cols[0], lw=.6,
             alpha=1)
    ax2.text(.1, .85, 'a', weight='semibold', fontsize=9,
             transform=ax2.transAxes)
    ax2.set_ylim(ylim61)
    ax2.set_xlim((-.3*24,.3*24))
    ax2.set_ylabel('Flux depth [%]', fontsize=8)
    ax2.set_xlabel('$\Delta t$ [hrs]', fontsize=8, labelpad=0)
    ax2.set_yticks(np.arange(.99,1.035,.01))
    ax2.set_yticklabels(['%i'%(i*1e2-1e2) for i in np.arange(.99,1.035,.01)],
                        fontsize=8)
    ax2.set_xticks(range(-5,6,5))
    ax2.set_xticklabels(['%i'%i for i in range(-5,6,5)], fontsize=8)
    ax2.set_facecolor('k')
    ax2.patch.set_alpha(.06)

    # zoom-in on second feature
    ax3.plot((self6.bjd-t0-1352.99)*24, self6.f, 'o',
             markerfacecolor=cols[1],
             markeredgecolor='k', markeredgewidth=.3, ms=1.5, alpha=1)    
    ax3.plot((self6.bjd-t0-1352.99)*24, self6.f, '-', c=cols[0], lw=.6,
             alpha=1)
    ax3.text(.1, .15, 'b', weight='semibold', fontsize=9,
             transform=ax3.transAxes)
    ax3.set_ylim(ylim62)
    ax3.set_xlim((-.14*24,.14*24))
    ax3.set_xlabel('$\Delta t$ [hrs]', fontsize=8, labelpad=0)
    ax3.set_yticks(np.arange(.94,1.025,.02))
    ax3.set_yticklabels(['%i'%(i*1e2-1e2) for i in np.arange(.94,1.025,.02)],
                        fontsize=8)
    ax3.set_xticks(np.arange(-2.5,3,2.5))
    ax3.set_xticklabels(['%.1f'%i for i in np.arange(-2.5,3,2.5)], fontsize=8)
    ax3.set_facecolor('k')
    ax3.patch.set_alpha(.06)

    # plot first full LC
    ax4.plot(self4.bjd-t0, self4.f, 'o', c=cols[0], ms=1, alpha=.3)    
    ax4.text(.6,.11,'TIC %i'%self4.tic, fontsize=7, weight='semibold',
             transform=ax4.transAxes)
    ax4.text(1363.2, .96, 'c', weight='semibold', fontsize=9)
    ax4.fill_between(bjdlim41, np.repeat(ylim4[0],2), np.repeat(ylim4[1],2),
                     color='k', alpha=.1)
    ax4.set_ylim(ylim4)
    ax4.set_xlim(bjdlim4)
    ax4.set_xlabel('Time [BJD-2,457,000]', fontsize=8, labelpad=0)
    ax4.set_ylabel('Normalized\nde-trended flux', fontsize=8)
    ax4.set_xticks(range(1355,1381,5))
    ax4.set_xticklabels(['%i'%i for i in range(1355,1381,5)], fontsize=8)
    ax4.set_yticks(np.arange(.92,1.03,.02))
    ax4.set_yticklabels(['%.2f'%i for i in np.arange(.92,1.03,.02)],
                        fontsize=8)
    
    # zoom-in on feature
    ax5.plot((self4.bjd-t0-1364.96)*24, self4.f, 'o',
             markerfacecolor=cols[1],
             markeredgecolor='k', markeredgewidth=.3, ms=1.5, alpha=1)    
    ax5.plot((self4.bjd-t0-1364.96)*24, self4.f, '-', c=cols[0], lw=.6,
             alpha=1)
    ax5.text(.1, .15, 'c', weight='semibold', fontsize=9,
             transform=ax5.transAxes)
    ax5.set_ylim(ylim4)
    ax5.set_xlim((-.5*24,.5*24))
    ax5.set_xlabel('$\Delta t$ [hrs]', fontsize=8, labelpad=0)
    ax5.set_ylabel('Flux depth [%]', fontsize=8)
    ax5.set_yticks(np.arange(.92,1.02,.02))
    ax5.set_yticklabels(['%i'%(i*1e2-1e2) for i in np.arange(.92,1.02,.02)],
                        fontsize=8)
    ax5.set_xticks(range(-10,11,5))
    ax5.set_xticklabels(['%i'%i for i in range(-10,11,5)], fontsize=8)
    ax5.set_facecolor('k')
    ax5.patch.set_alpha(.06)


    fig.subplots_adjust(bottom=.05, top=.98, right=.98, left=.19, hspace=.3)
    if label:
        plt.savefig('plots/wtfLCs.png')
    if pltt:
        plt.show()
    plt.close('all')


def plot_old_GAIA_Rs(self, pltt=True, label=False):
    # get old stellar radii
    fname = 'input_data/TESStargets/TICv7_Mdwarfsv1.csv'
    TICs, Rss, e_Rss = np.loadtxt(fname, delimiter=',', skiprows=5,
                                  usecols=(0,70,71)).T

    
    cols = ['#ff4700','#a80000','#ffff5f']

    # match stars
    _,g = np.unique(self.tics, return_index=True)
    s = np.argsort(self.tics[g])
    tics1, Rss1, e_Rss1 = self.tics[g][s], self.Rss[g][s], self.ehi_Rss[g][s]
    g = np.in1d(TICs, self.tics)
    s = np.argsort(TICs[g])
    tics0, Rss0, e_Rss0 = TICs[g][s].astype(int), Rss[g][s], e_Rss[g][s]
    assert tics0.size == tics1.size
    
    # plot new Rs vs old Rs
    fig = plt.figure(figsize=(3.3,3.3))
    ax = fig.add_subplot(111)
    #ax.errorbar(Rss0, Rss1, xerr=e_Rss0, yerr=e_Rss1, fmt='ko', ms=1,
    #            ecolor='k', elinewidth=.9, alpha=.4)
    print np.median((Rss1-Rss0)/Rss1)
    H, x_edges, y_edges = np.histogram2d(Rss0, Rss1, bins=np.linspace(.1,.7,30))
    cax = ax.pcolormesh(x_edges, y_edges, H.T,
                        cmap=truncate_cmap(plt.get_cmap('hot_r'),0,.95))
    cbar_axes = fig.add_axes([.13,.09,.84,.04])
    cbar = fig.colorbar(cax, cax=cbar_axes, orientation='horizontal')
    cbar.set_label('Number of TICs', fontsize=8, labelpad=.4)
    ax.plot([.1,.75], [.1,.75], 'k--', lw=2)

    ax.set_xlabel('TIC-7 Stellar Radii [R$_{\odot}$]', fontsize=8,
                  labelpad=1.5)
    ax.set_ylabel('GAIA-derived Stellar Radii [R$_{\odot}$]', fontsize=8,
                  labelpad=2)
    ax.set_xlim((.1,.7))
    ax.set_ylim((.1,.7))

    # add histogram of fractional uncertainties
    ax1 = fig.add_axes([.7,.35,.23,.23])
    #ratio_bins = np.linspace(0,.25,50)
    ratio_bins = np.logspace(0,np.log10(30),30)
    ax1.hist(e_Rss0/Rss0*1e2, bins=ratio_bins, histtype='step', color=cols[0],
             lw=.8, label='ticv7', log=1)
    ax1.hist(e_Rss0/Rss0*1e2, bins=ratio_bins, alpha=.5, color=cols[0], log=1)

    ax1.hist(e_Rss1/Rss1*1e2, bins=ratio_bins, histtype='step', color=cols[1],
             lw=.8, label='gaia', log=1)
    ax1.hist(e_Rss1/Rss1*1e2, bins=ratio_bins, alpha=.5, color=cols[1], log=1)

    print np.median(e_Rss1/Rss1*1e2)
    
    ax1.text(.58, 1.03, 'TIC-7', weight='semibold', fontsize=6, color=cols[0],
             transform=ax1.transAxes)
    ax1.arrow(.74, 1, .05, -.3, transform=ax1.transAxes, head_width=.02,
              color='k')
    ax1.text(.1, 1.03, 'GAIA', weight='semibold', fontsize=6, color=cols[1],
             transform=ax1.transAxes)
    ax1.arrow(.27, 1, .05, -.2, transform=ax1.transAxes, head_width=.02,
              color='k')
    
    ax1.set_xscale('log')
    ax1.set_xlim((ratio_bins[0],ratio_bins[-1]))
    ax1.set_ylim((.9,2e3))
    ax1.set_xlabel('Fractional R$_s$\nuncertainty [%]', fontsize=6, labelpad=0)
    ax1.set_ylabel('N', fontsize=6, labelpad=0)
    ax1.set_xticks([1,3,10,30])
    ax1.set_xticklabels(['%i'%i for i in [1,3,10,30]], fontsize=6)
    ax1.set_yticks([1,10,100,1000])
    ax1.set_yticklabels(['10$^%i$'%i for i in range(4)], fontsize=6)

    fig.subplots_adjust(bottom=.24, right=.97, top=.97, left=.13)
    if label:
        plt.savefig('plots/stellarradii.png')
    if pltt:
        plt.show()
    plt.close('all')


def plot_gaia_map_FPs(self, pltt=True, label=False):
    '''For TICs with PCs given a high FPP by vespa, query gaia and 
    plot the map to see if there are nearby sources.'''
    g = np.in1d(self.disposition_human, [-8,-7,-6,0,1,2,2.5])
    tics = self.tics[g]
    tics,inds = np.unique(tics, return_index=True)
    disps = self.disposition_human[g][inds]
    s = np.array([0,7,9,10,13,14,18,19,20,21,1,4,15,2,3,5,6,8,12,16,17,11])
    tics, disps = tics[s], disps[s]
    disp_dict = {-8:'FP',-7:'BEB2',-6:'BEB',0:'pPC',2.5:'pST',0:'pPC',1:'PC',
                 2:'ST'}
    NFPs = tics.size
    
    # black, dark orange, light orange, red, yellow
    cols = ['k','#ff4700','#fe9500','#a80000','#ffff5f']
    vmin, vmax = 9.5, 20
    Npix, pixscale = 5, 21  # pixscale in arcseconds per pixel
    ralim_arcsec = -Npix*pixscale-10, Npix*pixscale+10
    declim_arcsec = -Npix*pixscale-10, Npix*pixscale+10
    
    fig = plt.figure(figsize=(6.5,7.5))
    subplot_ind = 1
    for i in range(NFPs):

        # query GAIA around the input source
        g = self.tics == tics[i]
        tic, ra, dec = self.tics[g][0], self.ras[g][0], self.decs[g][0]

        ras, e_ras, decs, e_decs, Gmags, radius_arcsec = \
                                            query_nearby_gaia(tic, ra, dec,
                                                              Npixsearch=Npix,
                                                              pltt=False)

        try:
            FPP=np.loadtxt('PipelineResults_TIC/TIC_%i/FPPs_nobin.txt'%tics[i])
        except IOError:
            FPP = np.nan
        print tics[i], FPP, Gmags.min(), Gmags.max()

        ras *= 36e2
        decs *= 36e2
        ra *= 36e2
        dec *= 36e2

        ax = fig.add_subplot(5,5,subplot_ind)
        subplot_ind += 1
        #ax.plot(0, 0, 'x', ms=6, markerfacecolor='k', markeredgewidth=3,
        #        markeredgecolor='k')
        Gmag0 = 15
        s = np.argsort(Gmags)
        cax = ax.scatter(ras[s]-ra, decs[s]-dec, c=Gmags[s],
                         s=30*10**(-.4*(Gmags[s]-Gmag0)),
                         vmin=vmin, vmax=vmax, edgecolor='k', alpha=1,
                         cmap=truncate_cmap(plt.get_cmap('hot_r'),.1,.9))
        # plot edges such they overlap
        ax.scatter(ras[s]-ra, decs[s]-dec, s=30*10**(-.4*(Gmags[s]-Gmag0)),
                   vmin=vmin, vmax=vmax, edgecolor='k', facecolors='none',
                   alpha=1)
        #ax.plot(ras-ra, decs-dec, 'ko', ms=3)
        
        cbar_axes = fig.add_axes([.48,.16,.50,.04])
        cbar = fig.colorbar(cax, cax=cbar_axes, orientation='horizontal')
        cbar.set_label('$G$', labelpad=-2)
    
        if subplot_ind in range(2,19):
            ax.set_xticklabels('')
        if subplot_ind in [3,4,5,6,8,9,10,11,13,14,15,16,18,19,20,21,23]:
            ax.set_yticklabels('')
        
        # plot FWHM
        fname = 'PipelineResults_TIC/TIC_%i/FWHMs_arcsec.npy'%tic
        fwhm = np.nanmedian(np.load(fname))
        circ = patches.Circle((0,0), radius=0.5*fwhm, ls='--',
                              edgecolor='k', lw=2, facecolor='none')
        ax.add_patch(circ)
        ax.text(.06, .89, "%.1f''"%fwhm, transform=ax.transAxes,
                fontsize=8, weight='semibold')
        disp_label = disp_dict[disps[i]] if tic != 415969908 \
                     else 'PC+\nST'
        ax.text(.94, .97, '%s'%disp_label, transform=ax.transAxes,
                fontsize=8, weight='semibold', horizontalalignment='right',
                verticalalignment='top')
        ax.text(.06, .07, string.ascii_lowercase[i], transform=ax.transAxes,
                fontsize=8, weight='semibold')
        
        # plot TESS pixel
        if subplot_ind == NFPs+1:
            ax1 = fig.add_subplot(5,5,subplot_ind+1)
            rect = patches.Rectangle((-16,-80), pixscale, pixscale, lw=2,
                                     edgecolor='k', facecolor='k',
                                     alpha=.5)
            ax1.add_patch(rect)
            rect = patches.Rectangle((-16,-80), pixscale, pixscale, lw=2,
                                     edgecolor='k', facecolor='none')
            ax1.add_patch(rect)
            ax1.text(-6,-106, "TESS pixelscale = 21 arcsec", weight='semibold',
                     fontsize=8, horizontalalignment='center')
            ax1.set_ylim(declim_arcsec)
            ax1.set_xlim(ralim_arcsec)
            ax1.axis('off')
            
        ax.set_title('TIC %i'%tic, fontsize=8, weight='semibold', y=.96)        
        ax.set_ylim(declim_arcsec)
        ax.set_xlim(ralim_arcsec)


    fig.text(.01, .5, r'$\Delta \delta$ (J2000) [arcsec]', rotation=90,
             verticalalignment='center')
    fig.text(.5, .01, r'$\Delta \alpha$ (J2000) [arcsec]',
             horizontalalignment='center')
    fig.subplots_adjust(bottom=.06, left=.1, right=.99, top=.97,
                        wspace=.07, hspace=.13)
    if label:
        plt.savefig('plots/GAIAFPs.png')
    if pltt:
        plt.show()
    plt.close('all')



def plot_TIC300_lightcurvesOLD(pltt=True, label=False):

    # black, red, orange, warm orange
    cols = ['k','#a80000','#ff4700','#fe9500']
    t0 = 2457000
    
    fig = plt.figure(figsize=(3.3,4))
    ax1 = fig.add_subplot(111)

    # plot TIC 100103200 (PC favoured)
    d = loadpickle('PipelineResults_TIC/TIC_100103200/LC_-00099_8d4')
    P,T0, = d.params_optimized[0,:2]
    phase = foldAt(d.bjd, P, T0)
    phase[phase>.5] -= 1
    s = np.argsort(phase)
    tb, fb, efb = llnl.boxcar(phase[s], d.fcorr[s], d.ef[s], dt=.0015)
    ax1.plot(phase, d.fcorr, 'o', c=cols[0], alpha=.3, ms=1)
    ax1.plot(tb, fb, 'o', markerfacecolor=cols[1], markeredgecolor='k', ms=4)
    ax1.plot(phase[s], d.fmodels[0,s], '-', c=cols[1])
    snr = d.transit_condition_values[2,2]
    labelbase = 1.001
    dym = .001
    ax1.text(.021, labelbase-dym*2, 'S/N$_{transit}$ = %.1f'%snr,
             weight='semibold', fontsize=7)
    ax1.text(.021, labelbase, '%i.01'%d.tic, weight='semibold',
             fontsize=7)
    ax1.text(.021, labelbase-dym, 'P = %.1f days'%P, weight='semibold',
             fontsize=7)

    # plot TIC 100103201 (BEB favoured)
    d = loadpickle('PipelineResults_TIC/TIC_100103201/LC_-00099_8d4')
    i = 0
    P,T0, = d.params_optimized[i,:2]
    phase = foldAt(d.bjd, P, T0)
    phase[phase>.5] -= 1
    s = np.argsort(phase)
    tb, fb, efb = llnl.boxcar(phase[s], d.fcorr[s], d.ef[s], dt=.0015)
    dy = .007
    ax1.plot(phase, d.fcorr-dy, 'o', c=cols[0], alpha=.3, ms=1)
    ax1.plot(tb, fb-dy, 'o', markerfacecolor=cols[2], markeredgecolor='k', ms=4)
    ax1.plot(phase[s], d.fmodels[i,s]-dy, '-', c=cols[2])
    snr = d.transit_condition_values[i,2]
    ax1.text(.021, labelbase-dym*2-dy, 'S/N$_{transit}$ = %.1f'%snr,
             weight='semibold', fontsize=7)
    ax1.text(.021, labelbase-dy, '%i.01'%d.tic, weight='semibold',
             fontsize=7)
    ax1.text(.021, labelbase-dym-dy, 'P = %.1f days'%P, weight='semibold',
             fontsize=7)

    # plot 100103201.02
    i = 2
    P,T0, = d.params_optimized[i,:2]
    phase = foldAt(d.bjd, P, T0)
    phase[phase>.5] -= 1
    s = np.argsort(phase)
    tb, fb, efb = llnl.boxcar(phase[s], d.fcorr[s], d.ef[s], dt=.0015)
    ax1.plot(phase, d.fcorr-dy*2, 'o', c=cols[0], alpha=.3, ms=1)
    ax1.plot(tb, fb-dy*2, 'o', markerfacecolor=cols[3], markeredgecolor='k',
             ms=4)
    ax1.plot(phase[s], d.fmodels[i,s]-dy*2, '-', c=cols[3])
    snr = d.transit_condition_values[i,2]
    ax1.text(.021, labelbase-dym*2-dy*2, 'S/N$_{transit}$ = %.1f'%snr,
             weight='semibold', fontsize=7)
    ax1.text(.021, labelbase-dy*2, '%i.02'%d.tic, weight='semibold',
             fontsize=7)
    ax1.text(.021, labelbase-dym-dy*2, 'P = %.1f days'%P, weight='semibold',
             fontsize=7)

    ax1.set_xlim((-.02,.02))
    ax1.set_ylim((.982,1.003))
    ax1.set_ylabel('Normalized Flux + offset', fontsize=8, labelpad=0)
    ax1.set_xlabel('Orbital Phase', fontsize=8)
    

    fig.subplots_adjust(left=.19, right=.73, top=.98, bottom=.1)
    if label:
        plt.savefig('plots/TIC300_LC.png')
    if pltt:
        plt.show()
    plt.close('all')



    
def plot_TIC300_lightcurves(pltt=True, label=False):

    # black, red, orange, warm orange
    cols = ['k','#a80000','#ff4700','#fe9500','#ffff5f']
    t0 = 2457000
    
    fig = plt.figure(figsize=(3.3,3))
    ax1 = fig.add_subplot(111)

    # plot TIC 100103200 (PC favoured)
    d = loadpickle('PipelineResults_TIC/TIC_100103200/LC_-00099_8d4')
    bjdlim = d.bjd.min()-t0-.3, d.bjd.max()-t0+.3
    P,T0,_,D = d.params_guess[0]
    print P
    phase = foldAt(d.bjd, P, T0)
    phase[phase > .5] -= 1
    g = (phase*P >= -D) & (phase*P <= D)
    tb, fb, eb = llnl.boxcar(d.bjd, d.fcorr, d.ef, dt=30/24./60)
    ax1.plot(d.bjd-t0, d.fcorr, 'o', c=cols[0], alpha=.2, ms=1)
    ax1.plot(d.bjd[g]-t0, d.fcorr[g], 'o', markerfacecolor=cols[1], ms=3,
             markeredgecolor='k', markeredgewidth=.2)
    ax1.plot(tb-t0, fb, 'w-', lw=.3)
    ax1.text(1355, .9951, '100103200.01\n(P = %.1f days)'%P, color=cols[1],
             weight='bold', fontsize=5)


    # plot TIC 10010321 (BEB favoured)
    d = loadpickle('PipelineResults_TIC/TIC_100103201/LC_-00099_8d4')
    P,T0,_,D = d.params_guess[0]
    print P
    phase = foldAt(d.bjd, P, T0)
    phase[phase > .5] -= 1
    g = (phase*P >= -D) & (phase*P <= D)
    dy = .01
    ax1.plot(d.bjd-t0, d.fcorr-dy, 'o', c=cols[0], alpha=.3, ms=1)
    ax1.plot(d.bjd[g]-t0, d.fcorr[g]-dy, 'o', markerfacecolor=cols[2], ms=3,
             markeredgecolor='k', markeredgewidth=.2)
    ax1.text(1355., .9847, '100103201.01\n(P = %.1f days)'%P, color=cols[2],
             weight='bold', fontsize=5)
    ax1.plot(tb-t0, fb, 'w-', lw=.3)


    # single transit event
    P,T0,_,D = d.params_guess[2]
    P = 38.32060422 # single transit value
    print P
    phase = foldAt(d.bjd, P, T0)
    phase[phase > .5] -= 1
    g = (phase*P >= -D) & (phase*P <= D)
    dy = .01
    tb, fb, eb = llnl.boxcar(d.bjd, d.fcorr, d.ef, dt=30/24./60)
    ax1.plot(d.bjd[g]-t0, d.fcorr[g]-dy, 'o', markerfacecolor=cols[3], ms=3,
             markeredgecolor='k', markeredgewidth=.2)
    ax1.text(1371., .9948, '100103201.02\n(P = %.1f days)'%P, color=cols[3],
             weight='bold', fontsize=5)
    ax1.plot(tb-t0, fb-dy, 'w-', lw=.3)


    
    # mark responses in both light curves
    bjds = [1361.6, 1364, 1371.1, 1376, 1380.1]
    for i in bjds:
        ax1.axvline(i, color=cols[4], lw=3, ymax=.05)
        ax1.axvline(i, color=cols[4], lw=3, ymin=.95)
        ax1.axvline(i, color=cols[0], lw=1, ymax=.05)
        ax1.axvline(i, color=cols[0], lw=1, ymin=.95)
    
    ax1.set_xlim(bjdlim)
    ax1.set_ylim((.983,1.004))
    ax1.set_ylabel('Normalized Flux + offset', fontsize=8, labelpad=1)
    ax1.set_xlabel('Time [BJD - 2,457,000]', fontsize=8)
    ax1.set_yticks(np.arange(.985,1.003,.0025))
    ax1.set_yticklabels(['%.4f'%i for i in np.arange(.985,1.003,.0025)],
                        fontsize=7)
    ax1.set_xticks(np.arange(1355,1381,5))
    ax1.set_xticklabels(['%i'%i for i in np.arange(1355,1381,5)],
                        fontsize=7)

    fig.subplots_adjust(left=.18, right=.98, top=.98, bottom=.11)
    if label:
        plt.savefig('plots/TIC300_LC.png')
    if pltt:
        plt.show()
    plt.close('all')
