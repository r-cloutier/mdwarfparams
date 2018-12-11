from GAIAMdwarfs import *
import matplotlib.gridspec as gridspec
import matplotlib
from occurrencerateclass import *
import linear_lnlike as llnl
from TESS_search import *


matplotlib.rcParams.update({'ytick.labelsize': 9})
matplotlib.rcParams.update({'xtick.labelsize': 9})


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



def plot_1_planet_candidate_summary(self, SNRoffset=10, phaseylims=(.99,1.01),
                                    phasexlims=(-.1,.1), pltt=True,
                                    label=False):
    assert self.Ndet == 1

    # setup stuff
    cols = ['k','b','b']
    cols_linear = ['b','g','r']
    Z = self.params_guess[0,2]
    P, T0, aRs, rpRs, inc = self.params_optimized[0]
    rp = rvs.m2Rearth(rvs.Rsun2m(rpRs*self.Rs))
    b = rvs.impactparam_inc(P, self.Ms, self.Rs, inc)
    D = rvs.transit_width(P, self.Ms, self.Rs, rp, b)
    title_nums = self.Ndet, np.round(self.TESSmag,3), \
                 np.round(self.Teff), np.round(self.ehi_Teff), \
                 np.round(self.logg,3), np.round(self.ehi_logg,3), \
                 np.round(self.Rs,3), np.round(self.ehi_Rs,3)
    title = r'$\bf{'+'TIC'+'}$'+' '+r'$\bf{'+str(self.tic)+'}$'+' (%i planet candidate)\nT$_{mag}$ = %.3f,\tT$_{eff}$ = %i $\pm$ %i K,\t$\log{g}$ = %.3f $\pm$ %.3f,\tR$_s$ = %.3f $\pm$ %.3f R$_{\odot}$'%title_nums
    
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
        ax1.axvline(T0-t0+i*P, ymax=.03, lw=3, color=cols[2])
        ax1.axvline(T0-t0+i*P, ymin=.9, lw=3, color=cols[2])
    ax1.set_xlim(bjdlim)
    ax1.set_xticklabels('')
    ax1.set_ylabel('Normalized PDCSAP\nextracted flux', fontsize=9)
    ax1.set_title(title, fontsize=9, y=1)

    # plot the detrended light curve
    ax2 = fig.add_axes([.13,.7,.83,.12])
    ax2.plot(self.bjd-t0, self.fcorr, 'o', c=cols[0], ms=1, alpha=.5)
    for i in range(-50,50):
        ax2.axvline(T0-t0+i*P, ymax=.03, lw=3, color=cols[2])
        ax2.axvline(T0-t0+i*P, ymin=.9, lw=3, color=cols[2])
    ax2.set_xlim(bjdlim)
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
        ax3.axvline(T0-t0+i*P, ymax=.03, lw=3, color=cols[2])
        ax3.axvline(T0-t0+i*P, ymin=.95, lw=3, color=cols[2])
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
    ax5.plot(phase, self.fcorr, 'o', c=cols[0], ms=1, alpha=.2)
    ax5.errorbar(tb, fb, efb, fmt='o', c=cols[0], ms=4, elinewidth=1)
    ax5.plot(phase[s], self.fmodels[0][s], '-', c=cols[2], lw=2)
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
    e_rp = unp.std_devs(rvs.m2Rearth(rvs.Rsun2m(unp.uarray(rpRs,e_rpRs) * \
                                            unp.uarray(self.Rs,self.ehi_Rs))))
    ax5.text(1.48, .12, 'Derived parameters', fontsize=8, weight='semibold',
             transform=ax5.transAxes, horizontalalignment='center')
    ax5.text(1.05, 0,'r$_p$ = %.2f $\pm$ %.2f [R$_{\oplus}$]'%(np.round(rp,2),
                                                            np.round(e_rp,2)),
             fontsize=8, transform=ax5.transAxes)
    sma = rvs.semimajoraxis(unp.uarray(P,e_P),unp.uarray(self.Ms,self.ehi_Ms),0)
    Teq = unp.uarray(self.Teff,self.ehi_Teff) * \
          unp.sqrt(rvs.Rsun2m(unp.uarray(self.Rs,self.ehi_Rs)) / \
                   rvs.AU2m((2*sma))) * (1-.3)**(.25)
    Teq, e_Teq = unp.nominal_values(Teq), unp.std_devs(Teq)
    ax5.text(1.05,-.12,'T$_{eq}$ = %i $\pm$ %i [K] (Earth-like albedo: 0.3)'%(np.round(Teq),np.round(e_Teq)), fontsize=8, transform=ax5.transAxes)

    #fig.subplots_adjust(top=.93, right=.96, hspace=.2)
    if label:
        plt.savefig('plots/planetcandidatesummary_1_%i.png'%self.tic)
    if pltt:
        plt.show()
    plt.close('all')



def plot_linearsearch(self, dur_ind=0, pltt=True, label=''):

    cols = ['k','b','g']
    t0 = 2457000
    dur_ind = int(dur_ind)
    bjdlim = self.bjd.min()-t0-.3, self.bjd.max()-t0+.3

    fig = plt.figure(figsize=(3.7,2.7))
    ax = fig.add_subplot(111)    
    ax.plot(self.transit_times-t0, self.SNRs_linearsearch[:,dur_ind], '-', c=cols[dur_ind], lw=.4)
    ax.axhline(5, ls='--', color=cols[dur_ind], lw=1.5)
    ax.set_xlim(bjdlim)
    ax.set_xlabel('T$_0$ [BJD - 2,457,000]', fontsize=9)
    ax.set_ylabel('Linear search S/N =\n($\ln{\mathcal{L}}$ - median($\ln{\mathcal{L}}$)) / MAD($\ln{\mathcal{L}}$)',
                  fontsize=9)
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
    ax.plot(t0s[g1]-t0, snrs[g1], 'o', c=cols[1], ms=2)
    ax.plot(t0s[g2]-t0, snrs[g2], 'o', c=cols[2], ms=3)
    
    fig.subplots_adjust(bottom=.16, top=.92, right=.98, left=.16)
    if label:
        plt.savefig('plots/linearsearch_%i.png'%self.tic)
    if pltt:
        plt.show()
    plt.close('all')


def plot_GPdetrend(self, compute_thetaGP=True, pltt=True, label=False, ylim=(.996,1.006),
                   Npntsmin=5e2, Npntsmax=1e3, Nsig=3, medkernel=9):
    #assert np.all(self.quarters==0)
    thetaGPin_final, thetaGPout_final = self.thetaGPin[0], self.thetaGPout[0]
    cols = ['k','b','g','r']
    t0 = 2457000
    Nsects = np.unique(self.quarters).size
    #Nsects = 1 # TEMP
    
    # get thetaGPs
    if compute_thetaGP:
        thetain_comp2final, thetaout_comp2final, thetaGPsin, thetaGPsout = \
            _do_optimize_0_custom(self.bjd, self.f, self.ef, self.quarters)
        np.save('PipelineResults_TIC/thetain_comp2final_%i'%self.tic, thetain_comp2final)
        np.save('PipelineResults_TIC/thetaout_comp2final_%i'%self.tic, thetaout_comp2final)
        np.save('PipelineResults_TIC/thetaGPsin_%i'%self.tic, thetaGPsin)
        np.save('PipelineResults_TIC/thetaGPsout_%i'%self.tic, thetaGPsout)
    else:
        thetain_comp2final = np.load('PipelineResults_TIC/thetain_comp2final_%i.npy'%self.tic)
        thetaout_comp2final = np.load('PipelineResults_TIC/thetaout_comp2final_%i.npy'%self.tic)
        thetaGPsin = np.load('PipelineResults_TIC/thetaGPsin_%i.npy'%self.tic)
        thetaGPsout = np.load('PipelineResults_TIC/thetaGPsout_%i.npy'%self.tic)

    fig = plt.figure(figsize=(6.5,4))
    ax = fig.add_axes([.12,.3,.86,.66])
    ax2 = fig.add_axes([.12,.1,.86,.2])
    dts = np.zeros(Nsects)
    for j in range(Nsects):
        # get time binning from best-fit P_GP
        Npnts_per_timescale = 8.
        #timescale_to_resolve = np.exp(thetaGPout_final[3]) / Npnts_per_timescale
        timescale_to_resolve = np.exp(thetaGPsout[j,0,3]) / Npnts_per_timescale
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
        tbin, fbin, efbin = llnl.boxcar(self.bjd[s][g], medfilt(self.f[s][g],medkernel),
                                        self.ef[s][g], dt=dt)

        # plot best GP models
        ax.plot(self.bjd[s]-t0, self.f[s], 'o', c=cols[0], ms=1, alpha=.3)
        ax.plot(tbin-t0, fbin, 'o', ms=1.5, color=cols[2])
        # cover up pnt in second sector
        if j == 1:
            g3 = np.where(np.diff(tbin) > .5)[0][-1]
            #ax.plot(tbin[g3]-t0, fbin[g3], 'wo', ms=4)

        tmodel, mu, sig = _get_GP_model(thetaGPout_final, tbin, fbin, efbin,
                                        np.linspace(self.bjd[s].min(),self.bjd[s].max(),
                                                    1000))
        ax.fill_between(tmodel-t0, mu-sig, mu+sig, color=cols[1], alpha=.3)
        ax.plot(tmodel-t0, mu, '-', c=cols[1])
    
        # plot GP models from other iterations
        NGP = thetaGPsin.shape[1]
        assert NGP >= 9
        for i in range(9):
            tmodel, mu, sig = _get_GP_model(thetaGPsout[j,i], tbin, fbin, efbin, tmodel)
            if i == 0:
                ax.plot(tmodel-t0, mu, '-', c=cols[3], lw=.4, label='remaining GP models ')
            else:
                ax.plot(tmodel-t0, mu, '-', c=cols[3], lw=.4)

        # plot residuals
        trbin, rbin, erbin = llnl.boxcar(self.bjd[s][g], medfilt((self.f[s]-self.mu[s])[g]*1e6,medkernel),
                                        self.ef[s][g], dt=dt)
        ax2.plot(self.bjd[s]-t0, (self.f[s]-self.mu[s])*1e6, 'o', c=cols[0], ms=1, alpha=.3) 
        ax2.plot(trbin-t0, rbin, 'o', ms=1.5, color=cols[2])
        # cover up pnt in second sector
        if j == 1:
            g3 = np.where(np.diff(trbin) > .5)[0][-1]
            #ax2.plot(trbin[g3]-t0, rbin[g3], 'wo', ms=4)

        
    if Nsects > 1:
        x = self.bjd.mean()-t0+.2
        ax.axvline(x, ls='--', color='k', lw=2)
        ax2.axvline(x, ls='--', color='k', lw=2)
        ax.text(x-.7, 1.0052, 'sector 1', horizontalalignment='right', fontsize=6,
                weight='semibold')
        ax.text(x+.7, 1.0052, 'sector 2', horizontalalignment='left', fontsize=6,
                weight='semibold')
        
    ax2.set_xlabel('Time [BJD - 2,457,000]', fontsize=9)
    ax.set_ylabel('Normalized flux', fontsize=9)
    ax2.set_ylabel('De-trended flux\n[ppm]', fontsize=9, labelpad=.1)
    ax.set_ylim(ylim)
    bjdlim = self.bjd.min()-t0-.3, self.bjd.max()-t0+.3
    ax.set_xlim(bjdlim)
    ax2.set_xlim(bjdlim)

    # custom legend
    ax.plot(.05, .95, 'o', c=cols[0], ms=1, alpha=.8, transform=ax.transAxes)
    ax.text(.08, .95, 'raw photometry', fontsize=6, transform=ax.transAxes,
            verticalalignment='center')
    ax.plot(.05, .88, 'o', c=cols[2], ms=1.5, transform=ax.transAxes)
    endings = [', ',')']
    dt_labels = ''.join(['dt$_%i$ = %.2f hrs%s'%(i+1,dts[i],endings[i])
                         for i in range(Nsects)])
    ax.text(.08, .88, 'binned photometry\n(%s'%dt_labels, fontsize=6,
            transform=ax.transAxes, verticalalignment='center')
    ax.fill_between([.64,.67], [.94,.94], [.96,.96], color=cols[1], alpha=.3,
                    transform=ax.transAxes)
    ax.plot([.64,.67], [.95,.95], '-', color=cols[1], lw=1, transform=ax.transAxes)
    ax.text(.69, .95, 'max($\ln{\mathcal{L}}$) GP models mean and\nstandard deviation', fontsize=6,
            transform=ax.transAxes, verticalalignment='center')
    ax.plot([.64,.67], [.87,.87], '-', color=cols[3], lw=.7, transform=ax.transAxes)
    ax.text(.69, .87, 'unused GP models from remaining\niterations', fontsize=6,
            transform=ax.transAxes, verticalalignment='center')
    
    fig.subplots_adjust(bottom=.1, top=.95, right=.98, left=.12)
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


def plot_stellar_corner(dTIC, dTIC_sect1, dTIC_sect2, pltt=True, label=False):
    assert dTIC.shape == (93090,89)  #np.genfromtxt('input_data/TESStargets/TICv7_Mdwarfsv1.csv', delimiter=',', skip_header=5)
    assert dTIC_sect1.shape == (905, 39) #np.loadtxt('input_data/TESStargets/TESSMdwarfs_sector1_v2.csv', delimiter=',')
    assert dTIC_sect2.shape == (1062,39) #np.loadtxt('input_data/TESStargets/TESSMdwarfs_sector2_v2.csv', delimiter=',')

    cols = ['b','k']
    
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
                ax.hist(list1[i], bins=bins[i], color=cols[0], normed=1, alpha=.3)
                y2,_,_=ax.hist(list2[i], bins=bins[i], color=cols[1], normed=1,
                               histtype='step',lw=1.5)
                ax.hist(list2[i], bins=bins[i], color=cols[1], normed=1, alpha=.5)
                ax.set_xlim(lims[i])
                if j < 3: ax.set_xticklabels('')
                if j == 3:
                    ax.set_xlabel(labels[j], fontsize=7, labelpad=1)
                    ax.set_xticks(tickvals[j])
                    ax.set_xticklabels(ticklabels[j], fontsize=5)
                ax.set_yticklabels('')
                ax.plot(np.repeat(np.nanmedian(list1[i]),2),[0,np.max([y1,y2])], '--',
                        lw=.8,c=cols[0])
                ax.plot(np.repeat(np.nanmedian(list2[i]),2),[0,np.max([y1,y2])], '--',
                        lw=.8, c=cols[1])
                ax.text(.49, 1.05, titles1[i]%np.nanmedian(list1[i]), fontsize=7,
                        transform=ax.transAxes, horizontalalignment='right',
                        weight='semibold', color=cols[0])
                ax.text(.51, 1.05, titles2[i]%np.nanmedian(list2[i]), fontsize=7,
                        transform=ax.transAxes, horizontalalignment='left',
                        weight='semibold', color=cols[1])
                
            # 2d histogram
            elif i > j:
                ax.plot(list1[j], list1[i], '.', c=cols[0], ms=1, alpha=.1)
                ax.plot(list2[j], list2[i], '.', c=cols[1], ms=1.5, alpha=.9)
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
    assert hasattr(self, self.isTESSalert)
    g = self.disposition_human >= 0 
    g0 = self.disposition_human[g] == 0
    g1 = self.disposition_human[g] == 1
    
    cols = ['b','k']

    # plot planet candidates that passed human vetting (and maybe vespa?)
    fig = plt.figure(figsize=(6.5,3))

    # plot P
    ax1 = fig.add_subplot(131)
    ax1.errorbar(self.Ps[g0], self.rps[g0], self.ehi_rps[g0], fmt='o', ms=2,
                 color=cols[0], elinewidth=1, capsize=0)
    ax1.errorbar(self.Ps[g1], self.rps[g1], self.ehi_rps[g1], fmt='o', ms=2,
                 color=cols[1], elinewidth=1, capsize=0)
    ax1.set_xscale('log'), ax.set_yscale('log')
