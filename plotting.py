from GAIAMdwarfs import *
import matplotlib.gridspec as gridspec
import matplotlib
from occurrencerateclass import *
import linear_lnlike as llnl


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



def plot_1_planet_candidate_summary(self, SNRoffset=10, pltt=True, label=False):
    assert self.Ndet == 1

    # setup stuff
    cols = ['k','b','b','g']
    cols_linear = ['b','g','r']
    Z = self.params_guess[0,2]
    P, T0, aRs, rpRs, inc = self.params_optimized[0]
    rp = rvs.m2Rearth(rvs.Rsun2m(rpRs*self.Rs))
    b = rvs.impactparam_inc(P, self.Ms, self.Rs, inc)
    D = rvs.transit_width(P, self.Ms, self.Rs, rp, b)
    title_nums = self.tic, np.round(self.TESSmag,3), np.round(self.Teff), \
                 np.round(self.ehi_Teff), np.round(self.logg,3), \
                 np.round(self.ehi_logg,3), np.round(self.Rs,3), \
                 np.round(self.ehi_Rs,3)
    title = 'TIC %i\nT$_{mag}$ = %.3f,\tT$_{eff}$ = %i $\pm$ %i K,\t$\log{g}$ = %.3f $\pm$ %.3f,\tR$_s$ = %.3f $\pm$ %.3f R$_{\odot}$'%title_nums
    
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
        ax3.text(self.bjd.min()-t0+1, -i*SNRoffset+6,
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
    ax4.set_ylim((np.round(1-2.5*Z,3), np.round(1+2.5*Z,3)))
    ax4.set_ylabel('Normalized\nde-trended flux', fontsize=9)

    # zoom-in on transit
    ax5 = fig.add_axes([.13,.05,.45,.15])
    ax5.plot(phase, self.fcorr, 'o', c=cols[0], ms=1, alpha=.5)
    ax5.errorbar(tb, fb, efb, fmt='o', c=cols[2], ms=4, elinewidth=1)
    ax5.plot(phase[s], self.fmodels[0][s], '-', c=cols[3], lw=2)
    ax5.set_xlim((-2.9*D/P,2.9*D/P))
    ax5.set_ylim((np.round(1-2.5*Z,3), np.round(1+2.5*Z,3)))
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
    ax4.text(1.05, .34, '$\ln$a = %.3f'%(np.round(self.thetaGPout[0,0],3)),
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
    ax5.text(1.05,-.12,'T$_{eq}$ = %i $\pm$ %i [K] (Earth-like albedo; 0.3)'%(np.round(Teq),np.round(e_Teq)), fontsize=8, transform=ax5.transAxes)

    #fig.subplots_adjust(top=.93, right=.96, hspace=.2)
    if label:
        plt.savefig('plots/planetcandidatesummary_1_%i.png'%self.tic)
    if pltt:
        plt.show()
    plt.close('all')
