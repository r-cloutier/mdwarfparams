from GAIAMdwarfs import *
import matplotlib.gridspec as gridspec
import matplotlib
from occurrencerateclass import *


matplotlib.rcParams.update({'ytick.labelsize': 11})
matplotlib.rcParams.update({'xtick.labelsize': 11})


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
