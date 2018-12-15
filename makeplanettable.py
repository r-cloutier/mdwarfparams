from occurrencerateclass import *


def planetary_table(self, savetex=False):

    g = self.disposition_human >= 0
    Np = g.sum()
    disp_dict = {0:'pPC', 1:'PC', 2:'ST'}
    t0 = 2457000
    
    # make latex table of planet properties
    Ncol_notTIC = 8
    hdr = '\\begin{deluxetable*}{%s}\n'%(Ncol*'c')
    hdr += '\\tabletypesize{\small}\n'
    hdr += '\\tablecaption{Planet parameters of vetted objects-of-interest\label{table:starphot}}\n'

    # set planet column headings
    headings = ['TIC','$P$','$T_0$','$a/R_s$','$r_p/R_s$','$i$','$r_p$','$S$',
                'Disposition']
    units = ['','[days]','[BJD-2,457,000]','','','[deg]','[R$_{\\oplus}$]',
             '[S$_{\\oplus}$]','']
    assert len(headings) == len(units)
    assert len(headings) == Ncol_notTIC+1
    head1 = ' & '.join(headings)
    head2 = ' & '.join(units)
    hdr += '\\tablehead{%s \\\\ \n %s}\n'%(head1, head2)
    hdr += '\\startdata\n'

    # add parameter values
    for i in range(Np):
        d = loadpickle('PipelineResults_TIC/TIC_%i/LC_-00099_8d4'%(self.tics[g][i]))
        ind = np.isclose(d.params_optimized[:,0], self.Ps[g][i], rtol=.02)
        P, e_P = self.Ps[g][i], self.e_Ps[g][i]
        T0, e_T0 = d.params_optimized[ind,1]-t0, \
                   get_1sigma(d.params_results[ind,:,1].reshape(3))
        #aRs = d.params_optimized[ind,2]
        #p84 = d.params_results[ind,0,2] + d.params_results[ind,1,2]
        #p16 = d.params_results[ind,0,2] - d.params_results[ind,2,2]
        aRs, ehi_aRS, elo_aRs = np.zeros(3)
        samp_rpRs = d.params_optimized[ind,3], get_1sigma(params_res[:,3])
        rp, ehi_rp, wlo_rp
        
        hdr += '%i & %.4f $\\pm$ %.4f & %.4f $\\pm$ %.4f & %.3f$^{+%.3f}_{-%.3f}$ & %.3f$^{+%.3f}_{-%.3f} & %.2f $\\pm$ %.2f & %.3f $\\pm$ %.3f & %.1f $\\pm$ %.1f & %s \\\\ \n'%(d.tic, P, e_P, T0, e_T0,  aRa, ehi_aRs, elo_aRs, rpRs, ehi_rpRs, elo_rpRs, inc, ehi_inc, elo_inc, rp, ehi_rp, elo_rp, F, ehi_F, elo_F)
