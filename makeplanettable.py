from occurrencerateclass import *
from followup import *


def planetary_table(self, savetex=False):

    g = self.disposition_human >= 0
    Np = g.sum()
    assert Np == 16
    g = np.where(g)[0][np.argsort(self.tics_candidates[g])]
    disp_dict = {0:'pPC', 1:'PC', 2:'ST', 2.5:'pST'}
    t0 = 2457000
    
    # make latex table of planet properties
    Ncol_notTIC = 9
    hdr = '\\begin{longrotatetable}\n'
    hdr += '\\begin{deluxetable*}{c%s}\n'%(Ncol_notTIC*'c')
    hdr += '\\tabletypesize{\small}\n'
    hdr += '\\tablecaption{Planetary parameters for our 16 vetted candidates\label{table:planets}}\n'

    # set planet column headings
    headings = ['OI','TOI','$P$','$T_0$','$a/R_s$','$r_p/R_s$','$i$',
                '$r_p$','$S$','Disposition\\tablenotemark{a}']
    units = ['','','[days]','[BJD-2,457,000]','','','[deg]','[R$_{\\oplus}$]',
             '[S$_{\\oplus}$]','']
    assert len(headings) == len(units)
    assert len(headings) == Ncol_notTIC+1
    head1 = ' & '.join(headings)
    head2 = ' & '.join(units)
    hdr += '\\tablehead{%s \\\\ \n %s}\n'%(head1, head2)
    hdr += '\\startdata\n'

    # add parameter values
    for i in range(Np):

        # get values depding on if its a single transit or not
        if self.disposition_human[g][i] in [0,1]:
            Pentry = '%.4f $\\pm$ %.4f'%(self.Ps[g][i], self.e_Ps[g][i])
	    aRsentry = '%.2f$^{+%.2f}_{-%.2f}$'%(self.aRss[g][i],
                                                 self.ehi_aRss[g][i],
                                                 self.elo_aRss[g][i]) 
	    rpRsentry = '%.3f$^{+%.3f}_{-%.3f}$'%(self.rpRss[g][i],
                                                  self.ehi_rpRss[g][i],
                                                  self.elo_rpRss[g][i])
            incentry = '%.2f$^{+%.2f}_{-%.2f}$'%(self.incs[g][i],
                                                 self.ehi_incs[g][i],
                                                 self.elo_incs[g][i])
            rpentry = '%.2f$^{+%.2f}_{-%.2f}$'%(self.rps[g][i],
                                                self.ehi_rps[g][i],
                                                self.elo_rps[g][i])
            Fentry = '%.1f$\\pm$ %.1f'%(self.Fs[g][i], self.ehi_Fs[g][i])
            
        elif self.disposition_human[g][i] in [2,2.5]:
            Pentry = '%.2f$^{+%.2f}_{-%.2f}$'%(self.Ps_singletransit[g][i],
                                               self.ehi_Ps_singletransit[g][i],
                                               self.elo_Ps_singletransit[g][i])
            aRsentry = '%.2f$^{+%.2f}_{-%.2f}$'%(self.aRs_singletransit[g][i],
                                            self.ehi_aRs_singletransit[g][i],
                                            self.elo_aRs_singletransit[g][i])
            rpRsentry = '%.2f$^{+%.2f}_{-%.2f}$'%(self.rpRs_singletransit[g][i],
                                            self.ehi_rpRs_singletransit[g][i],
                                            self.elo_rpRs_singletransit[g][i])
            incentry = '%.2f$^{+%.2f}_{-%.2f}$'%(self.inc_singletransit[g][i],
                                            self.ehi_inc_singletransit[g][i],
                                            self.elo_inc_singletransit[g][i])
            rpentry = '%.2f$^{+%.2f}_{-%.2f}$'%(self.rps_singletransit[g][i],
                                            self.ehi_rps_singletransit[g][i],
                                            self.elo_rps_singletransit[g][i])
            Fentry = '%.1f$^{+%.1f}_{-%.1f}$'%(self.Fs_singletransit[g][i],
                                               self.ehi_Fs_singletransit[g][i],
                                               self.elo_Fs_singletransit[g][i])

        TOIentry = '%.2f'%self.tois[g][i] if self.tois[g][i] > 0 else '-'
            
        p = self.tics_candidates[g][i], TOIentry, Pentry, self.T0s[g][i]-t0, self.e_T0s[g][i], aRsentry, rpRsentry, incentry, rpentry, Fentry, disp_dict[self.disposition_human[g][i]]
        hdr += '%.2f & %s & %s & %.4f $\\pm$ %.4f & %s & %s & %s & %s & %s & %s \\\\ \n'%(p)

    # add footer
    hdr += '\\enddata\n'
    hdr += "\\tablenotetext{a}{Possible dispositions of objects of interest are a planet candidate (PC), a putative planet candidate (pPC), a single transit event (ST), or a putative single transit event (pST). The putative dispositions have FPP $\in [0.1,0.9]$ whereas the remaining candidates have FPP $<0.1$.}\n"
    hdr += '\\end{deluxetable*}\n\\end{longrotatetable}'

    # save tex table
    if savetex:
        f = open('TESSpaper/planetparams.tex', 'w')
        f.write(hdr)
        f.close()

        
    return hdr


def planetary_table_txt(self, savetxt=False):

    g = self.disposition_human >= 0
    Np = g.sum()
    assert Np == 16
    g = np.where(g)[0][np.argsort(self.tics_candidates[g])]
    disp_dict = {0:'pPC', 1:'PC', 2:'ST', 2.5:'pST'}
    t0 = 2457000
    
    # make latex table of planet properties
    Ncol_notTIC = 9
    hdr = '#OI,TOI,P,ehi_P,elo_P,T0,ehi_T0,elo_T0,aRs,ehi_aRs,elo_aRs,rpRs,ehi_rpRs,elo_rpRs,inc,ehi_inc,elo_inc,rp,ehi_rp,elo_rp,S,ehi_s,elo_S,Disposition\n'
    
    # add parameter values
    for i in range(Np):

        # get values depding on if its a single transit or not
        if self.disposition_human[g][i] in [0,1]:
            Pentry = '%.4f,%.4f,%.4f'%(self.Ps[g][i], self.e_Ps[g][i], self.e_Ps[g][i])
	    aRsentry = '%.2f,%.2f,%.2f'%(self.aRss[g][i],
                                         self.ehi_aRss[g][i],
                                         self.elo_aRss[g][i]) 
	    rpRsentry = '%.3f,%.3f,%.3f'%(self.rpRss[g][i],
                                          self.ehi_rpRss[g][i],
                                          self.elo_rpRss[g][i])
            incentry = '%.2f,%.2f,%.2f'%(self.incs[g][i],
                                         self.ehi_incs[g][i],
                                         self.elo_incs[g][i])
            rpentry = '%.2f,%.2f,%.2f'%(self.rps[g][i],
                                        self.ehi_rps[g][i],
                                        self.elo_rps[g][i])
            Fentry = '%.1f,%.1f,%.1f'%(self.Fs[g][i], self.ehi_Fs[g][i], self.ehi_Fs[g][i])
            
        elif self.disposition_human[g][i] in [2,2.5]:
            Pentry = '%.2f,%.2f,%.2f'%(self.Ps_singletransit[g][i],
                                       self.ehi_Ps_singletransit[g][i],
                                       self.elo_Ps_singletransit[g][i])
            aRsentry = '%.2f,%.2f,%.2f'%(self.aRs_singletransit[g][i],
                                         self.ehi_aRs_singletransit[g][i],
                                         self.elo_aRs_singletransit[g][i])
            rpRsentry = '%.2f,%.2f,%.2f'%(self.rpRs_singletransit[g][i],
                                          self.ehi_rpRs_singletransit[g][i],
                                          self.elo_rpRs_singletransit[g][i])
            incentry = '%.2f,%.2f,%.2f'%(self.inc_singletransit[g][i],
                                         self.ehi_inc_singletransit[g][i],
                                         self.elo_inc_singletransit[g][i])
            rpentry = '%.2f,%.2f,%.2f'%(self.rps_singletransit[g][i],
                                        self.ehi_rps_singletransit[g][i],
                                        self.elo_rps_singletransit[g][i])
            Fentry = '%.1f,%.1f,%.1f'%(self.Fs_singletransit[g][i],
                                       self.ehi_Fs_singletransit[g][i],
                                       self.elo_Fs_singletransit[g][i])

        TOIentry = '%.2f'%self.tois[g][i] if self.tois[g][i] > 0 else '-'
            
        p = self.tics_candidates[g][i], TOIentry, Pentry, self.T0s[g][i]-t0, self.e_T0s[g][i], self.e_T0s[g][i], aRsentry, rpRsentry, incentry, rpentry, Fentry, disp_dict[self.disposition_human[g][i]]
        hdr += '%.2f,%s,%s,%.4f,%.4f,%.4f,%s,%s,%s,%s,%s,%s\n'%(p)

    # save txt table
    if savetxt:
        f = open('TESSpaper/planetparams.csv', 'w')
        f.write(hdr)
        f.close()

        
    return hdr



def planetary_followup_table(self, savetex=False):

    g = self.disposition_human >= 0
    Np = g.sum()
    assert Np == 16
    g = np.where(g)[0][np.argsort(self.tics_candidates[g])]
    
    # make latex table of planet properties
    Ncol_notTIC = 10
    hdr = '\\begin{deluxetable*}{c%s}\n'%(Ncol_notTIC*'c')
    hdr += '\\tabletypesize{\small}\n'
    hdr += '\\tablecaption{Metric values indicating the feasibility of a variety of follow-up programs for our 16 vetted candidates\label{table:followup}}\n'

    # set planet column headings
    headings = ['OI','TOI','$J$','$P$','$r_p$','$m_p$\\tablenotemark{a}','$K$','$\\Omega$\\tablenotemark{b}',
                '$T_{\\text{eq}}$\\tablenotemark{c}',
                'TSM\\tablenotemark{d}','ESM\\tablenotemark{e}']
    units = ['','','','[days]','[R$_{\oplus}$]','[M$_{\\oplus}$]','[m/s]','','[K]','','']
    assert len(headings) == len(units)
    assert len(headings) == Ncol_notTIC+1
    head1 = ' & '.join(headings)
    head2 = ' & '.join(units)
    hdr += '\\tablehead{%s \\\\ \n %s}\n'%(head1, head2)
    hdr += '\\startdata\n'

    # add parameter values
    for i in range(Np):

        # compute followup parameters
        P = self.Ps_singletransit[g][i] if self.disposition_human[g][i] >= 2 \
            else self.Ps[g][i]
        rp = self.rps_singletransit[g][i] if self.disposition_human[g][i] >= 2 \
             else self.rps[g][i]
        Jmag = self.Jmags[g][i]
        Omega = rp / P**(1./3)
        exceedsOmega = does_exceed_Omega_cutoff(Jmag, Omega)
        Omegaentry = '\\textbf{%.2f}'%Omega if exceedsOmega else '%.2f'%Omega

        mp,K,Teq,SF,TSM = estimate_transmission_metric(P, rp, self.Jmags[g][i],
                                                       self.Mss[g][i],
                                                       self.Rss[g][i],
                                                       self.Teffs[g][i])
        exceedsTSM = does_exceed_TSM_cutoff(rp, TSM)
        TSMentry = '\\textbf{%.1f}'%TSM if exceedsTSM else '%.1f'%TSM
        mp,K,Teq,Br,ESM = estimate_eclipse_metric(P, rp, self.Kmags[g][i],
                                                  self.Mss[g][i],
                                                  self.Rss[g][i],
                                                  self.Teffs[g][i])
        exceedsESM = does_exceed_ESM_cutoff(ESM)
        ESMentry = '\\textbf{%.1f}'%ESM if exceedsESM else '%.1f'%ESM
        TOIentry = '%.2f'%self.tois[g][i] if self.tois[g][i] > 0 else '-'
        
        p = self.tics_candidates[g][i], TOIentry, Jmag, P, rp, mp, K, Omegaentry, Teq, \
            TSMentry, ESMentry
        hdr += '%.2f & %s & %.2f & %.2f & %.2f & %.2f & %.2f & %s & %i & %s & %s \\\\ \n'%(p)

    # add footer
    hdr += '\\enddata\n'
    hdr += '\\tablecomments{Bolded values are indicative of candidates that exceed threshold values of that parameter (see \\citealt{cloutier18b} for $\\Omega$ and \\citealt{kempton18} for the TSM and ESM) and should be strongly considered for rapid confirmation and follow-up.}\n'
    hdr += '\\tablenotetext{a}{Planet masses are estimated from the planet radius using the deterministic version of the mass-radius relation from \cite{chen17}.}'
    hdr += "\\tablenotetext{b}{$\\Omega$ is a diagnostic metric that is indicative of the observing time required to characterize a planet's RV mass \\citep{cloutier18b}. $\\Omega= r_p/P^{1/3}$ where $r_p$ is given in Earth radii and $P$ in days.}\n"
    hdr += "\\tablenotetext{c}{Planetary equilibrium temperature is calculated assuming zero albedo and full heat redistribution via $T_{\\text{eq}} = T_{\\text{eff}} \\sqrt{R_s/2a}$.}\n"
    hdr += "\\tablenotetext{d}{The transmission spectroscopy metric from \\citep{kempton18}. See Sect.~\\ref{sect:atmospheres} for the definition.}\n"
    hdr += "\\tablenotetext{e}{The emission spectroscopy metric from \\citep{kempton18}. See Sect.~\\ref{sect:atmospheres} for the definition.}\n"
    hdr += '\\end{deluxetable*}'

    # save tex table
    if savetex:
        f = open('TESSpaper/planetfollowup.tex', 'w')
        f.write(hdr)
        f.close()
        
    return hdr



def does_exceed_Omega_cutoff(Jmag, Omega):
    return Omega >= .14*Jmag - .35


def does_exceed_TSM_cutoff(rp, TSM):
    if rp < 1.5:
        cutoff = 12
    elif 1.5 < rp < 2.75:
        cutoff = 92
    elif 2.75 < rp < 4:
        cutoff = 84
    elif 4 < rp < 10:
        cutoff = 96

    return TSM >= cutoff


def does_exceed_ESM_cutoff(ESM):
    #_,_,_,_,ESM_gj1132 = estimate_eclipse_metric(1.6289, 1.15, 8.322,
    #                                             .181, .201, 3270)
    ESM_gj1132 = 7.5
    print 'gj1132 ESM = %.5f'%ESM_gj1132
    return ESM >= ESM_gj1132
