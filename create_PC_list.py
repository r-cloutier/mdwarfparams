from occurrencerateclass import *
import warnings, datetime
from uncertainties import unumpy as unp


def get_PCs(self):
    '''Given a folder with an OccurrenceRateclass object, get all the planet 
    candidates from orion and with human dispositions, so that they can be 
    added to the list of planet candidates.'''

    # get the results object
    #self = loadpickle('%s/TIC_results_0_10000_det'%folder)
    assert hasattr(self, 'tics')
    assert hasattr(self, 'disposition_human')

    # isolate planet candidates, putative PCs, and single transits
    g01 = np.in1d(self.disposition_human, [0,1])  # PCs and pPCs
    g2 = np.in1d(self.disposition_human , [2,2.5])  # STs and pSTs
    g = (g01 + g2).astype(bool)
    assert g.sum() == g01.sum() + g2.sum()

    # get each TIC's sectors
    sector_dict = get_sectors(self.tics[g])
    
    # get columns of interest
    # TIC,sectors,ra,dec,Tmag,Jmag,Hmag,Kmag,dist,Ms,Rs,Teff,P,T0,rp/Rs,rp,disp
    outstr = '#TIC,sectors,ra,dec,Tmag,Jmag,Hmag,Kmag,dist,Ms,Rs,Teff,P,T0,rp/Rs,rp,disposition\n'
    for i in np.where(g)[0]:
        outstr += '%i;%s;%.5f;%.5f;'%(self.tics[i], sector_dict[self.tics[i]],
                                      self.ras[i], self.decs[i])
        outstr += '%.3f;%.3f;%.3f;%.3f;'%(self.TESSmags[i], self.Jmags[i],
                                          self.Hmags[i], self.Kmags[i])
        outstr += '%.3f;%.3f;%.3f;%i;'%(self.dists[i], self.Mss[i],
                                        self.Rss[i], self.Teffs[i])

        if self.disposition_human[i] < 2:
            p, rpRs, rp = self.Ps[i], self.rpRss[i], self.rps[i]
        elif self.disposition_human[i] >= 2:
            p = self.Ps_singletransit[i]
            rpRs = self.rpRs_singletransit[i]
            rp = self.rps_singletransit[i]
        else:
            raise ValueError('disposition %.1f is not valid.'%self.disposition_human[i])
            
        outstr += '%.6f;%.6f;%.5f;%.2f;%s\n'%(p, self.T0s[i], rpRs, rp,
                                              get_disposition(self.disposition_human[i]))

    return outstr



def get_CTOIs(self, fname_index):
    '''Given a folder with an OccurrenceRateclass object, get all the planet 
    candidates from orion with human dispositions and create a list of CTOIs using the 
    CTOI format from ExoFOP 
    (https://exofop.ipac.caltech.edu/tess/templates/params_planet_YYYYMMDD_001.txt)'''

    assert hasattr(self, 'tics')
    assert hasattr(self, 'disposition_human')
    assert np.all(np.in1d(fname_index, np.arange(1,1e3,dtype=int)))

    # isolate planet candidates, putative PCs, and single transits
    g1 = np.in1d(self.disposition_human, [1])  # PCs and pPCs
    g2 = np.in1d(self.disposition_human , [2])  # STs and pSTs
    g = (g1 + g2).astype(bool)
    assert g.sum() == g1.sum() + g2.sum()
   
    # get string
    date = datetime.datetime.now()
    date_str = '%.4d%.2d%.2d'%(date.year, date.month, date.day)
 
    # get columns of interest including uncertainties (NaNs where appropriate)
    # TIC,flag,disp,P,T0,Z,D,inc,b,rp/Rs,a/Rs,rp,mp,Teq,S,rho_star,sma,ecc,omega,tau_peri,K_RV,tag,group,prop_P,notes
    # for definitions see ORION_PCs/ctoi_header.txt 
    f = open('ORION_PCs/ctoi_header.txt', 'r')
    outstr = f.read()
    f.close()
    for i in np.where(g)[0]:

        # check if potential CTOIs are already a TOI or a CTOI
        if (not is_TOI(self.tics[i])) and (not is_CTOI(self.tics[i])):

        
            # get multiple transit parameters
            if self.disposition_human[i] < 2:
                P, eP = self.Ps[i], self.e_Ps[i]
                inc, einc = self.incs[i], np.mean([self.ehi_incs[i],self.elo_incs[i]])
                rpRs, erpRs = self.rpRss[i], np.mean([self.ehi_rpRss[i],self.elo_rpRss[i]])
                aRs, eaRs = self.aRss[i], np.mean([self.ehi_aRss[i],self.elo_aRss[i]])
                rp, erp = self.rps[i], np.mean([self.ehi_rps[i],self.elo_rps[i]])
                sma, esma = self.smas[i], np.mean([self.ehi_smas[i],self.elo_smas[i]])
                notes = 'new CTOI from ORION (https://arxiv.org/abs/1812.08145)'

            # get single transit parameters
            elif self.disposition_human[i] >= 2:
                P, eP = np.nan, np.nan#self.Ps_singletransit[i], np.mean([self.ehi_Ps_singletransit[i],self.elo_Ps_singletransit[i]])
                inc, einc = self.inc_singletransit[i], np.mean([self.ehi_inc_singletransit[i],self.elo_inc_singletransit[i]])
                rpRs, erpRs = self.rpRs_singletransit[i], np.mean([self.ehi_rpRs_singletransit[i],self.elo_rpRs_singletransit[i]])
                aRs, eaRs = self.aRs_singletransit[i], np.mean([self.ehi_aRs_singletransit[i],self.elo_aRs_singletransit[i]])
                rp, erp = self.rps_singletransit[i], np.mean([self.ehi_rps_singletransit[i],self.elo_rps_singletransit[i]])
                #sampMs = np.random.randn(1000)*self.ehi_Mss[i] + self.Mss[i]
                #_,_,sampP = get_samples_from_percentiles(P, self.ehi_Ps_singletransit[i], self.elo_Ps_singletransit[i])
                #sampsma = rvs.semimajoraxis(sampP, sampMs, 0)
                #v = np.percentile(sampsma, (16,50,84))
                sma, esma = np.nan, np.nan#rvs.semimajoraxis(P, self.Mss[i], 0), np.mean([v[2]-v[1], v[1]-v[0]])
                notes = 'new single transit CTOI from ORION (https://arxiv.org/abs/1812.08145)'

            else:
                raise ValueError('disposition %.1f is not valid.'%self.disposition_human[i])

        
            # get other planet params
            T0, eT0 = self.T0s[i], self.e_T0s[i]
            Z, eZ = self.depths[i], unp.std_devs(unp.uarray(rpRs,erpRs)**2)
            ub = rvs.impactparam_inc_aRs(unp.uarray(aRs,eaRs), unp.uarray(inc,einc))
            uD = rvs.transit_width_aRs(unp.uarray(P,eP), unp.uarray(aRs,eaRs), unp.uarray(Z,eZ), ub)
            uTeq = unp.uarray(self.Teffs[i],self.ehi_Teffs[i]) * unp.sqrt(rvs.Rsun2m(unp.uarray(self.Rss[i],self.ehi_Rss[i])) \
                                                                          / rvs.AU2m(2*unp.uarray(sma,esma)))
            uS = unp.uarray(self.Rss[i],self.ehi_Rss[i])**2 * (unp.uarray(self.Teffs[i],self.ehi_Teffs[i])/5777)**4 * (1./unp.uarray(sma,esma))**2
            D, eD = unp.nominal_values(uD), unp.std_devs(uD)
            b, eb = unp.nominal_values(ub), unp.std_devs(ub)
            Teq, eTeq = unp.nominal_values(uTeq), unp.std_devs(uTeq)
            S, eS = unp.nominal_values(uS), unp.std_devs(uS)

            # add ancillary stuff
            flag, disp = 'newctoi', 'PC'
            tag = '%s_cloutier_orion_%.5d'%(date_str, fname_index)
            group, prop_P = '', 0
        
            # add planet parameters
            outstr += 'TIC%i.01|%s|%s|'%(self.tics[i], flag, disp)
            outstr += '%.6f|%.6f|%.5f|%.5f|'%(P, eP, T0, eT0)
            outstr += '%.1f|%.1f|%.3f|%.3f|%.2f|%.2f|'%(Z, eZ, D, eD, inc, einc)
            outstr += '%.3f|%.3f|%.4f|%.4f|%.1f|%.1f|'%(b, eb, rpRs, erpRs, aRs, eaRs)
            outstr += '%.2f|%.2f|%.2f|%.2f|%.1f|%.1f|'%(rp, erp, np.nan, np.nan, Teq, eTeq)
            outstr += '%.1f|%.1f|%.2f|%.2f|%.4f|%.4f|'%(S, eS, np.nan, np.nan, sma, esma)
            outstr += '%.2f|%.2f|%.2f|%.2f|%.2f|%.2f|'%(0, 0, np.nan, np.nan, np.nan, np.nan)
            outstr += '%.2f|%.2f|'%(np.nan, np.nan)
            outstr += '%s|%s|%i|%s\n'%(tag, group, prop_P, notes)

    # replace NaNs
    outstr = outstr.replace('nan','')
    
    # remind me to check stuff:
    # TIC name is not a duplicate
    # TO is correct for single transit events
    warnings.warn('\nThe parameters here (e.g. T0 for single transits) should be confirmed by-eye.') 

    # save file
    fname_out = 'ORION_PCs/params_planet_%s_%.3d.txt'%(date_str, fname_index)
    f = open(fname_out, 'w')
    f.write(outstr)
    f.close()
    
    return outstr



def is_TOI(tic):
    '''check if the input TIC is already a TOI. The up-to-date list of TOIs can be downloaded from
    list https://exofop.ipac.caltech.edu/tess/index.php'''
    # get published list of all TOIs
    toi_tics = np.loadtxt('ORION_PCs/exofop_tess_tois.csv', delimiter=',', skiprows=3, usecols=(0), dtype='|S50')
    # check if the input tic is in the list of tois
    return '"%s"'%tic in toi_tics 


def is_CTOI(tic):
    '''check if the input TIC is already a CTOI. The up-to-date list of CTOIs can be downloaded from
    list https://exofop.ipac.caltech.edu/tess/index.php'''
    # get published list of all CTOIs
    ctoi_tics = np.loadtxt('ORION_PCs/exofop_tess_ctois.csv', delimiter=',', skiprows=3, usecols=(0), dtype='|S50')
    # check if the input tic is in the list of ctois
    return '"%s"'%tic in ctoi_tics 

    

def get_sectors(tics):
    '''Given a list of TICs, read-in the individual TESS sector targets lists 
    and identify which sector the TIC was observed in.'''
    # get mapping between tics and sectors
    fs = np.array(glob.glob('input_data/TESStargets/all_targets_S*'))
    tics_tot, sectors_tot = np.zeros(0), np.zeros(0)
    for f in fs:
        if '.txt' in f:
            tics_thissector = np.loadtxt(f)[:,0]
        elif '.csv' in f:
            tics_thissector = np.loadtxt(f, delimiter=',', skiprows=6)[:,0]
        else:
            raise ValueError('Cannot read this file: %s'%f)

        # map tics to their sector
        tics_tot = np.append(tics_tot, tics_thissector)
        sect = int(f.split('_')[-2].split('S')[1])
        sectors_tot = np.append(sectors_tot,
                                np.repeat(sect, tics_thissector.size))

    # save each tic's sector(s)
    assert tics_tot.size == sectors_tot.size
    sector_dict = {}
    for tic in tics:
        sector_dict[tic] = str(list(np.unique(sectors_tot[np.in1d(tics_tot,tic)]).astype(int))).replace(' ','')
        
    return sector_dict


def get_disposition(disposition_num):
    '''Map the disposition to a str'''
    if disposition_num == 0:
        return 'pPC'
    elif disposition_num == 1:
        return 'PC'
    elif disposition_num == 2:
        return 'ST'
    elif disposition_num == 2.5:
        return 'pST'
    else:
        return np.nan
    
    
def append_to_list(planet_str):
    '''Append the PCs in the input planet string to the master list of planet 
    candidates from all sectors.'''
    # get file to append to
    fname = 'orion_PCs.txt'
    f = open(fname, 'r')
    g = f.read()
    f.close()

    # append 
    f = open(fname, 'w')
    f.write(g + ''.join(planet_str.split('disposition\n')[1:]))  # removes the hdr
    f.close()


if __name__ == '__main__':
    fname = 'PipelineResults_TIC_sector10/TIC_results_0_10000_det'
    self = loadpickle(fname)
    #append_to_list(get_PCs(self))
    get_CTOIs(self, 10)
