from occurrencerateclass import *


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
    fname = 'PipelineResults_TIC_sector9/TIC_results_0_10000_det'
    self = loadpickle(fname)
    append_to_list(get_PCs(self))
