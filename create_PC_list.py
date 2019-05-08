from OccurrenceRateclass import *


def get_PCs(folder):
    '''Given a folder with an OccurrenceRateclass object, get all the planet 
    candidates from orion and with human dispositions, so that they can be 
    added to the list of planet candidates.'''

    # get the results object
    self = loadpickle('%s/TIC_results_0_10000_det'%folder)
    assert hasattr(self, 'tics')
    assert hasattr(self, 'disposition_human')

    # isolate planet candidates, putative PCs, and single transits
    g01 = np.in1d(self.disposition_human, [0,1])  # PCs and pPCs
    g2 = self.disposition_human >= 2  # STs and pSTs

    # get columns of interest
    # sector,TIC,ra,dec,Tmag,Jmag,Hmag,Kmag,dist,Ms,Rs,Teff,P,T0,rp/Rs,rp
    outarr = 


    
def append_to_list(planet_arr):
    '''     append them to the 
    list of planet candidates from all sectors.'''
