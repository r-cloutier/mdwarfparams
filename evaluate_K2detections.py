'''
Given the results of pipeline runs from "compile_results.py", compare the 
detected planets to the known K2 transiting planets to get a sense of the 
pipeline performance.
'''
from compile_results import *


def get_known_K2():
    _,epicnames,Ps = np.loadtxt('input_data/K2targets/' + \
                                'K2knownMdwarfplanetsv2.csv').T
    return epicnames, Ps


def compare_det_and_FPs(K2resultsclass):
    # get answers
    epicnames_known, Ps_known = get_known_K2()

    self = K2resultsclass
    
    # check which planets are detected
    det = np.zeros(Ps_known.size, dtype=bool)
    for i in range(Ps_known.size):

        g = self.epicnames == epicnames_known[i]
        print epicnames_known[i]
        if g.sum() == 0:
            det[i] = np.nan

        Ps_det = self.params_guess[g][1:,0]  # skip first entry which is nan

        if np.any(np.isclose(Ps_det, Ps_known[i], rtol=.05)):
            det[i] = True

        else:
            det[i] = False
            

    # now check for any false positives
    NFPs = np.zeros(Ps_known.size)
    for i in range(Ps_known.size):

        g = self.epicnames == epicnames_known[i]
        if g.sum() == 0:
            Ndet = np.nan
        else:
            Ndet = self.params_guess[g][1:,0].size
        
        g = epicnames_known == epicnames_known[i]
        Ndet_correct = det[g].sum()

        NFPs[i] = Ndet - Ndet_correct

    return det, NFPs

