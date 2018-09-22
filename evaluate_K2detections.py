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
    
    # check which planets are detected and which are FPs
    det, FP = np.zeros(0, dtype=bool), np.zeros(0, dtype=bool) 
    epicnamesout, cond_vals = np.zeros(0, dtype=int), np.zeros((0, 5))
    for i in self.unique_inds:

        print self.epicnames[i]
        g = self.epicnames == self.epicnames[i]
        Ndet = g.sum() - 1  # remove nan entry
        if Ndet == 0:
            det = np.append(det, False)
            FP = np.append(FP, False)
            epicnamesout = np.append(epicnamesout, self.epicnames[i])
            cond_vals = np.append(cond_vals, np.repeat(np.nan,5).reshape(1,5),
                                  0)
        
        Ps_det = self.params_guess[g][1:,0]  # skip the first nan entry

        g2 = epicnames_known == self.epicnames[i]
        assert g2.sum() > 0
        
        for j in range(Ndet):
                    
            isdet = np.isclose(Ps_det[j], Ps_known[g2], rtol=.05)

            if np.any(isdet):
                det = np.append(det, True)
                FP = np.append(FP, False)
            else:
                det = np.append(det, False)
                FP = np.append(FP, True)

            epicnamesout = np.append(epicnamesout, self.epicnames[i])
            cond_vals = np.append(cond_vals,
                                  self.cond_vals[g][j+1].reshape(1,5), 0)

    return epicnamesout, det, FP, cond_vals
