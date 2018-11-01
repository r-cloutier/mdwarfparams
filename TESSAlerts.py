from planetsearch import *
from TESSLCclass import *


def get_TESS_alerts():
    '''get list of TESS alerts'''
    fname = 'hlsp_tess-data-alerts_tess_phot_alert-summary-s01+s02_tess_v3_spoc'
    d = np.loadtxt('input_data/%s.csv'%fname, delimiter=',',
                   usecols=np.delete(np.arange(13),2))
    TICid,TOIid,ra,dec,Tmag,eTmag,T0,eT0,P,eP,duration,eduration = d.T
    g = np.isfinite(P) & np.isfinite(T0) & np.isfinite(duration)
    TESSalert_dict = {'TICids': TICid[g], 'TOIids': TOIid[g], 'ras':ra[g],
                      'decs':dec[g], 'Tmags':Tmag[g], 'e_Tmags':eTmag[g],
                      'T0s':T0[g], 'e_T0s':eT0[g], 'Ps':P[g], 'e_Ps':eP[g],
                      'durations':duration, 'e_durations':eduration}
    return TESSalert_dict


def get_star(hdr):
    '''reader header info on the star'''
    Teff, logg, Rs = hdr['TEFF'], hdr['LOGG'], hdr['RADIUS']
    Ms = rvs.kg2Msun(10**(logg) * 1e-2 * rvs.Rsun2m(Rs)**2 / 6.67e-11)
    star_dict = {'ra':hdr['RA_OBJ'], 'dec':hdr['DEC_OBJ'],
                 'Tmag':hdr['TESSMAG'], 'Teff':Teff, 'logg':logg,
                 'M_H':hdr['MH'], 'Rs':Rs, 'Ms':Ms}
    return star_dict




def get_TESS_alert_LC(TICid, sectors=range(1,3)):
    '''download TESS alert LC from one of the sectors'''
    # create folders
    try:
        os.mkdir('MAST')
    except OSError:
        pass
    try:
        folder = 'MAST/TESSAlerts'
        os.mkdir(folder)
    except OSError:
        pass
    try:
        folder2 = '%s/TIC%i'%(folder, TICid)
        os.mkdir(folder2)
    except OSError:
        pass
        
    url = 'https://archive.stsci.edu/hlsps/tess-data-alerts'
    for i in sectors:
        fname = 'hlsp_tess-data-alerts_tess_phot_%.11d-s%.2d_tess_v1_lc.fits'%(TICid, i)
        os.system('wget %s/%s'%(url, fname))
        if os.path.exists(fname):
            os.system('mv %s %s'%(fname, folder2))
            break

    # read fits file
    hdu = fits.open('%s/%s'%(folder2, fname))
    assert len(hdu) > 1
    star_dict = get_star(hdu[0].header)
    bjd = hdu[1].data['TIME']
    f   = hdu[1].data['PDCSAP_FLUX']
    ef  = hdu[1].data['PDCSAP_FLUX_ERR']

    # clean: remove bad points and normalize
    g = (f != 0) & np.isfinite(bjd) & np.isfinite(f) & np.isfinite(ef)
    bjd, f, ef = bjd[g], f[g], ef[g]
    ef /= np.nanmedian(f) 
    f /= np.nanmedian(f)

    return star_dict, bjd, f, ef
    
    
def run_TESSalert_planetsearch(TICid, P, T0):
    '''given a TIC id, get the light and run the planet search algorithm'''
    # save stellar data and time-series
    star_dict, bjd, f, ef = get_TESS_alert_LC(TICid)
    self = TESSLC(TICid, -99)  # -99 is unique to planet_search
    self.Ptrue, self.T0true = P, T0
    self.bjd, self.f, self.ef = bjd, f, ef
    for attr in star_dict.keys():
        setattr(self, attr, star_dict[attr])
    self.DONE = False
    self._pickleobject()

    # fit initial GP
    thetaGPall, resultsGPall, thetaGPin, thetaGPout = do_optimize_0(bjd, f, ef)
    self.thetaGPall, self.resultsGPall = thetaGPall, resultsGPall
    self.thetaGPin, self.thetaGPout = thetaGPin, thetaGPout
    self._pickleobject()

    # search for transits in the corrected LC and get the transit parameters
    # guesses
    print 'Searching for transit-like events...\n'
    params, EBparams, maybeEBparams = find_transits(self, bjd, f, ef,
                                                    thetaGPout)
    self.params_guess = params
    self.params_guess_labels = np.array(['Ps', 'T0s', 'depths [Z]', \
                                         'durations [D]'])
    self.EBparams_guess, self.maybeEBparams_guess = EBparams, maybeEBparams
    self._pickleobject()

    # fit transit models to light curve
    run_mcmc(self)
    
    # is TESS alert planet detected
    self.is_TESSalert_detected = False
    if np.any(np.isclose(self.params_optimized[:,0], self.Ptrue, rtol=.02)):
        self.is_TESSalert_detected = True

    self.DONE = True
    self._pickleobject()


def do_i_run_this_star(TICid):
    # check if star is already done
    fname = 'PipelineResults/TIC_%i/TESSLC_-00099'%TICid
    if os.path.exists(fname):
        return not loadpickle(fname).DONE
    else:
        return True


if __name__ == '__main__':
    index = int(sys.argv[1])
    tess_dict = get_TESS_alerts()
    TICid = tess_dict['TICids'][index]
    P     = tess_dict['Ps'][index]
    T0    = tess_dict['T0s'][index]
    if do_i_run_this_star(TICid):
        run_TESSalert_planetsearch(TICid, P, T0)
