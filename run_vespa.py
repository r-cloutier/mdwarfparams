from LCclass import *
from vespa import FPPCalculation
import requests
from bs4 import BeautifulSoup

print('Ensure to load python v3.6 ($ source activate py36)')

global TESS_pixel_scale_arcsec
TESS_pixel_scale_arcsec = 20.69369686271056


def run_vespa_on_a_TIC(self, FWHMarcsec=0):
    # setup star file
    _setup_star_input(self)
    print('\nComputing the maximum separation to search for blends...')
    maxrad = get_EB_maxrad_condition(self) if FWHMarcsec == 0 else FWHMarcsec
    
    # run vespa on each planet candidate
    Nplanets = self.Ndet
    #vespa_results = np.zeros((Nplanets, 24))
    FPPs = np.zeros(Nplanets)
    for i in range(Nplanets):

        # setup light curve data file
        P,T0,rpRs,D = _get_planet_params(self,i)
        _setup_transitlightcurve_file(self,P,T0,D,i)

        # setup false positive probability file
        name = '%s_%i'%(self.object_name,i+1)
        maxoccdepth = get_EB_maxoccdepth_condition(self, D)
        _setup_fpp_input(self, name, P, rpRs, maxrad, maxoccdepth)  

        # run vespa
        print('\nRunning vespa on TIC %i...'%self.tic)
        os.system('starfit --all %s'%self.folder_full)
        os.system('calcfpp %s'%self.folder_full)

        # get read in results and compute the FPP of the planet candidate
        #vespa_results[i] = np.loadtxt('results.txt', skiprows=1)
        #FPPs[i] = vespa_results[i,-1]
        fpp = FPPCalculation.load(self.folder_full)
        FPPs[i] = fpp.FPP()
        
    # save FPPs
    np.save('%s/FPPs'%self.folder, FPPs)

    return FPPs


def _setup_star_input(self):
    # estimate imag from TESS
    assert hasattr(self, 'TESSmag')
    imag = _TESSmag2imag(self.TESSmag, self.Kmag)
    
    # read-in tempate star input file
    f = open('star_template.ini', 'r')
    g = f.read()
    f.close()

    # add stellar parameters
    g = g.replace('<<Teff>>', '%i'%self.Teff)
    g = g.replace('<<eTeff>>', '%i'%np.mean([self.ehi_Teff, self.elo_Teff]))
    g = g.replace('<<logg>>', '%.3f'%self.logg)
    g = g.replace('<<elogg>>', '%.3f'%np.mean([self.ehi_logg, self.elo_logg]))
    g = g.replace('<<par>>', '%.5f'%self.par)
    g = g.replace('<<epar>>', '%.5f'%self.e_par)
    g = g.replace('<<Jmag>>', '%.3f'%self.Jmag)
    g = g.replace('<<eJmag>>', '%.3f'%self.e_Jmag)
    g = g.replace('<<Hmag>>', '%.3f'%self.Hmag)
    g = g.replace('<<eHmag>>', '%.3f'%self.e_Hmag)
    g = g.replace('<<Kmag>>', '%.3f'%self.Kmag)
    g = g.replace('<<eKmag>>', '%.3f'%self.e_Kmag)
    g = g.replace('<<imag>>', '%.3f'%imag)
    
    # write new star input file
    h = open('%s/star.ini'%self.folder_full, 'w')
    h.write(g)
    h.close()


    
def _TESSmag2imag(TESSmag, Kmag):
    '''Convert the input TESS and K-band magnitudes to i band using Eq 7 from 
    Muirhead+2018, https://arxiv.org/abs/1710.00193.'''
    a, b, c = .6262, 1.0004, -.3326
    imag = (TESSmag - a + c*Kmag) / (b+c)
    return imag



def _setup_fpp_input(self, name, P, rpRs, maxrad, maxoccdepth):
    # read-in tempate fpp input file
    f = open('fpp_template.ini', 'r')
    g = f.read()
    f.close()

    # add stellar parameters
    g = g.replace('<<name>>', '%s'%name)
    g = g.replace('<<ra>>', '%.4f'%self.ra)
    g = g.replace('<<dec>>', '%.4f'%self.dec)
    g = g.replace('<<P>>', '%.6f'%P)
    g = g.replace('<<rpRs>>', '%.5f'%rpRs)
    g = g.replace('<<maxrad>>', '%.6f'%maxrad)
    g = g.replace('<<maxoccdepth>>', '%.6f'%maxoccdepth)
    
    # write new star input file
    h = open('%s/fpp.ini'%self.folder_full, 'w')
    h.write(g)
    h.close()



def _setup_transitlightcurve_file(self, P, T0, D, index):
    # phase fold the data in units of days from T0
    phase = foldAt(self.bjd, P, T0)
    phase[phase>.5] -= 1
    phase_days = phase*P

    # remove other planetary transits
    assert index < self.Ndet
    if self.Ndet > 1:
        inds = np.delete(np.arange(self.Ndet), int(index))
        other_transit_fmodel = np.product(self.fmodels[inds], 0)
        f = self.fcorr / other_transit_fmodel
    else:
        f = self.fcorr
        
    # focus on the LC close to transit
    g = (phase_days >= -3*D) & (phase_days <= 3*D)
    t, f, ef = phase_days[g], f[g], self.ef[g]
    
    # write light curve file
    s = np.argsort(t)
    np.savetxt('%s/transit.txt'%self.folder_full, np.array([t,f,ef]).T[s],
               delimiter='\t')



def _get_planet_params(self, index):
    assert index < self.Ndet
    i = int(index)
    P,T0,_,rpRs,_ = self.params_optimized[i]
    D = self.params_guess[i,-1]
    return P, T0, rpRs, D



def get_EB_maxrad_condition(self):
    '''
    Search target pixel files to get the FWHM to define the maximum separation 
    for vespa to consider when searching for blended stars.
    '''
    # read in target pixel files to find the average PSF FWHM
    bjds, tpfs = _read_TESS_TPF(self.tic)
    Ntimes, NpixX, NpixY = tpfs.shape
    x, y = np.meshgrid(np.arange(NpixX), np.arange(NpixY))
    FWHMs = np.zeros(0)
    for i in range(Ntimes):
        if np.any(np.isclose(self.bjd, bjds[i], atol=2.5/60/24)):
            img = tpfs[i]
            if np.all(np.isfinite(img)):
                gx, gy = np.where(img == np.nanmax(img))
                p0 = np.nanmax(img), float(gx), float(gy), 2, 2
                bnds_low = 0, 0, 0, 0, 0
                bnds_upp = p0[0]*10, NpixX, NpixY, NpixX, NpixY
                popt,_ = curve_fit(gaussian2D, (x,y), img.ravel(), p0=p0,
                                   bounds=(bnds_low,bnds_upp))
                sig = np.nanmean(popt[-2:])
                FWHMs = np.append(FWHMs, sig*2*np.sqrt(2*np.log(2)))

    # set maximum search radius to 1 PSF FWHM in arcseconds 
    FWHMs *= TESS_pixel_scale_arcsec
    np.save('%s/FWHMs_arcsec'%self.folder_full, FWHMs)
    maxrad = np.nanmedian(FWHMs)
    return maxrad


def _read_TESS_TPF(tic):
    # make directories
    try:
        os.mkdir('MAST')
    except OSError:
        pass
    try:
        os.mkdir('MAST/TESS')
    except OSError:
        pass
    folder2 = 'MAST/TESS/TIC%i'%tic
    try:
        os.mkdir(folder2)
    except OSError:
        pass
    
    # setup directories to fetch the data
    tic_str = '%.16d'%int(tic)
    tid1 = tic_str[:4]
    tid2 = tic_str[4:8]
    tid3 = tic_str[8:12]
    tid4 = tic_str[12:]

    # download tar file if not already
    fnames = np.array(glob.glob('%s/tess*tp.fits'%folder2))
    if fnames.size == 0:
        sectors = [1,2]
        for j in range(len(sectors)):
            sctr = 's%.4d'%sectors[j]
            url = 'https://archive.stsci.edu/missions/tess/tid/'
            folder = '%s/%s/%s/%s/%s/'%(sctr,tid1,tid2,tid3,tid4)
            fs = _listFD(url+folder, ext='_tp.fits')
            for i in range(fs.size):
                os.system('wget %s'%fs[i])
                fname = str(fs[i]).split('/')[-1]
                if os.path.exists(fname):
                    os.system('mv %s %s'%(fname, folder2))

    # read-in TPFs
    fname = np.array(glob.glob('%s/tess*tp.fits'%folder2))[0]
    hdu = fits.open(fname)[1].data
    bjds, TPFs = hdu['TIME']+2457000, hdu['FLUX']
    return bjds, TPFs
    

def _listFD(url, ext=''):
    '''List the contents of an https directory.'''
    page = requests.get(url).text
    soup = BeautifulSoup(page, 'html.parser')
    return np.array([url + '/' + node.get('href') \
                     for node in soup.find_all('a') \
                     if node.get('href').endswith(ext)])


def gaussian2D(arrs, A, mux, muy, sigx, sigy):
    xarr, yarr = arrs
    f = A*np.exp(-0.5*(((mux-xarr) / sigx)**2 + ((muy-yarr) / sigy)**2))
    return f.ravel()



def get_EB_maxoccdepth_condition(self, D, occdepth_upper_percentile=.95):
    '''use advice from crossfield 2016 on how to constrain eclipsing binaries 
    occultation depths (sect 6.1) 
    (http://adsabs.harvard.edu/abs/2016ApJS..226....7C)
    '''
    # take the resulting depths from the linear search that happen to be
    # outside of transit, and use the distribution of depths to set an eclipse
    # upper limit
    # first get out of transit depths from the linear search over durations
    fint = interp1d(self.bjd, np.product(self.fmodels,0))
    fullfmodel_at_transit_times = fint(self.transit_times)
    depths_out_transit=self.depths_linearsearch[fullfmodel_at_transit_times==1]
    
    # get the best duration
    g = abs(self.durations-D)==np.min(abs(self.durations-D))

    # return the Xth percentile of the out of transit depths 
    #plt.hist(depths_out_transit[:,g].flatten(),bins=1000,cumulative=1,normed=1)
    #plt.show()
    maxoccdepth = np.percentile(depths_out_transit[:,g].flatten(),
                                occdepth_upper_percentile)
    return maxoccdepth


if __name__ == '__main__':
    fs = np.array(glob.glob('PipelineResults_TIC/TIC_*/LC_-00099'))
    for i in range(fs.size):
        self = loadpickle(fs[i])
        print(i, self.tic)
        fwhm_fname = '%s/FWHMs_arcsec.npy'%self.folder_full
        FWHMarcsec = np.nanmedian(np.load(fwhm_fname)) if os.path.exists(fwhm_fname) else 0
        FPPs = run_vespa_on_a_TIC(self, FWHMarcsec=FWHMarcsec)
