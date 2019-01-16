from LCclass import *
from vespa import FPPCalculation
import requests
from bs4 import BeautifulSoup

print('Ensure to load python v3.6 ($ source activate py36)')

global TESS_pixel_scale_arcsec
TESS_pixel_scale_arcsec = 20.69369686271056


# set single transit stuff manually
#TICin,Pin,T0in,Din,planet_indin = 49678165,35.63695787,2458371.194515521,.3,0
#TICin,Pin,T0in,Din,planet_indin = 92444219,39.52254392,2458342.47285214,.35,2
#TICin,Pin,T0in,Din,planet_indin = 307210830,7.45111,2458355.2864,.04,0
#TICin,Pin,T0in,Din,planet_indin = 47484268,20.281385881512993,2458378.7845721757,.08,1
#TICin,Pin,T0in,Din,planet_indin = 278661431,17.63173026405275,2458343.929916539,.21,0
#TICin,Pin,T0in,Din,planet_indin = 100103201,7.4029000187292695,2458371.1104931794,0.09158793659674,0
TICin,Pin,T0in,Din,planet_indin = 231279823,5.976069768622068,2458361.641027355,.14,0
#TICin,Pin,T0in,Din,planet_indin = 415969908,15.811009507697745,2458381.0697858157,0.185690306,2

def run_vespa_on_a_TIC(self, planet_indices=[0], FWHMarcsec=0, binLC=False,
                       singletransit=False):
    # setup star file
    _setup_star_input(self)
    print('\nComputing the maximum separation to search for blends...')
    maxrad = get_EB_maxrad_condition(self) if FWHMarcsec == 0 else FWHMarcsec
    
    # run vespa on each planet candidate
    planet_indices = np.array([planet_indices]).reshape(len(planet_indices))
    Nplanets = planet_indices.size
    FPPs = np.zeros(Nplanets)
    for i in range(Nplanets):

        print(self.tic, i)

        # setup light curve data file
        if singletransit:
            assert TICin == self.tic
            P,T0,D = Pin,T0in,Din
        else:
            P,T0,_,D = self.params_guess[i]
        print(P, T0-2457000)
        _setup_transitlightcurve_file(self,P,T0,D,i,binLC)

        # setup false positive probability file
        name = '%s_%i'%(self.object_name,i+1)
        maxoccdepth = get_EB_maxoccdepth_condition(self, D)
        _,_,_,rpRs,_ = self.params_optimized[i]
        _setup_fpp_input(self, name, P, rpRs, maxrad, maxoccdepth)  

        # run vespa
        print('\nRunning vespa on TIC %i...'%self.tic)
        os.system('starfit --all %s'%self.folder_full)
        os.system('calcfpp %s'%self.folder_full)

        # get read in results and compute the FPP of the planet candidate
        #vespa_results[i] = np.loadtxt('results.txt', skiprows=1)
        #FPPs[i] = vespa_results[i,-1]
        try:
            fpp = FPPCalculation.load(self.folder_full)
            FPPs[i] = fpp.FPP()
            print('FPP', fpp.FPP())
        except FileNotFoundError:
            FPPs[i] = np.nan
        
    # save FPPs
    label = 'bin' if binLC else 'nobin'
    np.savetxt('%s/FPPs_%s.txt'%(self.folder_full,label), FPPs)

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



def _setup_transitlightcurve_file(self, P, T0, D, index, binLC=False):
    # bin the light curve?
    if binLC:
        bjd, fcorr, ef = boxcar(self.bjd, self.fcorr, self.ef, dt=D/6)
        label = 'binLC'
    else:
        bjd, fcorr, ef = self.bjd, self.fcorr, self.ef
        label = 'nobin'
        
    # phase fold the data in units of days from T0
    print('orbtial P', P)
    phase = foldAt(bjd, P, T0)
    phase[phase>.5] -= 1
    phase_days = phase*P

    # remove other planetary transits
    #assert index < self.Ndet
    #if self.Ndet > 1:
    #    inds = np.delete(np.arange(self.Ndet), int(index))
    #    other_transit_fmodel = np.product(self.fmodels[inds], 0)
    #    f = fcorr / other_transit_fmodel
    #else:
    f = fcorr
        
    # focus on the LC close to transit
    g = (phase_days >= -5*D) & (phase_days <= 5*D)
    #g = np.arange(f.size)
    tr, fr, efr = phase_days[g], f[g], ef[g]

    # plot LC
    plt.plot(phase_days, f, 'b.', alpha=.3)
    plt.errorbar(tr, fr, efr, fmt='k.', alpha=.5, elinewidth=.2)
    plt.xlabel('Phase in days [P=%.5f days]'%P)
    plt.ylabel('Normalized flux')
    plt.title('T0 = %.5f'%T0)
    plt.savefig('%s/vespaLC%i_%s.png'%(self.folder_full,index,label))
    plt.close('all')

    plt.plot(phase_days, f, 'b.', alpha=.3)
    plt.errorbar(tr, fr, efr, fmt='k.', alpha=.5, elinewidth=.2)
    plt.xlabel('Phase in days [P=%.5f days]'%P)
    plt.ylabel('Normalized flux')
    plt.title('T0 = %.5f'%T0)
    plt.xlim((-4*D,4*D))
    plt.savefig('%s/vespaLC_zoom%i_%s.png'%(self.folder_full,index,label))
    plt.close('all')
    
    # write light curve file
    s = np.argsort(tr)
    np.savetxt('%s/transit.txt'%self.folder_full, np.array([tr,fr,efr]).T[s],
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
    bjds, tpfs = read_TESS_TPF(self.tic)
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
                try:
                    popt,_ = curve_fit(gaussian2D, (x,y), img.ravel(), p0=p0,
                                       bounds=(bnds_low,bnds_upp))
                    sig = np.nanmean(popt[-2:])
                    FWHMs = np.append(FWHMs, sig*2*np.sqrt(2*np.log(2)))
                except RuntimeError:
                    fwhm = np.nanmedian(FWHMs) if np.any(np.isfinite(FWHMs)) \
                           else 50.
                    FWHMs = np.append(FWHMs, fwhm)

    # set maximum search radius to 1 PSF FWHM in arcseconds 
    FWHMs *= TESS_pixel_scale_arcsec
    np.save('%s/FWHMs_arcsec'%self.folder_full, FWHMs)
    maxrad = np.nanmedian(FWHMs)
    return maxrad


def read_TESS_TPF(tic):
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


def boxcar(t, f, ef, dt=.2, include_edges=False, tfull=np.zeros(0)):
    '''Boxbar bin the light curve.'''
    # check that the desired binning is coarser than the input sampling
    if np.diff(t).mean() > dt:
        if include_edges:
            assert tfull.size > 1
            tbin = np.append(np.append(tfull.min(), t), tfull.max())
            fbin = np.append(np.append(np.median(f[:10]), f), \
                             np.median(f[-10:]))
            efbin = np.append(np.append(np.median(ef[:10]), ef),
                              np.median(ef[-10:]))
            return tbin, fbin, efbin
        else:
            return t, f, ef

    Nbounds = int(np.floor((t.max()-t.min()) / dt))
    tbin, fbin, efbin = np.zeros(Nbounds-1), np.zeros(Nbounds-1), \
                        np.zeros(Nbounds-1)
    for i in range(Nbounds-1):
        inds = np.arange(t.size/Nbounds).astype(int) + int(i*t.size/Nbounds)
        tbin[i]  = t[inds].mean()
        fbin[i]  = np.median(f[inds])
        efbin[i] = np.mean(ef[inds]) / np.sqrt(inds.size)
    if include_edges:
        assert tfull.size > 1
        tbin = np.append(np.append(tfull.min(), tbin), tfull.max())
        fbin = np.append(np.append(np.median(f[:10]), fbin), np.median(f[-10:]))
        efbin = np.append(np.append(np.median(ef[:10]), efbin),
                          np.median(ef[-10:]))
    return tbin, fbin, efbin



def clean_vepsa_dir(directory):
    assert 'LC' not in directory
    os.system('rm %s/mist*'%directory)
    os.system('rm %s/FPP*'%directory)
    os.system('rm %s/calc*'%directory)
    os.system('rm %s/fpp.ini'%directory)
    os.system('rm %s/chains/*'%directory)
    os.system('rm %s/lhood*'%directory)
    os.system('rm %s/*.h5'%directory)
    os.system('rm %s/star*'%directory)
    os.system('rm %s/transit*'%directory)
    os.system('rm %s/trap*'%directory)
    os.system('rm %s/trsig*'%directory)
    os.system('rm %s/vespa*'%directory) 


if __name__ == '__main__':
    # get dispositions
    tics, disp, Ps = np.loadtxt('PipelineResults_TIC/PipelineResults_TIC_8d4_dispgeq0.txt').T
    tics2, Ps2 = np.loadtxt('PipelineResults_TIC/PipelineResults_TIC_8d4_planetPs.txt').T
    ticsSS, PsSS = np.loadtxt('PipelineResults_TIC/PipelineResults_TIC_8d4_planetPsSS.txt').T

    binLC = False
    
    #fs = np.array(glob.glob('PipelineResults_TIC/TIC_*/LC_-00099_8d4'))[23:]
    fs = np.array(['PipelineResults_TIC/TIC_231279823/LC_-00099_8d4'])
    for i in range(len(fs)):
        print(i, fs[i].split('/')[1].split('_')[-1])
        #if not os.path.exists('/'.join(fs[i].split('/')[:2])+'/FPPs.npy'):
    
        # get light curve
        self = loadpickle(fs[i])
        self.folder_full = self.folder_full.replace('PipelineResults_TIC_8d4','PipelineResults_TIC')
        self.fname_full = self.fname_full.replace('PipelineResults_TIC_8d4','PipelineResults_TIC')
        self.fname_full = self.fname_full.replace('LC_-00099','LC_-00099_8d4')

        # ensure that this system has a PC
        if np.any(disp[tics==self.tic] >= 0):

            # plot light curves for get the planet indices
            #print('all "detected" Ps', self.params_guess[:,0], 
		  #'\ncandidate periods', Ps[tics==self.tic],
                  #'\ncandidate period dispostions (should be >= 0)', disp[tics==self.tic])
            #plt.plot(self.bjd, self.f, '.')
            #plt.show()
            # enter 9 to skip
            #planet_indices = np.array(list(input('What are the planet indices? '))).astype(int)
	    #match_planet = np.isclose(self.params_guess[:,0],Ps2[tics2==self.tic],rtol=.02)
            #if np.any(match_planet):
            #    planet_indices = np.where(match_planet)[0]
            #else:
	    #    planet_indices = np.where(PsSS[ticsSS==self.tic])[0]
                
	    #singletransit = 2 in disp[tics==self.tic]
            #if singletransit:
            #    planet_indices = np.array([planet_indin])
            singletransit = True
            planet_indices = np.zeros(1) +2  #TEMP

            if 9 not in planet_indices:
                print('Running VESPA...')
                fwhm_fname = '%s/FWHMs_arcsec.npy'%self.folder_full
                FWHMarcsec = np.nanmedian(np.load(fwhm_fname)) \
                             if os.path.exists(fwhm_fname) else 0
                print(planet_indices)
                clean_vepsa_dir(self.folder_full)
                FPPs = run_vespa_on_a_TIC(self, planet_indices,
                                          FWHMarcsec=FWHMarcsec,  # TEMP
                                          binLC=binLC,
                                          singletransit=singletransit)

        else:
            pass
