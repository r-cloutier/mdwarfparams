# use the results from the Kepler-GAIA cross-match here: https://gaia-kepler.fun
from imports import *
from query_gaia_2MASS import *
from priors import get_results


def get_initial_KepID_data(fout):
    '''Get the available data for the Kepler-GAIA cross-matched stars.'''
    # get cross-matach results with different search radii
    fs = np.sort(glob.glob('input_data/Keplertargets/kepler_dr2_*fits'))

    # get data from the cross-match
    KepIDs = np.zeros(0)
    ras = np.zeros(0)
    decs = np.zeros(0)
    ls = np.zeros(0)
    bs = np.zeros(0)
    GBPmags = np.zeros(0)
    GRPmags = np.zeros(0) 
    e_GBPmags = np.zeros(0)
    e_GRPmags = np.zeros(0) 
    Kepmags = np.zeros(0)
    pars = np.zeros(0)
    e_pars = np.zeros(0)
    dists = np.zeros(0)
    ehi_dists = np.zeros(0)
    elo_dists = np.zeros(0)
    dist_modality = np.zeros(0)
    Teff = np.zeros(0)
    e_Teff = np.zeros(0)
    logg = np.zeros(0)
    e_logg = np.zeros(0)
    Rs = np.zeros(0)
    e_Rs = np.zeros(0)
    Jmags = np.zeros(0)
    Hmags = np.zeros(0)
    Kmags = np.zeros(0)
    for i in range(fs.size):
        
        f = fits.open(fs[i])[1].data
        KepIDs = np.append(KepIDs, f['kepid'])
        ras = np.append(ras, f['ra'])
        decs = np.append(decs, f['dec'])
        ls = np.append(ls, f['l'])
        bs = np.append(bs, f['b'])

        GBPmags = np.append(GBPmags, f['phot_bp_mean_mag'])
        FBP = f['phot_bp_mean_flux']
        eFBP = f['phot_bp_mean_flux_error']
        e_GBPmags = np.append(e_GBPmags, -2.5*np.log10(FBP / (FBP+eFBP)))

        GRPmags = np.append(GRPmags, f['phot_rp_mean_mag'])
        FRP = f['phot_rp_mean_flux']
        eFRP = f['phot_rp_mean_flux_error']
        e_GRPmags = np.append(e_GRPmags, -2.5*np.log10(FRP / (FRP+eFRP)))

        Kepmags = np.append(Kepmags, f['kepmag'])
        pars = np.append(pars, f['parallax']) + .029  # systemic correction
        e_pars = np.append(e_pars, f['parallax_error'])
        dists = np.append(dists, f['r_est'])
        ehi_dists = np.append(ehi_dists, f['r_hi']-f['r_est'])
        elo_dists = np.append(elo_dists, f['r_est']-f['r_lo'])
        dist_modality = np.append(dist_modality, f['r_modality_flag'])
        
        Teff = np.append(Teff, f['teff'])
        e_Teff = np.append(e_Teff, f['teff_err1'])
        logg = np.append(logg, f['logg'])
        e_logg = np.append(e_logg, f['logg_err1'])
        Rs = np.append(Rs, f['radius'])
        e_Rs = np.append(e_Rs, f['radius_err1'])
        Jmags = np.append(Jmags, f['jmag'])
        Hmags = np.append(Hmags, f['hmag'])
        Kmags = np.append(Kmags, f['kmag'])

    # remove duplicates from each search radius and make initial target list
    _,g1 = np.unique(KepIDs, return_index=True)
    print 'Number of stars in Kepler-GAIA catalog = %i'%g1.size
    g = (np.in1d(np.arange(Teff.size), g1)) & (Teff-e_Teff <= 4e3) & \
        (logg+e_logg > 3.5) & (Rs-e_Rs < .75)
    print 'Number of preliminary M dwarfs in Kepler-GAIA catalog = %i'%g.sum()
    KepIDs, ras, decs, ls, bs = KepIDs[g], ras[g], decs[g], ls[g], bs[g]
    GBPmags, GRPmags, e_GBPmags, e_GRPmags = GBPmags[g], GRPmags[g], \
                                             e_GBPmags[g], e_GRPmags[g]
    Kepmags, pars, e_pars = Kepmags[g], pars[g], e_pars[g]
    dists, ehi_dists, elo_dists = dists[g], ehi_dists[g], elo_dists[g]
    dist_modality = dist_modality[g]
    Teff, e_Teff = Teff[g], e_Teff[g]
    logg, e_logg = logg[g], e_logg[g]
    Rs, e_Rs = Rs[g], e_Rs[g]
    Jmags, Hmags, Kmags = Jmags[g], Hmags[g], Kmags[g]

    # get 2MASS uncertainties which are not included in the cross-match data
    e_Jmags, e_Hmags, e_Kmags = get_2MASS_Kep(ras, decs, Jmags, Hmags, Kmags)
    hdr = 'KepID,ra_deg,dec_deg,ls_deg,bs_deg,Kepmag,GBPmag,e_GBPmag,GRPmag,'+ \
          'e_GRPmag,parallax_mas,e_parallax,Jmag,e_Jmag,Hmag,e_Hmag,Kmag,'+ \
          'e_Kmag'
    outarr_tmp = np.array([KepIDs, ras, decs, ls, bs, Kepmags, GBPmags,
                           e_GBPmags, GRPmags, e_GRPmags, pars, e_pars, Jmags,
                           e_Jmags, Hmags, e_Hmags, Kmags, e_Kmags]).T
    np.savetxt(fout, outarr_tmp, delimiter=',', header=hdr, fmt='%.8e')   
    
    # get distance posteriors using Bailer-Jones R script
    distpost_success = save_posteriors(KepIDs, pars, e_pars, ls, bs, Kep=True)
    g = (distpost_success == True) & (dist_modality == 1)
    print 'Number of M dwarfs with reliable GAIA distances = %i'%g.sum()
    KepIDs, ras, decs, ls, bs = KepIDs[g], ras[g], decs[g], ls[g], bs[g]
    GBPmags, GRPmags, e_GBPmags, e_GRPmags = GBPmags[g], GRPmags[g], \
                                             e_GBPmags[g], e_GRPmags[g]
    Kepmags, pars, e_pars = Kepmags[g], pars[g], e_pars[g]
    dists, ehi_dists, elo_dists = dists[g], ehi_dists[g], elo_dists[g]
    #dist_modality = dist_modality[g]
    Teff, e_Teff = Teff[g], e_Teff[g]
    logg, e_logg = logg[g], e_logg[g]
    Rs, e_Rs = Rs[g], e_Rs[g]
    Jmags, Hmags, Kmags = Jmags[g], Hmags[g], Kmags[g]
    e_Jmags, e_Hmags, e_Kmags = e_Jmags[g], e_Hmags[g], e_Kmags[g]

    # compute parameter posteriors
    p  = compute_posterior_pdfs(KepIDs, ls, bs, GBPmags, e_GBPmags, GRPmags,
                                e_GRPmags, Jmags, e_Jmags, Hmags, e_Hmags,
                                Kmags, e_Kmags, Kep=True)
    mus, ehi_mus, elo_mus, dists, ehi_dists, elo_dists, AKs, e_AKs, MKs, \
        ehi_MKs, elo_MKs, Rss, ehi_Rss, elo_Rss, Teffs, ehi_Teffs, elo_Teffs, \
        Mss, ehi_Mss, elo_Mss, loggs, ehi_loggs, elo_loggs = p

    # identify M dwarfs
    g = (MKs>4.6) & (MKs < 9.8)
    print 'Number of M dwarfs in Kepler-GAIA catalog = %i'%g.sum()
    KepIDs, ras, decs, ls, bs = KepIDs[g], ras[g], decs[g], ls[g], bs[g]
    GBPmags, GRPmags, e_GBPmags, e_GRPmags = GBPmags[g], GRPmags[g], \
                                             e_GBPmags[g], e_GRPmags[g]
    Jmags, Hmags, Kmags = Jmags[g], Hmags[g], Kmags[g]
    e_Jmags, e_Hmags, e_Kmags = e_Jmags[g], e_Hmags[g], e_Kmags[g]
    Kepmags, pars, e_pars = Kepmags[g], pars[g], e_pars[g]
    dists, ehi_dists, elo_dists = dists[g], ehi_dists[g], elo_dists[g]
    mus, ehi_mus, elo_mus = mus[g], ehi_mus[g], elo_mus[g]
    AKs, e_AKs = AKs[g], e_AKs[g]
    MKs, ehi_MKs, elo_MKs = MKs[g], ehi_MKs[g], elo_MKs[g]
    Rss, ehi_Rss, elo_Rss = Rss[g], ehi_Rss[g], elo_Rss[g]
    Teffs, ehi_Teffs, elo_Teffs = Teffs[g], ehi_Teffs[g], elo_Teffs[g]
    Mss, ehi_Mss, elo_Mss = Mss[g], ehi_Mss[g], elo_Mss[g]
    loggs, ehi_loggs, elo_loggs = loggs[g], ehi_loggs[g], elo_loggs[g]
    
    # save results
    hdr = 'KepID,ra_deg,dec_deg,GBPmag,e_GBPmag,GRPmag,e_GRPmag,Kepmag,'+ \
          'Jmag,e_Jmag,Hmag,e_Hmag,Kmag,e_Kmag,parallax_mas,e_parallax,'+ \
          'dist_pc,ehi_dist,elo_dist,mu,ehi_mu,elo_mu,AK,e_AK,MK,ehi_MK,'+ \
          'elo_MK,Rs_RSun,ehi_Rs,elo_Rs,Teff_K,ehi_Teff,elo_Teff,Ms_MSun,'+ \
          'ehi_Ms,elo_Ms,logg_dex,ehi_logg,elo_logg'
    outarr = np.array([KepIDs, ras, decs, GBPmags, e_GBPmags, GRPmags,
                       e_GRPmags, Kepmags, Jmags, e_Jmags, Hmags, e_Hmags,
                       Kmags, e_Kmags, pars, e_pars, dists, ehi_dists,
                       elo_dists, mus, ehi_mus, elo_mus, AKs, e_AKs, MKs,
                       ehi_MKs, elo_MKs, Rss, ehi_Rss, elo_Rss, Teffs,
                       ehi_Teffs, elo_Teffs, Mss, ehi_Mss, elo_Mss, loggs,
                       ehi_loggs, elo_loggs]).T
    np.savetxt(fout, outarr, delimiter=',', header=hdr, fmt='%.8e')


def get_initial_EPIC_data(fout):
    '''Get the available data for the K2-GAIA cross-matched stars.'''
    # get cross-matach results with different search radii
    fs = np.sort(glob.glob('input_data/K2targets/k2_dr2_*fits'))

    # get data from the cross-match
    EPICs = np.zeros(0)
    ras = np.zeros(0)
    decs = np.zeros(0)
    ls = np.zeros(0)
    bs = np.zeros(0)
    GBPmags = np.zeros(0)
    GRPmags = np.zeros(0) 
    e_GBPmags = np.zeros(0)
    e_GRPmags = np.zeros(0) 
    Kepmags = np.zeros(0)
    pars = np.zeros(0)
    e_pars = np.zeros(0)
    dists = np.zeros(0)
    ehi_dists = np.zeros(0)
    elo_dists = np.zeros(0)
    dist_modality = np.zeros(0)
    Teff = np.zeros(0)
    e_Teff = np.zeros(0)
    logg = np.zeros(0)
    e_logg = np.zeros(0)
    Rs = np.zeros(0)
    e_Rs = np.zeros(0)
    for i in range(fs.size):
        
        f = fits.open(fs[i])[1].data
        EPICs = np.append(EPICs, f['epic_number'])
        ras = np.append(ras, f['ra'])
        decs = np.append(decs, f['dec'])
        ls = np.append(ls, f['l'])
        bs = np.append(bs, f['b'])

        GBPmags = np.append(GBPmags, f['phot_bp_mean_mag'])
        FBP = f['phot_bp_mean_flux']
        eFBP = f['phot_bp_mean_flux_error']
        e_GBPmags = np.append(e_GBPmags, -2.5*np.log10(FBP / (FBP+eFBP)))

        GRPmags = np.append(GRPmags, f['phot_rp_mean_mag'])
        FRP = f['phot_rp_mean_flux']
        eFRP = f['phot_rp_mean_flux_error']
        e_GRPmags = np.append(e_GRPmags, -2.5*np.log10(FRP / (FRP+eFRP)))

        Kepmags = np.append(Kepmags, f['k2_kepmag'])
        pars = np.append(pars, f['parallax']) + .029  # systemic correction
        e_pars = np.append(e_pars, f['parallax_error'])
        dists = np.append(dists, f['r_est'])
        ehi_dists = np.append(ehi_dists, f['r_hi']-f['r_est'])
        elo_dists = np.append(elo_dists, f['r_est']-f['r_lo'])
        dist_modality = np.append(dist_modality, f['r_modality_flag'])
        
        Teff = np.append(Teff, f['k2_teff'])
        e_Teff = np.append(e_Teff, f['k2_tefferr1'])
        logg = np.append(logg, f['k2_logg'])
        e_logg = np.append(e_logg, f['k2_loggerr1'])
        Rs = np.append(Rs, f['k2_rad'])
        e_Rs = np.append(e_Rs, f['k2_raderr1'])

    # remove duplicates from each search radius and make initial target list
    _,g1 = np.unique(EPICs, return_index=True)
    print 'Number of stars in K2-GAIA catalog = %i'%g1.size
    g = (np.in1d(np.arange(Teff.size), g1)) & (Teff-e_Teff <= 4e3) & \
        (logg+e_logg > 3.5) & (Rs-e_Rs < .75)
    print 'Number of preliminary M dwarfs in K2-GAIA catalog = %i'%g.sum()
    EPICs, ras, decs, ls, bs = EPICs[g], ras[g], decs[g], ls[g], bs[g]
    GBPmags, GRPmags, e_GBPmags, e_GRPmags = GBPmags[g], GRPmags[g], \
                                             e_GBPmags[g], e_GRPmags[g]
    Kepmags, pars, e_pars = Kepmags[g], pars[g], e_pars[g]
    dists, ehi_dists, elo_dists = dists[g], ehi_dists[g], elo_dists[g]
    dist_modality = dist_modality[g]
    Teff, e_Teff = Teff[g], e_Teff[g]
    logg, e_logg = logg[g], e_logg[g]
    Rs, e_Rs = Rs[g], e_Rs[g]

    # get 2MASS uncertainties which are not included in the cross-match data
    Jmags,e_Jmags,Hmags,e_Hmags,Kmags,e_Kmags = get_2MASS_K2(ras, decs, GBPmags, GRPmags)
    hdr = 'EPIC,ra_deg,dec_deg,ls_deg,bs_deg,Kepmag,GBPmag,e_GBPmag,GRPmag,'+ \
          'e_GRPmag,parallax_mas,e_parallax,Jmag,e_Jmag,Hmag,e_Hmag,Kmag,'+ \
          'e_Kmag'
    outarr_tmp = np.array([EPICs, ras, decs, ls, bs, Kepmags, GBPmags,
                           e_GBPmags, GRPmags, e_GRPmags, pars, e_pars, Jmags,
                           e_Jmags, Hmags, e_Hmags, Kmags, e_Kmags]).T
    np.savetxt(fout, outarr_tmp, delimiter=',', header=hdr, fmt='%.8e')   
    
    # get distance posteriors using Bailer-Jones R script
    distpost_success = save_posteriors(EPICs, pars, e_pars, ls, bs, K2=True)
    g = (distpost_success == True) & (dist_modality == 1)
    print 'Number of M dwarfs with reliable GAIA distances = %i'%g.sum()
    EPICs, ras, decs, ls, bs = EPICs[g], ras[g], decs[g], ls[g], bs[g]
    GBPmags, GRPmags, e_GBPmags, e_GRPmags = GBPmags[g], GRPmags[g], \
                                             e_GBPmags[g], e_GRPmags[g]
    Kepmags, pars, e_pars = Kepmags[g], pars[g], e_pars[g]
    dists, ehi_dists, elo_dists = dists[g], ehi_dists[g], elo_dists[g]
    #dist_modality = dist_modality[g]
    Teff, e_Teff = Teff[g], e_Teff[g]
    logg, e_logg = logg[g], e_logg[g]
    Rs, e_Rs = Rs[g], e_Rs[g]
    Jmags, Hmags, Kmags = Jmags[g], Hmags[g], Kmags[g]
    e_Jmags, e_Hmags, e_Kmags = e_Jmags[g], e_Hmags[g], e_Kmags[g]

    # compute parameter posteriors
    p  = compute_posterior_pdfs(EPICs, ls, bs, GBPmags, e_GBPmags, GRPmags,
                                e_GRPmags, Jmags, e_Jmags, Hmags, e_Hmags,
                                Kmags, e_Kmags, K2=True)
    mus, ehi_mus, elo_mus, dists, ehi_dists, elo_dists, AKs, e_AKs, MKs, \
        ehi_MKs, elo_MKs, Rss, ehi_Rss, elo_Rss, Teffs, ehi_Teffs, elo_Teffs, \
        Mss, ehi_Mss, elo_Mss, loggs, ehi_loggs, elo_loggs = p

    # identify M dwarfs
    g = (MKs>4.6) & (MKs < 9.8)
    print 'Number of M dwarfs in K2-GAIA catalog = %i'%g.sum()
    EPICs, ras, decs, ls, bs = EPICs[g], ras[g], decs[g], ls[g], bs[g]
    GBPmags, GRPmags, e_GBPmags, e_GRPmags = GBPmags[g], GRPmags[g], \
                                             e_GBPmags[g], e_GRPmags[g]
    Jmags, Hmags, Kmags = Jmags[g], Hmags[g], Kmags[g]
    e_Jmags, e_Hmags, e_Kmags = e_Jmags[g], e_Hmags[g], e_Kmags[g]
    Kepmags, pars, e_pars = Kepmags[g], pars[g], e_pars[g]
    dists, ehi_dists, elo_dists = dists[g], ehi_dists[g], elo_dists[g]
    mus, ehi_mus, elo_mus = mus[g], ehi_mus[g], elo_mus[g]
    AKs, e_AKs = AKs[g], e_AKs[g]
    MKs, ehi_MKs, elo_MKs = MKs[g], ehi_MKs[g], elo_MKs[g]
    Rss, ehi_Rss, elo_Rss = Rss[g], ehi_Rss[g], elo_Rss[g]
    Teffs, ehi_Teffs, elo_Teffs = Teffs[g], ehi_Teffs[g], elo_Teffs[g]
    Mss, ehi_Mss, elo_Mss = Mss[g], ehi_Mss[g], elo_Mss[g]
    loggs, ehi_loggs, elo_loggs = loggs[g], ehi_loggs[g], elo_loggs[g]
    
    # save results
    hdr = 'EPIC,ra_deg,dec_deg,GBPmag,e_GBPmag,GRPmag,e_GRPmag,Kepmag,'+ \
          'Jmag,e_Jmag,Hmag,e_Hmag,Kmag,e_Kmag,parallax_mas,e_parallax,'+ \
          'dist_pc,ehi_dist,elo_dist,mu,ehi_mu,elo_mu,AK,e_AK,MK,ehi_MK,'+ \
          'elo_MK,Rs_RSun,ehi_Rs,elo_Rs,Teff_K,ehi_Teff,elo_Teff,Ms_MSun,'+ \
          'ehi_Ms,elo_Ms,logg_dex,ehi_logg,elo_logg'
    outarr = np.array([EPICs, ras, decs, GBPmags, e_GBPmags, GRPmags,
                       e_GRPmags, Kepmags, Jmags, e_Jmags, Hmags, e_Hmags,
                       Kmags, e_Kmags, pars, e_pars, dists, ehi_dists,
                       elo_dists, mus, ehi_mus, elo_mus, AKs, e_AKs, MKs,
                       ehi_MKs, elo_MKs, Rss, ehi_Rss, elo_Rss, Teffs,
                       ehi_Teffs, elo_Teffs, Mss, ehi_Mss, elo_Mss, loggs,
                       ehi_loggs, elo_loggs]).T
    np.savetxt(fout, outarr, delimiter=',', header=hdr, fmt='%.8e')

        
def get_2MASS_Kep(ras_deg, decs_deg, Jmags, Hmags, Kmags,
                  radius_deg=.017, phot_rtol=.02):
    '''Match Kepler stars with GAIA data to the 2MASS point-source catlog to 
    retrieve photometric uncertainties.'''
    # get 2MASS data
    d = np.load('input_data/Keplertargets/fp_2mass.fp_psc12298.npy')
    inds = np.array([0,1,3,5,6,8,9,11])
    ras2M, decs2M, J2M, eJ2M, H2M, eH2M, K2M, eK2M = d[:,inds].T
    
    # match each star individually
    Nstars = ras_deg.size
    e_Jmags, e_Hmags, e_Kmags = np.zeros(Nstars), np.zeros(Nstars), \
                                np.zeros(Nstars)
    print 'Getting 2MASS photometry...'
    for i in range(Nstars):

        if i % 1e2 == 0:
            print float(i) / Nstars
        
        # get matching photometry between Kepler-GAIA and 2MASS
        g = (ras2M >= ras_deg[i] - radius_deg) & \
            (ras2M <= ras_deg[i] + radius_deg) & \
            (decs2M >= decs_deg[i] - radius_deg) & \
            (decs2M <= decs_deg[i] + radius_deg) & \
            np.isclose(J2M, Jmags[i], rtol=phot_rtol) & \
            np.isclose(H2M, Hmags[i], rtol=phot_rtol) & \
            np.isclose(K2M, Kmags[i], rtol=phot_rtol)

        if g.sum() > 0:
            g2 = (abs(J2M[g]-Jmags[i]) == np.min(abs(J2M[g]-Jmags[i]))) & \
                 (abs(K2M[g]-Kmags[i]) == np.min(abs(K2M[g]-Kmags[i])))
            e_Jmags[i] = eJ2M[g][g2]
            e_Hmags[i] = eH2M[g][g2]
            e_Kmags[i] = eK2M[g][g2]

        else:
            e_Jmags[i], e_Hmags[i], e_Kmags[i] = np.repeat(np.nan, 3)

    return e_Jmags, e_Hmags, e_Kmags 



def get_2MASS_K2(ras_deg, decs_deg, GBPmags, GRPmags,
                 radius_deg=.017):
    '''Match K2 stars with GAIA data to the 2MASS point-source catlog to 
    retrieve photometry and uncertainties.'''
    # get 2MASS point source catalog
    d = np.load('input_data/Keplertargets/fp_2mass.fp_psc12298.npy')
    inds = np.array([0,1,3,5,6,8,9,11])
    ras2M, decs2M, J2M, eJ2M, H2M, eH2M, K2M, eK2M = d[:,inds].T
    
    # match each star individually
    Nstars = ras_deg.size
    GBP_GRP = GBPmags - GRPmags
    Jmags, Hmags, Kmags = np.zeros(Nstars), np.zeros(Nstars), np.zeros(Nstars)
    e_Jmags, e_Hmags, e_Kmags = np.zeros(Nstars), np.zeros(Nstars), \
                                np.zeros(Nstars)
    print 'Getting 2MASS photometry...'
    for i in range(Nstars):

        if i % 1e2 == 0:
            print float(i) / Nstars
        
        # get matching photometry between Kepler-GAIA and 2MASS
        H_K = GAIAcolor2HK(GBP_GRP[i])
        g = (ras2M >= ras_deg[i] - radius_deg) & \
            (ras2M <= ras_deg[i] + radius_deg) & \
            (decs2M >= decs_deg[i] - radius_deg) & \
            (decs2M <= decs_deg[i] + radius_deg) & \
            np.isclose(H2M-K2M, H_K, atol=3*.2411) # evans18 color relation

        if g.sum() > 0:
            g2 = (abs((H2M-K2M)[g]-H_K) == np.min(abs((H2M-K2M)[g]-H_K)))
            Jmags[i] = J2M[g][g2]
            Hmags[i] = H2M[g][g2]
            Kmags[i] = K2M[g][g2]
            e_Jmags[i] = eJ2M[g][g2]
            e_Hmags[i] = eH2M[g][g2]
            e_Kmags[i] = eK2M[g][g2]

        else:
            Jmags[i], Hmags[i], Kmags[i] = np.repeat(np.nan, 3)
            e_Jmags[i], e_Hmags[i], e_Kmags[i] = np.repeat(np.nan, 3)

    return Jmags, e_Jmags, Hmags, e_Hmags, e_Kmags, Kmags 


def GAIAcolor2HK(GBP_GRP):
    '''Use color relation from Evans18 to map GBP-GRP to H_K
    (https://www.aanda.org/articles/aa/pdf/2018/08/aa32756-18.pdf)'''
    a, b, c = 1.991, 6.098, .4238 - GBP_GRP
    H_K = (-b + np.sqrt(b*b - 4*a*c)) / (2*a) 
    return H_K


def save_posteriors(IDnums, pars, e_pars, ls, bs, Kep=False, K2=False,
                    TESS=False):
    '''Go to the directory with the Bailor-Jones + 2018 R scripts to compute 
    the distance posterior for a single source. Save the posterior with a 
    unique file name.'''
    assert IDnums.size == pars.size
    assert IDnums.size == e_pars.size
    assert IDnums.size == ls.size
    assert IDnums.size == bs.size

    if Kep:
        prefix = 'KepID'
    elif K2:
        prefix = 'EPIC'
    elif TESS:
        prefix = 'TIC'
    
    cwd = os.getcwd()
    os.chdir('%s/Gaia-DR2-distances_custom'%cwd)
    cmd_prefix = 'Rscript get_dist_post.R %s'%prefix
    distpost_success = np.zeros(IDnums.size).astype(bool)
    for i in range(IDnums.size):
        print float(i) / IDnums.size
        cmd = '%s_%i %.6e %.6e %.6f %.6f'%(cmd_prefix, IDnums[i], pars[i],
                                           e_pars[i], ls[i], bs[i])
        # run if not done already
        fout = 'DistancePosteriors/%s_%i.csv'%(prefix, IDnums[i])
        if not os.path.exists(fout):
            os.system(cmd)
        if os.path.exists(fout):
            distpost_success[i] = True
            
    os.chdir(cwd)
    return distpost_success



def compute_posterior_pdfs(IDnums, ls, bs, GBPmags, e_GBPmags, GRPmags,
                           e_GRPmags, Jmags, e_Jmags, Hmags, e_Hmags, Kmags,
                           e_Kmags, Nsamp=1e3, Kep=False, K2=False):
    '''Read-in the distance posterior given the identifier and compute 
    parameters of interest given input photometry.'''
    if Kep:
        prefix = 'KepID'
    elif K2:
        prefix = 'EPIC'
    
    N = IDnums.size
    mus, ehi_mus, elo_mus = np.zeros(N), np.zeros(N), np.zeros(N)
    dists, ehi_dists, elo_dists = np.zeros(N), np.zeros(N), np.zeros(N)
    AKs, e_AKs = np.zeros(N), np.zeros(N)
    MKs, ehi_MKs, elo_MKs = np.zeros(N), np.zeros(N), np.zeros(N)
    Rss, ehi_Rss, elo_Rss = np.zeros(N), np.zeros(N), np.zeros(N)
    Teffs, ehi_Teffs, elo_Teffs = np.zeros(N), np.zeros(N), np.zeros(N)
    Mss, ehi_Mss, elo_Mss = np.zeros(N), np.zeros(N), np.zeros(N)
    loggs, ehi_loggs, elo_loggs = np.zeros(N), np.zeros(N), np.zeros(N)
    for i in range(N):

        print float(i)/N
        # get dist pdf
        try:
            fname='Gaia-DR2-distances_custom/DistancePosteriors/%s_%i.csv'%(prefix,
                                                                            IDnums[i])
            x_dist, pdf_dist = np.loadtxt(fname, delimiter=',', skiprows=1,
                                          usecols=(1,2)).T
        except IOError:
            pass

        # sample parameter distributions
        Nsamp = int(Nsamp)
        samp_dist = np.random.choice(x_dist, Nsamp, p=pdf_dist/pdf_dist.sum())
        dists[i], ehi_dists[i], elo_dists[i] = \
                                        get_results(samp_dist.reshape(Nsamp,1))
        samp_mu = 5*np.log10(samp_dist) - 5
        AK_tmp = compute_AK_mwdust(ls[i], bs[i], dists[i],
                                   np.mean([ehi_dists[i],elo_dists[i]]))
        AKs[i], e_AKs[i] = unp.nominal_values(AK_tmp), unp.std_devs(AK_tmp)
        samp_AK = np.random.normal(AKs[i], e_AKs[i], Nsamp)
        samp_Kmag = np.random.normal(Kmags[i], e_Kmags[i], Nsamp)
        samp_MK = samp_Kmag - samp_mu - samp_AK
        samp_Rs = sample_Rs_from_MK(samp_MK)
        samp_GBPmag = np.random.normal(GBPmags[i], e_GBPmags[i], Nsamp)
        samp_GRPmag = np.random.normal(GRPmags[i], e_GRPmags[i], Nsamp)
        samp_Jmag = np.random.normal(Jmags[i], e_Jmags[i], Nsamp)
        samp_Hmag = np.random.normal(Hmags[i], e_Hmags[i], Nsamp)
        samp_Teff = sample_Teff_from_colors(samp_GBPmag, samp_GRPmag,
                                            samp_Jmag, samp_Hmag)
        samp_Ms = sample_Ms_from_MK(samp_MK)
        samp_logg = sample_logg(samp_Ms, samp_Rs)

        # get point estimates
        mus[i], ehi_mus[i], elo_mus[i] = get_results(samp_mu.reshape(Nsamp,1))
        MKs[i], ehi_MKs[i], elo_MKs[i] = get_results(samp_MK.reshape(Nsamp,1))
        Rss[i], ehi_Rss[i], elo_Rss[i] = get_results(samp_Rs.reshape(Nsamp,1))
        Teffs[i], ehi_Teffs[i], elo_Teffs[i] = \
                                        get_results(samp_Teff.reshape(Nsamp,1))
        Mss[i], ehi_Mss[i], elo_Mss[i] = get_results(samp_Ms.reshape(Nsamp,1))
        loggs[i], ehi_loggs[i], elo_loggs[i] = \
                                        get_results(samp_logg.reshape(Nsamp,1))

        # save posteriors
        hdr = 'GBPmag,GRPmag,Jmag,Hmag,Kmag,dist_pc,mu,AK,MK,Rs,Teff,Ms,logg'
        outarr = np.array([samp_GBPmag, samp_GRPmag, samp_Jmag, samp_Hmag,
                           samp_Kmag, samp_dist, samp_mu, samp_AK, samp_MK,
                           samp_Rs, samp_Teff, samp_Ms, samp_logg]).T
        fname='Gaia-DR2-distances_custom/DistancePosteriors/KepID_allpost_%i'%IDnums[i]
        np.savetxt(fname, outarr, delimiter=',', header=hdr, fmt='%.8e')
        
    return mus, ehi_mus, elo_mus, dists, ehi_dists, elo_dists, AKs, e_AKs, \
        MKs, ehi_MKs, elo_MKs, Rss, ehi_Rss, elo_Rss, Teffs, ehi_Teffs, \
        elo_Teffs, Mss, ehi_Mss, elo_Mss, loggs, ehi_loggs, elo_loggs


def sample_Rs_from_MK(samp_MK):
    '''Use relation from Mann+2015 (table 1)'''
    a, b, c, Rs_sigma_frac = 1.9515, -.3520, .01680, .0289
    p = np.poly1d((c,b,a))
    samp_MK_tmp = np.copy(samp_MK)
    samp_MK_tmp[(samp_MK<=4.6) | (samp_MK>=9.8)] = np.nan
    samp_Rs = p(samp_MK_tmp)
    samp_Rs += np.random.normal(0, samp_Rs*Rs_sigma_frac, samp_MK.size)
    return samp_Rs


def sample_Ms_from_MK(samp_MK):
    '''Use relation from Benedict+2016'''
    c0 = np.random.normal(.2311, 4e-4, samp_MK.size)
    c1 = np.random.normal(-.1352, 7e-4, samp_MK.size)
    c2 = np.random.normal(.04, 5e-4, samp_MK.size)
    c3 = np.random.normal(.0038, 2e-4, samp_MK.size)
    c4 = np.random.normal(-.0032, 1e-4, samp_MK.size)
    samp_MK_tmp = np.copy(samp_MK)
    samp_MK_tmp[(samp_MK<=4.6) | (samp_MK>10)] = np.nan
    samp_MK_tmp[samp_MK>=10] = np.nan
    dMK = samp_MK_tmp - 7.5
    samp_Ms = c0 + c1*dMK + c2*dMK**2 + c3*dMK**3 + c4*dMK**4
    return samp_Ms
    

def sample_Teff_from_colors(samp_GBPmag, samp_GRPmag, samp_Jmag, samp_Hmag,
                            Teff_scatter=49):
    '''Use the relation from Mann+2015 (table 2)'''
    a, b, c, d, e, f, g = 3.172, -2.475, 1.082, -.2231, .01738, .08776, -.04355
    pG = np.poly1d((e,d,c,b,a))
    p2 = np.poly1d((g,f,0))
    samp_Teff = 35e2 * (pG(samp_GBPmag-samp_GRPmag) + p2(samp_Jmag-samp_Hmag)) \
                + np.random.normal(0, Teff_scatter, samp_Jmag.size)
    return samp_Teff


def sample_logg(samp_Ms, samp_Rs):
    G = 6.67e-11
    samp_logg = np.log10(G*rvs.Msun2kg(samp_Ms)*1e2 / rvs.Rsun2m(samp_Rs)**2)
    return samp_logg


if __name__ == '__main__':
    #fout = 'input_data/Keplertargets/KepMdwarfsv10.csv'
    #get_initial_KepID_data(fout)
    fout = 'input_data/K2targets/K2Mdwarfsv10.csv'
    get_initial_EPIC_data(fout)
