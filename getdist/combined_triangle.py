from getdist.gaussian_mixtures import GaussianND
from getdist import plots
from getdist_plot import *


def combined_triangle(case1,case2):

    if case1 == 'pessimistic':
        case1_cap='Pessimistic'
        case1_short='pess'
    elif case1 == 'optimistic':
        case1_cap='Optimistic'
        case1_short='opt'

    if case2 == 'pessimistic':
        case2_cap='Pessimistic'
        case2_short='pess'
    elif case2 == 'optimistic':
        case2_cap='Optimistic'
        case2_short='opt'

    ## Fishers
    cf1_path = '../results/cosmicfish_internal/photometric/%s/CosmicFish_v1.0_nulcdm_internal_class-%s-3PT_WLGCph_fishermatrix'%(case1,case1_cap)
    cf2_path = '../results/cosmicfish_internal/spectroscopic/%s/CosmicFish_v1.0_nulcdm_internal_class-%s-own_GCsp_fishermatrix'%(case2,case2_cap)
    cf3_path = '../results/cosmicfish_internal/combined/%s_%s/CosmicFish_v1.0_nulcdm_internal_class-combined-%s-%s_fishermatrix'%(case1_short,case2_short,case1_short,case2_short)

    cf1 = Fisher(cf1_path+'.txt')
    cf2 = Fisher(cf2_path+'.txt')
    cf3 = Fisher(cf3_path+'.txt')

    mp1_path = '../results/montepython_fisher/photometric/%s_HP/fisher'%(case1)
    mp2_path = '../results/montepython_fisher/spectroscopic/%s_HP/fisher'%(case2)
    mp3_path = '../results/montepython_fisher/combined/%s_%s/fisher'%(case1_short,case2_short)

    mp1 = Fisher(mp1_path+'.mat')
    mp2 = Fisher(mp2_path+'.mat')
    mp3 = Fisher(mp3_path+'.txt')

    ## Load Gaussians distributions correponding to the Fishers
    gauss1=GaussianND(mean=cf1.fiducials, cov=cf1.fishermatrix, names=cf1.paramnames,labels=cf1.labels,is_inv_cov=True)
    gauss2=GaussianND(mean=cf2.fiducials, cov=cf2.fishermatrix, names=cf2.paramnames,labels=cf2.labels,is_inv_cov=True)
    gauss3=GaussianND(mean=cf3.fiducials, cov=cf3.fishermatrix, names=cf3.paramnames,labels=cf3.labels,is_inv_cov=True)

    gauss4=GaussianND(mean=cf1.fiducials, cov=mp1.fishermatrix, names=mp1.paramnames,labels=cf1.labels,is_inv_cov=True)
    gauss5=GaussianND(mean=cf2.fiducials, cov=mp2.fishermatrix, names=mp2.paramnames,labels=cf2.labels,is_inv_cov=True)
    gauss6=GaussianND(mean=cf3.fiducials, cov=mp3.fishermatrix, names=mp3.paramnames,labels=cf3.labels,is_inv_cov=True)

    ## Cosmo Triangle
    g = plots.get_subplot_plotter(settings=plot_settings)
    g.triangle_plot([gauss1,gauss2,gauss3,gauss4,gauss5,gauss6],filled=[False,False,False,False,False,False],params=cf1.cosmonames,legend_labels=['CosmicFish: photo/%s'%(case1_short),'CosmicFish: spectro/%s'%(case2_short),'CosmicFish: combined','MontePython/Fisher: photo/%s'%(case1_short),'MontePython/Fisher: spectro/%s'%(case2_short),'MontePython/Fisher: combined'],contour_lws=[1,1,2,1,1,1,1],contour_ls=[':',':',':','-','-','-'],contour_colors=['k','b','r','k','b','r'])
    g.export('./combined_%s_%s_cosmo.pdf'%(case1_short,case2_short))
    print('Combined %s %s triangle finished'%(case1_short,case2_short))
    del g

    return

combined_triangle('pessimistic','pessimistic')
combined_triangle('pessimistic','optimistic')
combined_triangle('optimistic','pessimistic')
combined_triangle('optimistic','optimistic')
