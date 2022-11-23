from getdist.gaussian_mixtures import GaussianND
from getdist import plots
from getdist_plot import *
import os
os.chdir(os.path.dirname(os.path.realpath(__file__)))
### MCMC
# When providing the paths as a list remember to remove the extensions, and the number specifying the chain
## E.g.
# Suppose the first chain has the filename 2022-08-30_100000__1.txt, then only 2022-08-30_100000_ must be provided
# This is done for compatibility with getdist

#pathlist = ['../results/montepython_mcmc/nulcdm_photo_opt/2022-09-11_100000_','../results/montepython_mcmc/nulcdm_photo_opt/2022-09-10_100000_' ]
#mcmc1 = mcmc(pathlist=pathlist)

## Fishers
# Fisher CF
cf1 = Fisher('../results/cosmicfish_internal/photometric/optimistic/CosmicFish_v1.0_nulcdm_internal_class-Optimistic-3PT_WLGCph_fishermatrix.txt')

cf2 = Fisher('../results/cosmicfish_internal/photometric/optimistic/CosmicFish_v1.0_nulcdm_internal_camb-Optimistic-3PT_WLGCph_fishermatrix.txt')

mp1 = Fisher('../results/montepython_fisher/photometric/optimistic_HP/fisher.mat')

## Load Gaussians distributions correponding to the Fishers
gauss2=GaussianND(mean=cf2.fiducials, cov=cf2.fishermatrix,names=cf2.paramnames,labels=cf2.labels,is_inv_cov=True)
gauss1=GaussianND(mean=cf1.fiducials, cov=cf1.fishermatrix,names=cf1.paramnames,labels=cf1.labels,is_inv_cov=True)
gauss3=GaussianND(mean=mp1.fiducials, cov=mp1.fishermatrix,names=mp1.paramnames,labels=mp1.labels,is_inv_cov=True)

## Cosmo Triangle
g = plots.get_subplot_plotter(settings=plot_settings)
#g.triangle_plot([mcmc1.samples,gauss1,gauss2],filled=[True,False,False],params=cf1.cosmonames,legend_labels=['MCMC','CF-class','CF-camb'],contour_lws=[1,2,2])
g.triangle_plot([gauss1,gauss2,gauss3],filled=[False,False,False],params=cf1.cosmonames,legend_labels=['CF-class','CF-camb','MP-Fish'],contour_lws=[2,2,2])
g.export('./WLxGCPh_Opt_cosmo.pdf')
print('Cosmo triangle finished')
del g

## Nuisance Triangle
g = plots.get_subplot_plotter(settings=plot_settings)
#g.triangle_plot([mcmc1.samples,gauss1,gauss2],filled=[True,False,False],params=cf1.nuisancenames,legend_labels=['MCMC','CF-class','CF-camb'],contour_lws=[1,2,2])
g.triangle_plot([gauss1,gauss2,gauss3],filled=[False,False,False],params=cf1.nuisancenames,legend_labels=['CF-class','CF-camb','MP-Fish'],contour_lws=[2,2,2])
g.export('./WLxGCPh_Opt_nuisance.pdf')
print('Nuisance triangle finished')
del g

## Cosmo-Nuisance Rectangle
g = plots.get_subplot_plotter(settings=plot_settings)
#g.rectangle_plot(cf1.cosmonames,cf1.nuisancenames,roots=[mcmc1.samples,gauss1,gauss2],filled=[True,False,False],legend_labels=['MCMC','CF-class','CF-camb'],contour_lws=[1,2,2])
g.rectangle_plot(cf1.cosmonames,cf1.nuisancenames,roots=[gauss1,gauss2,gauss3],filled=[False,False,False],legend_labels=['CF-class','CF-camb','MP-Fish'],contour_lws=[2,2,2])
g.export('./WLxGCPh_Opt_cross.pdf')
print('Cross-terms rectangle finished')
del g
