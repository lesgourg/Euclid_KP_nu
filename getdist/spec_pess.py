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
pathlist = ['../results/montepython_mcmc/w0wa_spec_pess/2022-09-09_100000_' ]
mcmc1 = mcmc(pathlist=pathlist)


## Fishers
# Fisher CF
cf1 = Fisher('../results/cosmicfish_internal/spectroscopic/pessimistic/CosmicFish_v1.0_nulcdm_internal_class-Pessimistic-own_GCsp_fishermatrix.txt')

# Fisher MP
mp1 = Fisher('../results/montepython_fisher/spectroscopic/pessimistic_HP/fisher.mat')


## Load Gaussians distributions correponding to the Fishers
gauss=GaussianND(mean=mp1.fiducials, cov=mp1.fishermatrix ,names=mp1.paramnames,labels=mp1.labels,is_inv_cov=True)
gauss2=GaussianND(mean=cf1.fiducials, cov=cf1.fishermatrix,names=cf1.paramnames,labels=cf1.labels,is_inv_cov=True)

## Cosmo Triangle
g = plots.get_subplot_plotter(settings=plot_settings)
g.triangle_plot([mcmc1.samples,gauss,gauss2],filled=[True,False,False],params=cf1.cosmonames,legend_labels=['MCMC','MP-Fisher','CF-Fisher'],contour_lws=[1,2,2])
g.export('./GCsp_Pess_cosmo.pdf')
print('Cosmo triangle finished')
del g

## Nuisance Triangle
g = plots.get_subplot_plotter(settings=plot_settings)
g.triangle_plot([mcmc1.samples,gauss,gauss2],filled=[True,False,False],params=cf1.nuisancenames,legend_labels=['MCMC','MP-Fisher','CF-Fisher'],contour_lws=[1,2,2])
g.export('./GCsp_Pess_nuisance.pdf')
print('Nuisance triangle finished')
del g

## Cosmo-Nuisance Rectangle
g = plots.get_subplot_plotter(settings=plot_settings)
g.rectangle_plot(cf1.cosmonames,cf1.nuisancenames,roots=[mcmc1.samples,gauss,gauss2],filled=[True,False,False],legend_labels=['MCMC','MP-Fisher','CF-Fisher'],contour_lws=[1,2,2])
g.export('./GCsp_Pess_cross.pdf')
print('Cross-terms rectangle finished')
del g
