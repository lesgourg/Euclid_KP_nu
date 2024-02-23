import sys, os
import numpy as np
import re
from scipy.linalg import inv

from getdist import paramnames
from getdist import loadMCSamples
from getdist import plots
from getdist.gaussian_mixtures import GaussianND

import seaborn as sns
colors = [sns.color_palette('colorblind')[ii] for ii in [7,1,2]]
#colors = ['#648FFF','#DC267F','#FFB000','#785EF0','#FE6100']
#colors = ['#E69F00','#56B4E9','#009E73','#F0E442','#0072B2','#D55E00','#CC79A7']

import matplotlib.pyplot as plt

sys.path.append('../../cosmicfish_reloaded/')
from cosmicfishpie.analysis import fisher_matrix as cfm
from cosmicfishpie.analysis import fisher_operations as cfo

################# Getdist default plot settings #############

# plot_settings = plots.GetDistPlotSettings()
# plot_settings.legend_fontsize = 20
# plot_settings.legend_loc = 'upper right'
# plot_settings.axes_labelsize = 30
# plot_settings.axes_fontsize = 24

plt.style.use('../../Euclid_KP_nu/plots/plot-style-triangle.txt')
print('exit')
#############################################################

cosmo_pars = ['h','sigma8','ns','mnu','Neff','w0','wa']

##### Dictionaries for conversions of names and labels ######
names_mptocf={
'wa_fld' : 'wa', 'w0_fld' : 'w0','Omega_b':'Omegab', 'h' : 'h','n_s' : 'ns','sigma8' : 'sigma8','Omega_m_camb' : 'Omegam',
'bias_1' : 'b1','bias_2' : 'b2', 'bias_3' : 'b3', 'bias_4' : 'b4', 'bias_5' : 'b5','bias_6' : 'b6','bias_7' : 'b7','bias_8' : 'b8','bias_9' : 'b9','bias_10' : 'b10',
'aIA' : 'AIA', 'etaIA' :'etaIA',
'lnbsigma8_0' : 'lnbgs8_1', 'lnbsigma8_1' : 'lnbgs8_2' , 'lnbsigma8_2' : 'lnbgs8_3' , 'lnbsigma8_3' : 'lnbgs8_4',
'P_shot0' : 'Ps_1', 'P_shot1' : 'Ps_2' , 'P_shot2' : 'Ps_3', 'P_shot3' : 'Ps_4',
'N_eff_camb':'Neff' , 'm_nu_camb':'mnu', 'lgf_R0':'logfR0',
'sigma_p0':'sigmap_0',
'sigma_p1':'sigmap_1',
'sigma_p2':'sigmap_2',
'sigma_p3':'sigmap_3',
'sigma_v0':'sigmav_0',
'sigma_v1':'sigmav_1',
'sigma_v2':'sigmav_2',
'sigma_v3':'sigmav_3',
'N_ur':'N_ur' ,
}

labels_dict={
'wa':'w_a', 'w0':'w_0','Omegab':'\Omega_{\mathrm{b},0}', 'h':'h','ns':'n_s','sigma8':'\sigma_8','Omegam' : '\Omega_{\mathrm{m},0}',
'b1' : 'b_1','b2' : 'b_2', 'b3' : 'b_3', 'b4' : 'b_4', 'b5' : 'b_5','b6' : 'b_6','b7' : 'b_7','b8' : 'b_8','b9' : 'b_9','b10' : 'b_{10}',
'AIA' : 'A_{\mathrm{IA}}', 'etaIA' :'\eta_\mathrm{IA}', 'betaIA' : '\beta_\mathrm{IA}',
'lnbgs8_1' : '\ln(b_g \sigma_8)_1', 'lnbgs8_2' : '\ln(b_g \sigma_8)_2', 'lnbgs8_3' : '\ln(b_g \sigma_8)_3', 'lnbgs8_4' : '\ln(b_g \sigma_8)_4',
'Ps_1'  :  'P_{S1}', 'Ps_2'  :  'P_{S2}','Ps_3'  :  'P_{S3}','Ps_4'  :  'P_{S4}',
'Neff' : 'N_\mathrm{eff}', 'mnu':r'\sum m_\nu \mathrm{[eV]}','logfR0':r'\log_{10}(f_{R0})',        
        'sigmap_0': r'\sigma_{p\,0}',
        'sigmap_1': r'\sigma_{p\,1}',
        'sigmap_2': r'\sigma_{p\,2}',
        'sigmap_3': r'\sigma_{p\,3}',
        'sigmav_0': r'\sigma_{v\,0}',
        'sigmav_1': r'\sigma_{v\,1}',
        'sigmav_2': r'\sigma_{v\,2}',
        'sigmav_3': r'\sigma_{v\,3}',
        'N_ur':r'\Delta N_\mathrm{eff}'
}

############################## Fisher matrix loading class ############################

class Fisher:
    def __init__(self,path):
        self.path = os.path.abspath(path)
        cf1 = cfm.fisher_matrix(file_name = self.path)
        cf1 = cfo.reshuffle(cf1,cosmo_pars)
        self.fiducials = cf1.get_param_fiducial().copy()
        self.paramnames = cf1.get_param_names().copy()
        self.labels = [labels_dict[i] for i in self.paramnames ]
        self.fishermatrix = cf1.get_fisher_matrix().copy()
        self.cosmonames = cosmo_pars.copy()

############################## MCMC chains loading class ############################

class mcmc():
    def __init__(self,pathlist,probe=None):
        self.pathlist = pathlist
        self.samples = self.loadChains()

    def loadChainsParamnames(self,path):
        path = path + '.paramnames'
        params = np.genfromtxt(path,dtype=str,usecols=0)
        mpnames = [names_mptocf[i] for i in params]
        mpnames
        mplabels = [labels_dict[i] for i in mpnames]
        return paramnames.ParamNames(names=mpnames,labels=mplabels)

    def loadChains(self) :
        samples = loadMCSamples(self.pathlist[0],
                            settings = {'mult_bias_correction_order':1,'smooth_scale_2D':0.5, 'smooth_scale_1D':0.5}
                           )
        samples.setRanges({'mnu':[0.005,None],
                                   'wa':[-1,1],
                                   'w0':[-1.5,-0.5],
                                   'N_ur':[0,None]})
        samples.setParamNames(self.loadChainsParamnames(self.pathlist[0]))
        for i in self.pathlist[1:] :
            samples1=loadMCSamples(i)
            samples1.setParamNames(self.loadChainsParamnames(i))
            samples = samples.getCombinedSamplesWithSamples(samples1)
            del samples1
        return samples

############
# mnu+Neff #
############
cosmo_pars = ['h','Omegam','Omegab','sigma8','ns','mnu', 'Neff']

###################photo################
all_pars =['h','Omegam','Omegab','sigma8','ns','mnu','Neff','b1','b2','b3','b4','b5','b6','b7','b8','b9','b10','AIA','etaIA']

#MCMC
pathlist = ['/home/sefa/P_MN_opt_old/2023-08-21_100000_']
mcmc1 = mcmc(pathlist=pathlist)
#CF/CAMB
cf1 = cfm.fisher_matrix(file_name='../../Euclid_KP_nu/results/cosmicfish_internal/photometric/optimistic/CosmicFish_v1.0_wCDM+mnu+Neff_internal_camb-Optimistic-3PT_WLGCph_fishermatrix.txt')
cf1= cfo.reshuffle(cf1,all_pars)
cf1= cfo.marginalise(cf1, cosmo_pars)

cf2 = cfm.fisher_matrix(file_name='../../Euclid_KP_nu/results/montepython_fisher/photo_opt_mnu+Neff/fisher.mat')
cf2= cfo.reshuffle(cf2,all_pars)
cf2= cfo.marginalise(cf2, cosmo_pars)

## Load Gaussians distributions correponding to the Fishers
gauss1=GaussianND(mean=cf1.get_param_fiducial(), cov=cf1.get_fisher_matrix() ,names=cf1.get_param_names(),is_inv_cov=True)
gauss2=GaussianND(mean=cf2.get_param_fiducial(), cov=cf2.get_fisher_matrix() ,names=cf2.get_param_names(),is_inv_cov=True)

mean = mcmc1.samples.mean(cosmo_pars)
sigma= mcmc1.samples.std(cosmo_pars)

ticks_array=np.array([mean-2*sigma,mean+2*sigma])
ticks_array[:,5]=np.array([mean,mean+3.2*sigma])[:,5]
tick_array=[]

for i in range(len(cosmo_pars)):
    ttick = []
    for j in range(2):
        if i!=5 or ticks_array[j,i]>0.005:
            ttick.append(float('{:.2e}'.format(ticks_array[j,i])))
    tick_array.append(ttick)

tick_array= [[0.63, 0.73], [0.31, 0.32], [0.043, 0.055], [0.80, 0.82], [0.95, 0.98], [0.06, 0.2], [2, 4]]
## Cosmo Triangle
g = plots.get_subplot_plotter(rc_sizes=True,subplot_size = 6,subplot_size_ratio= 1,width_inch=6)
g.triangle_plot([mcmc1.samples,gauss1,gauss2],filled=[True,False,False],params=cosmo_pars,legend_labels=[r'${\tt MP/MCMC}$',r'${\tt CF/CAMB}$',r'${\tt MP/Fish}$'],contour_lws=[1,1,1],colors=colors,    contour_colors=colors)

for i in range(len(cosmo_pars)):
    for j in range(len(cosmo_pars)):
        if j > i:
            continue
        # g.subplots[i,j].minorticks_on()
        if i == j:
            g.subplots[i,i].xaxis.set_ticks(tick_array[i],labelsize=20)
            g.subplots[i,j].xaxis.set_tick_params('both',labelsize=9)
            continue
        if i > 0 and j == 0:
            g.subplots[i,j].yaxis.set_ticks(tick_array[i],labelsize=15)
            g.subplots[i,j].yaxis.set_tick_params('both',labelsize=9)
        if i== len(cosmo_pars)-1:
            g.subplots[i,j].xaxis.set_ticks(tick_array[j],labelsize=15)
            g.subplots[i,j].xaxis.set_tick_params('both',labelsize=9)

g.fig.align_ylabels(axs=None)

g.export('./Validation_Plots/P_MN_validation.pdf')
print('Cosmo triangle finished')
del g

#####################spectro######################
all_pars =['Omegam','Omegab','h','ns','sigma8','mnu','Neff','lnbgs8_1','lnbgs8_2','lnbgs8_3','lnbgs8_4','Ps_1','Ps_2','Ps_3','Ps_4']
cosmo_pars = ['h','Omegam','Omegab','sigma8','ns','mnu', 'Neff']

#MCMC
pathlist = ['/home/sefa/chains/spec_opt_mnu+Neff/2023-06-09_10000_']
mcmc1 = mcmc(pathlist=pathlist)
#CF/CAMB
cf1 = cfm.fisher_matrix(file_name='../../Euclid_KP_nu/results/cosmicfish_internal/spectroscopic/optimistic/CosmicFish_v1.0_wCDM+mnu+Neff_internal_camb-Optimistic-3PT_GCsp_fishermatrix.txt')
cf1= cfo.reshuffle(cf1,all_pars)
cf1= cfo.marginalise(cf1, cosmo_pars)

cf2 = cfm.fisher_matrix(file_name='../../Euclid_KP_nu/results/montepython_fisher/spec_opt_mnu+Neff/fisher.mat')
cf2= cfo.reshuffle(cf2,all_pars)
cf2= cfo.marginalise(cf2, cosmo_pars)

## Load Gaussians distributions correponding to the Fishers
gauss1=GaussianND(mean=cf1.get_param_fiducial(), cov=cf1.get_fisher_matrix() ,names=cf1.get_param_names(),is_inv_cov=True)
gauss2=GaussianND(mean=cf2.get_param_fiducial(), cov=cf2.get_fisher_matrix() ,names=cf2.get_param_names(),is_inv_cov=True)

mean = mcmc1.samples.mean(cosmo_pars)
sigma= mcmc1.samples.std(cosmo_pars)

ticks_array=np.array([mean-2*sigma,mean+2*sigma])
ticks_array[:,5]=np.array([mean,mean+3.2*sigma])[:,5]

tick_array=[]

for i in range(len(cosmo_pars)):
    ttick = []
    for j in range(2):
        if i!=5 or ticks_array[j,i]>0.005:
            ttick.append(float('{:.2e}'.format(ticks_array[j,i])))
    tick_array.append(ttick)

tick_array= [[0.67, 0.68], [0.31, 0.33], [0.045, 0.053], [0.79, 0.82], [0.93, 1.00], [0.12, 0.48], [2.5, 3.5]]

## Cosmo Triangle
g = plots.get_subplot_plotter(rc_sizes=True,subplot_size = 6,subplot_size_ratio= 1,width_inch=6)
g.triangle_plot([mcmc1.samples,gauss1,gauss2],filled=[True,False,False],params=cosmo_pars,legend_labels=[r'${\tt MP/MCMC}$',r'${\tt CF/CAMB}$',r'${\tt MP/Fish}$'],contour_lws=[1,1,1],colors=colors,    contour_colors=colors)

for i in range(len(cosmo_pars)):
    for j in range(len(cosmo_pars)):
        if j > i:
            continue
        # g.subplots[i,j].minorticks_on()
        if i == j:
            g.subplots[i,i].xaxis.set_ticks(tick_array[i],labelsize=20)
            g.subplots[i,j].xaxis.set_tick_params('both',labelsize=9)
            continue
        if i > 0 and j == 0:
            g.subplots[i,j].yaxis.set_ticks(tick_array[i],labelsize=15)
            g.subplots[i,j].yaxis.set_tick_params('both',labelsize=9)
        if i== len(cosmo_pars)-1:
            g.subplots[i,j].xaxis.set_ticks(tick_array[j],labelsize=15)
            g.subplots[i,j].xaxis.set_tick_params('both',labelsize=9)

g.fig.align_ylabels(axs=None)

g.export('./Validation_Plots/S_MN_validation.pdf')
print('Cosmo triangle finished')
del g


##########
# mnu+w0 #
##########
cosmo_pars = ['h','Omegam','Omegab','sigma8','ns','mnu', 'w0']

###################photo################
all_pars =['h','Omegam','Omegab','sigma8','ns','w0','mnu','b1','b2','b3','b4','b5','b6','b7','b8','b9','b10','AIA','etaIA']

#MCMC
pathlist = ['/home/sefa/P_w0M_opt_old/2023-08-21_100000_']
mcmc1 = mcmc(pathlist=pathlist)
#CF/CAMB
cf1 = cfm.fisher_matrix(file_name='../../Euclid_KP_nu/results/cosmicfish_internal/photometric/optimistic/CosmicFish_v1.0_wCDM+mnu+Neff_internal_camb-Optimistic-3PT_WLGCph_fishermatrix.txt')
cf1= cfo.reshuffle(cf1,all_pars)
cf1= cfo.marginalise(cf1, cosmo_pars)

cf2 = cfm.fisher_matrix(file_name='../../Euclid_KP_nu/results/montepython_fisher/photo_opt_mnu+w0/fisher.mat')
cf2= cfo.reshuffle(cf2,all_pars)
cf2= cfo.marginalise(cf2, cosmo_pars)

## Load Gaussians distributions correponding to the Fishers
gauss1=GaussianND(mean=cf1.get_param_fiducial(), cov=cf1.get_fisher_matrix() ,names=cf1.get_param_names(),is_inv_cov=True)
gauss2=GaussianND(mean=cf2.get_param_fiducial(), cov=cf2.get_fisher_matrix() ,names=cf2.get_param_names(),is_inv_cov=True)

mean = mcmc1.samples.mean(cosmo_pars)
sigma= mcmc1.samples.std(cosmo_pars)

ticks_array=np.array([mean-2*sigma,mean+2*sigma])
ticks_array[:,5]=[0.06,0.24]
ticks_array[:,3]=[0.8,0.818]
ticks_array[:,6]=[-1.03,-0.97]
ticks_array[:,1]=[0.311,0.319]


tick_array=[]

for i in range(len(cosmo_pars)):
    ttick = []
    for j in range(2):
        if i!=5 or ticks_array[j,i]>0.005:
            ttick.append(float('{:.2e}'.format(ticks_array[j,i])))
    tick_array.append(ttick)

tick_array= [[0.63, 0.72], [0.31, 0.32], [0.043, 0.057], [0.8, 0.82], [0.94, 0.98], [0.06, 0.24], [-1.03, -0.97]]

## Cosmo Triangle
g = plots.get_subplot_plotter(rc_sizes=True,subplot_size = 6,subplot_size_ratio= 1,width_inch=6)
g.triangle_plot([mcmc1.samples,gauss1,gauss2],filled=[True,False,False],params=cosmo_pars,legend_labels=[r'${\tt MP/MCMC}$',r'${\tt CF/CAMB}$',r'${\tt MP/Fish}$'],contour_lws=[1,1,1],colors=colors,    contour_colors=colors)

for i in range(len(cosmo_pars)):
    for j in range(len(cosmo_pars)):
        if j > i:
            continue
        # g.subplots[i,j].minorticks_on()
        if i == j:
            g.subplots[i,i].xaxis.set_ticks(tick_array[i],labelsize=20)
            g.subplots[i,j].xaxis.set_tick_params('both',labelsize=9)
            continue
        if i > 0 and j == 0:
            g.subplots[i,j].yaxis.set_ticks(tick_array[i],labelsize=15)
            g.subplots[i,j].yaxis.set_tick_params('both',labelsize=9)
        if i== len(cosmo_pars)-1:
            g.subplots[i,j].xaxis.set_ticks(tick_array[j],labelsize=15)
            g.subplots[i,j].xaxis.set_tick_params('both',labelsize=9)

g.fig.align_ylabels(axs=None)

g.export('./Validation_Plots/P_w0M_validation.pdf')
print('Cosmo triangle finished')
del g

#####################spectro######################
all_pars =['Omegam','Omegab','h','ns','sigma8','w0','mnu','lnbgs8_1','lnbgs8_2','lnbgs8_3','lnbgs8_4','Ps_1','Ps_2','Ps_3','Ps_4']
cosmo_pars = ['h','Omegam','Omegab','sigma8','ns','mnu', 'w0']

#MCMC
pathlist = ['/home/sefa/chains/spec_opt_mnu+w0/2023-06-09_10000_']
mcmc1 = mcmc(pathlist=pathlist)
#CF/CAMB
cf1 = cfm.fisher_matrix(file_name='../../Euclid_KP_nu/results/cosmicfish_internal/spectroscopic/optimistic/CosmicFish_v1.0_wCDM+mnu+Neff_internal_camb-Optimistic-3PT_GCsp_fishermatrix.txt')
cf1= cfo.reshuffle(cf1,all_pars)
cf1= cfo.marginalise(cf1, cosmo_pars)

cf2 = cfm.fisher_matrix(file_name='../../Euclid_KP_nu/results/montepython_fisher/spec_opt_mnu+w0/fisher.mat')
cf2= cfo.reshuffle(cf2,all_pars)
cf2= cfo.marginalise(cf2, cosmo_pars)

## Load Gaussians distributions correponding to the Fishers
gauss1=GaussianND(mean=cf1.get_param_fiducial(), cov=cf1.get_fisher_matrix() ,names=cf1.get_param_names(),is_inv_cov=True)
gauss2=GaussianND(mean=cf2.get_param_fiducial(), cov=cf2.get_fisher_matrix() ,names=cf2.get_param_names(),is_inv_cov=True)

mean = mcmc1.samples.mean(cosmo_pars)
sigma= mcmc1.samples.std(cosmo_pars)

ticks_array=np.array([mean-2*sigma,mean+2*sigma])
ticks_array[:,5]=[0.12,0.36]
tick_array=[]

for i in range(len(cosmo_pars)):
    ttick = []
    for j in range(2):
        if i!=5 or ticks_array[j,i]>0.005:
            ttick.append(float('{:.2e}'.format(ticks_array[j,i])))
    tick_array.append(ttick)

tick_array= [[0.67, 0.68], [0.3, 0.33], [0.047, 0.052], [0.79, 0.83], [0.95, 1.], [0.12, 0.48], [-1.07, -0.93]]

## Cosmo Triangle
g = plots.get_subplot_plotter(rc_sizes=True,subplot_size = 6,subplot_size_ratio= 1,width_inch=6)
g.triangle_plot([mcmc1.samples,gauss1,gauss2],filled=[True,False,False],params=cosmo_pars,legend_labels=[r'${\tt MP/MCMC}$',r'${\tt CF/CAMB}$',r'${\tt MP/Fish}$'],contour_lws=[1,1,1],colors=colors,    contour_colors=colors)

for i in range(len(cosmo_pars)):
    for j in range(len(cosmo_pars)):
        if j > i:
            continue
        # g.subplots[i,j].minorticks_on()
        if i == j:
            g.subplots[i,i].xaxis.set_ticks(tick_array[i],labelsize=20)
            g.subplots[i,j].xaxis.set_tick_params('both',labelsize=9)
            continue
        if i > 0 and j == 0:
            g.subplots[i,j].yaxis.set_ticks(tick_array[i],labelsize=15)
            g.subplots[i,j].yaxis.set_tick_params('both',labelsize=9)
        if i== len(cosmo_pars)-1:
            g.subplots[i,j].xaxis.set_ticks(tick_array[j],labelsize=15)
            g.subplots[i,j].xaxis.set_tick_params('both',labelsize=9)

g.fig.align_ylabels(axs=None)

g.export('./Validation_Plots/S_w0M_validation.pdf')
print('Cosmo triangle finished')
del g

#########
# w0+wa #
#########
cosmo_pars = ['h','Omegam','Omegab','sigma8','ns','w0', 'wa']

###################photo################
all_pars =['h','Omegam','Omegab','sigma8','ns','w0','wa','b1','b2','b3','b4','b5','b6','b7','b8','b9','b10','AIA','etaIA']


#MCMC
pathlist = ['/home/sefa/P_w0wa_opt_old/2023-08-21_100000_']
mcmc1 = mcmc(pathlist=pathlist)
#CF/CAMB
cf1 = cfm.fisher_matrix(file_name='../../Euclid_KP_nu/results/cosmicfish_internal/photometric/optimistic/CosmicFish_v1.0_wCDM+mnu+Neff_internal_camb-Optimistic-3PT_WLGCph_fishermatrix.txt')
cf1= cfo.reshuffle(cf1,all_pars)
cf1= cfo.marginalise(cf1, cosmo_pars)

cf2 = cfm.fisher_matrix(file_name='../../Euclid_KP_nu/results/montepython_fisher/photo_opt_wCDM/fisher.mat')
cf2= cfo.reshuffle(cf2,all_pars)
cf2= cfo.marginalise(cf2, cosmo_pars)

## Load Gaussians distributions correponding to the Fishers
gauss1=GaussianND(mean=cf1.get_param_fiducial(), cov=cf1.get_fisher_matrix() ,names=cf1.get_param_names(),is_inv_cov=True)
gauss2=GaussianND(mean=cf2.get_param_fiducial(), cov=cf2.get_fisher_matrix() ,names=cf2.get_param_names(),is_inv_cov=True)

mean = mcmc1.samples.mean(cosmo_pars)
sigma= mcmc1.samples.std(cosmo_pars)

ticks_array=np.array([mean-2*sigma,mean+2*sigma])
tick_array=[]

for i in range(len(cosmo_pars)):
    ttick = []
    for j in range(2):
        if True:
            ttick.append(float('{:.2e}'.format(ticks_array[j,i])))
    tick_array.append(ttick)
tick_array= [[0.63, 0.73], [0.31, 0.32], [0.044, 0.055], [0.805, 0.815], [0.95, 1.], [-1.06, -0.94], [-0.25, 0.25]]

## Cosmo Triangle
g = plots.get_subplot_plotter(rc_sizes=True,subplot_size = 6,subplot_size_ratio= 1,width_inch=6)
g.triangle_plot([mcmc1.samples,gauss1,gauss2],filled=[True,False,False],params=cosmo_pars,legend_labels=[r'${\tt MP/MCMC}$',r'${\tt CF/CAMB}$',r'${\tt MP/Fish}$'],contour_lws=[1,1,1],colors=colors,    contour_colors=colors)

for i in range(len(cosmo_pars)):
    for j in range(len(cosmo_pars)):
        if j > i:
            continue
        # g.subplots[i,j].minorticks_on()
        if i == j:
            g.subplots[i,i].xaxis.set_ticks(tick_array[i],labelsize=20)
            g.subplots[i,j].xaxis.set_tick_params('both',labelsize=9)
            continue
        if i > 0 and j == 0:
            g.subplots[i,j].yaxis.set_ticks(tick_array[i],labelsize=15)
            g.subplots[i,j].yaxis.set_tick_params('both',labelsize=9)
        if i== len(cosmo_pars)-1:
            g.subplots[i,j].xaxis.set_ticks(tick_array[j],labelsize=15)
            g.subplots[i,j].xaxis.set_tick_params('both',labelsize=9)

g.fig.align_ylabels(axs=None)

g.export('./Validation_Plots/P_w0wa_validation.pdf')
print('Cosmo triangle finished')
del g

#####################spectro######################
all_pars =['Omegam','Omegab','h','ns','sigma8','w0','wa','lnbgs8_1','lnbgs8_2','lnbgs8_3','lnbgs8_4','Ps_1','Ps_2','Ps_3','Ps_4']
cosmo_pars = ['h','Omegam','Omegab','sigma8','ns','w0', 'wa']

#MCMC
pathlist = ['/home/sefa/chains/spec_opt_wCDM/2023-06-09_10000_']
mcmc1 = mcmc(pathlist=pathlist)
#CF/CAMB
cf1 = cfm.fisher_matrix(file_name='../../Euclid_KP_nu/results/cosmicfish_internal/spectroscopic/optimistic/CosmicFish_v1.0_wCDM+mnu+Neff_internal_camb-Optimistic-3PT_GCsp_fishermatrix.txt')
cf1= cfo.reshuffle(cf1,all_pars)
cf1= cfo.marginalise(cf1, cosmo_pars)

cf2 = cfm.fisher_matrix(file_name='../../Euclid_KP_nu/results/montepython_fisher/spec_opt/fisher.mat')
cf2= cfo.reshuffle(cf2,all_pars)
cf2= cfo.marginalise(cf2, cosmo_pars)

## Load Gaussians distributions correponding to the Fishers
gauss1=GaussianND(mean=cf1.get_param_fiducial(), cov=cf1.get_fisher_matrix() ,names=cf1.get_param_names(),is_inv_cov=True)
gauss2=GaussianND(mean=cf2.get_param_fiducial(), cov=cf2.get_fisher_matrix() ,names=cf2.get_param_names(),is_inv_cov=True)

mean = mcmc1.samples.mean(cosmo_pars)
sigma= mcmc1.samples.std(cosmo_pars)

ticks_array=np.array([mean-2*sigma,mean+2*sigma])
tick_array=[]

for i in range(len(cosmo_pars)):
    ttick = []
    for j in range(2):
        if True:
            ttick.append(float('{:.2e}'.format(ticks_array[j,i])))
    tick_array.append(ttick)

tick_array= [[0.67, 0.68], [0.30, 0.34], [0.045, 0.053], [0.79, 0.83], [0.94, 1.], [-1.2, -0.8], [-0.6, 0.6]]

## Cosmo Triangle
g = plots.get_subplot_plotter(rc_sizes=True,subplot_size = 6,subplot_size_ratio= 1,width_inch=6)
g.triangle_plot([mcmc1.samples,gauss1,gauss2],filled=[True,False,False],params=cosmo_pars,legend_labels=[r'${\tt MP/MCMC}$',r'${\tt CF/CAMB}$',r'${\tt MP/Fish}$'],contour_lws=[1,1,1],colors=colors,    contour_colors=colors)

for i in range(len(cosmo_pars)):
    for j in range(len(cosmo_pars)):
        if j > i:
            continue
        # g.subplots[i,j].minorticks_on()
        if i == j:
            g.subplots[i,i].xaxis.set_ticks(tick_array[i],labelsize=20)
            g.subplots[i,j].xaxis.set_tick_params('both',labelsize=9)
            continue
        if i > 0 and j == 0:
            g.subplots[i,j].yaxis.set_ticks(tick_array[i],labelsize=15)
            g.subplots[i,j].yaxis.set_tick_params('both',labelsize=9)
        if i== len(cosmo_pars)-1:
            g.subplots[i,j].xaxis.set_ticks(tick_array[j],labelsize=15)
            g.subplots[i,j].xaxis.set_tick_params('both',labelsize=9)

g.fig.align_ylabels(axs=None)

g.export('./Validation_Plots/S_w0wa_validation.pdf')
print('Cosmo triangle finished')
del g