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

cosmo_pars = ['h','sigma8','ns','mnu','Neff','Omegam','Omegab']
cosmo_pars = ['h','Omegam','sigma8','mnu','b4', 'b8']

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
'Neff' : 'N_\mathrm{eff}', 'mnu':r'\sum m_\nu\, \mathrm{[eV]}','logfR0':r'\log_{10}(f_{R0})',
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
                                   #'wa':[-1,1],
                                   #'w0':[-1.5,-0.5],
                                   #'N_ur':[0,None]
                                    })
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

###################photo################
all_pars =['h','Omegam','Omegab','sigma8','ns','mnu','Neff','b1','b2','b3','b4','b5','b6','b7','b8','b9','b10','AIA','etaIA']

#MCMC
pathlist = ['/home/sefa/P_MN_opt_fiting_Halofit/2023-09-14_100000_']
mcmc1 = mcmc(pathlist=pathlist)

pathlist = ['/home/sefa/P_MN_opt/2023-08-31_100000_']
mcmc2 = mcmc(pathlist=pathlist)

mean = mcmc1.samples.mean(cosmo_pars)
sigma= mcmc1.samples.std(cosmo_pars)

ticks_array=np.array([mean-2.8*sigma, mean, mean+2.8*sigma])
tick_array=[]

for i in range(len(cosmo_pars)):
    ttick = []
    for j in range(3):
        if i!=5 or ticks_array[j,i]>0.005:
            ttick.append(float('{:.2e}'.format(ticks_array[j,i])))
    tick_array.append(ttick)


plot_pars = ['h','Omegam','sigma8','mnu','b4', 'b8']

tick_array= [[0.6, 0.7], [0.31, 0.317], [0.81, 0.82], [0.1, 0.2], [1.28, 1.31], [1.46,1.5]]

g = plots.get_subplot_plotter(rc_sizes=True,subplot_size = 6,subplot_size_ratio= 1,width_inch=6)
g.triangle_plot([mcmc1.samples,mcmc2.samples],filled=[False,False],params=plot_pars,legend_labels=[r'Data ${\tt HALOFIT}$ model ${\tt HMcode}$',r'Data ${\tt HMcode}$ model ${\tt HMcode}$'], contour_lws=[1,1],colors=colors,    contour_colors=colors)
g.legend.set_title(r'Photometric optimistic')
plt.style.use('../../Euclid_KP_nu/plots/plot-style-triangle.txt')

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

g.export('./bias_plots/P_Halofit_v_P_Hmcode_bias.pdf')
print('Bias triangle finished')
del g
