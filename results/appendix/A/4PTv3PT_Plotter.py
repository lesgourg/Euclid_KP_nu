from copy import deepcopy as copy
import os, sys
sys.path.append('../')

#Importing main modules
from cosmicfishpie.fishermatrix import cosmicfish
from cosmicfishpie.cosmology import cosmology as cosmo
import seaborn as sns
snscolors=sns.color_palette("Set1")
from numpy import array
import matplotlib.pyplot as plt

from cosmicfishpie.analysis import fisher_plotting as fpp
from cosmicfishpie.analysis import fisher_matrix as fm
from cosmicfishpie.analysis import fisher_operations as fo

A =fm.fisher_matrix(file_name='CosmicFish_v1.0_Camb_Spectro-Opt_3PT_epsmu_0.1_GCsp_fishermatrix.txt')
B =fm.fisher_matrix(file_name='CosmicFish_v1.0_Camb_Spectro-Opt_4PT_FWD_epsmu_0.1_GCsp_fishermatrix.txt')

plt.style.use('../../../plots/plot-style.txt')

transform_latex_dict ={'M_\\nu':'\\sum m_\\nu', 
 'P_{S1}':r'P_{\mathrm{s}1}',
 'P_{S2}':r'P_{\mathrm{s}2}',
 'P_{S3}':r'P_{\mathrm{s}3}',
 'P_{S4}':r'P_{\mathrm{s}4}'}


for i in range(10):
    transform_latex_dict['b{}'.format(i+1)] = 'b_{{{}}}'.format(i+1)

for i in range(4):
   transform_latex_dict['\\ln(b_g \\sigma_8)_{}'.format(i+1)] = 'lbs_{{{}}}'.format(i+1)

fishers_list= array([A,B])
all_pars =['Omegam','Omegab','h','ns','sigma8','Neff','w0','wa','lnbgs8_1','lnbgs8_2','lnbgs8_3','lnbgs8_4','Ps_1','Ps_2','Ps_3','Ps_4']
all_pars =['Omegam','Omegab','h','ns','sigma8','Neff','mnu','w0','wa']
fishers_list= array([fo.reshuffle(fish,all_pars) for fish in fishers_list])

plot_options = {'fishers_list':fishers_list, 
          'colors': snscolors,
          'fish_labels': [
                    r'3PT',
                    r'4PT_FWD',
                    ],
          'plot_pars': fishers_list[1].get_param_names()[:],
          'axis_custom_factors': {'all':3},  ## Axis limits cover 3-sigma bounds of first Fisher matrix
          'plot_method': 'Gaussian',
          'file_format': '.pdf',   ##file format for all the plots
          'outpath' : './plots/',  ## directory where to store the files, if non-existent, it will be created
          'outroot':'fix_bias_3PT_v_4PTFWD',  ## file name root for all the plots, extra names can be added individually
          'legend_title':r'${\tt CF/CAMB}$ Spectroscopic Optimistic',
          'legend_title_fontsize':30,
          'dots_legend_fontsize':30,
          'xtickfontsize' :45,
          'ytickfontsize' :35,
          'xticklabsize'  :45,
          'yticklabsize'  :35,
          'ylabelfontsize':35,
          'patches_legend_fontsize':35,
          'yrang':[-1,1],
          'transform_latex_dict':transform_latex_dict,
          'figure_title':r'Impact of Derivative Method'
          ,'xticksrotation':45
          } 
import matplotlib.pyplot as plt
plt.style.use('../../../plots/plot-style.txt')

fish_plotter = fpp.fisher_plotting(**plot_options)
fish_plotter.load_gaussians()
fish_plotter.compare_errors(plot_options)

fish_plotter.matrix_ratio()
