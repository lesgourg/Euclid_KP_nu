from copy import deepcopy as copy
import os, sys
sys.path.append('../')

#Importing main modules
from cosmicfishpie.fishermatrix import cosmicfish
from cosmicfishpie.cosmology import cosmology as cosmo
import seaborn as sns
snscolors=sns.color_palette("Set1")
from numpy import array

from cosmicfishpie.analysis import fisher_plotting as fpp
from cosmicfishpie.analysis import fisher_matrix as fm
from cosmicfishpie.analysis import fisher_operations as fo

A =fm.fisher_matrix(file_name='./CosmicFish_v1.0_Class_P_opt_MN_clustering_GCphWL_fishermatrix.txt')
B =fm.fisher_matrix(file_name='./CosmicFish_v1.0_Class_P_opt_MN_matter_GCphWL_fishermatrix.txt')

transform_latex_dict ={'m_\\nu':'\\sum m_\\nu', 
 'P_{S1}':r'P_{\mathrm{s}1}',
 'P_{S2}':r'P_{\mathrm{s}2}',
 'P_{S3}':r'P_{\mathrm{s}3}',
 'P_{S4}':r'P_{\mathrm{s}4}'}

for i in range(10):
    transform_latex_dict['b{}'.format(i+1)] = 'b_{{{}}}'.format(i+1)

for i in range(4):
   transform_latex_dict['\\ln(b_g \\sigma_8)_{}'.format(i+1)] = 'lbs_{{{}}}'.format(i+1)

fishers_list= array([A,B])
all_pars =['h','Omegam','Omegab','sigma8','ns','mnu','Neff','b1','b2','b3','b4','b5','b6','b7','b8','b9','b10','AIA','etaIA']
fishers_list=[fo.reshuffle(fish,all_pars) for fish in fishers_list]
pars_to_plot =['h','Omegam','Omegab','sigma8','ns','mnu','Neff','b1','b2','b3','b4','b5','b6','b7','b8','b9','b10','AIA','etaIA']

#pars_to_plot=['h','Omegam','sigma8','mnu','b4','b8']
#fishers_list=[fo.marginalise(fish,pars_to_plot) for fish in fishers_list]
plot_options = {'fishers_list':fishers_list[:], 
          'colors': snscolors,
          'fish_labels': [
                '$P_\mathrm{cc}$',
                '$P_\mathrm{mm}$',
                    ],
          'plot_pars': pars_to_plot,
          'axis_custom_factors': {'all':3},  ## Axis limits cover 3-sigma bounds of first Fisher matrix
          'plot_method': 'Gaussian',
          'file_format': '.pdf',   ##file format for all the plots
          'outpath' : './plots/',  ## directory where to store the files, if non-existent, it will be created
          'outroot':'cosmopars_Pmm_v_Pcb',  ## file name root for all the plots, extra names can be added individually
          'legend_title':r'${\tt CF/CLASS}$ Photometric Optimistic',
          'legend_title_fontsize':30,
          'dots_legend_fontsize':30,
          'xtickfontsize' :45,
          'ytickfontsize' :35,
          'xticklabsize'  :45,
          'yticklabsize'  :35,
          'ylabelfontsize':35,
          'patches_legend_fontsize':35,
          'yrang':[-40,40],
          'transform_latex_dict':transform_latex_dict,
          'figure_title':r'Impact of bias prescription'
          ,'xticksrotation':45
          }

import matplotlib.pyplot as plt
plt.style.use('../../../plots/plot-style.txt')

fish_plotter = fpp.fisher_plotting(**plot_options)
fish_plotter.load_gaussians()
fish_plotter.compare_errors(plot_options)

