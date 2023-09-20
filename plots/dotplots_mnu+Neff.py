# %%
%matplotlib inline
%load_ext autoreload
%autoreload 2-
import sys, os
from cosmicfishpie.analysis import fisher_plotting as fpp
from cosmicfishpie.analysis import fisher_matrix as fm
from cosmicfishpie.analysis import fisher_operations as fo
import seaborn as sns
snscolors=sns.color_palette("Set1")

# %%
envkey = 'OMP_NUM_THREADS'
# Set this environment variable to the number of available cores in your machine, 
# to get a fast execution of the Einstein Boltzmann Solver
print("The value of {:s} is: ".format(envkey), os.environ.get(envkey))
os.environ[envkey] = str(8)
os.environ[envkey] = str(8)
print("The value of {:s} is: ".format(envkey), os.environ.get(envkey))

# %%
transform_latex_dict ={'M_\\nu':'\sum m_\\nu', 
 'P_{S1}':r'P_{\mathrm{s}1}',
 'P_{S2}':r'P_{\mathrm{s}2}',
 'P_{S3}':r'P_{\mathrm{s}3}',
 'P_{S4}':r'P_{\mathrm{s}4}'}


for i in range(10):
    transform_latex_dict['b{}'.format(i+1)] = 'b_{{{}}}'.format(i+1)

for i in range(4):
   transform_latex_dict['\\ln(b_g \\sigma_8)_{}'.format(i+1)] = 'lbs_{{{}}}'.format(i+1)


# %%
A= fm.fisher_matrix(file_name='../../Euclid_KP_nu/results/cosmicfish_internal/photometric/optimistic/CosmicFish_v1.0_wCDM+mnu+Neff_internal_camb-Optimistic-3PT_WLGCph_fishermatrix.txt')
B= fm.fisher_matrix(file_name='../../Euclid_KP_nu/results/cosmicfish_internal/photometric/optimistic/CosmicFish_v1.0_wCDM+mnu+Neff_internal_class-Optimistic-3PT_WLGCph_fishermatrix.txt')
C= fm.fisher_matrix(file_name='../../Euclid_KP_nu/results/montepython_fisher/photo_opt_mnu+Neff/fisher.mat')
#C= fm.fisher_matrix(file_name='/home/sefa/fishers/photo_opt_wCDM/fisher.mat')

#fishers_list=[A,B,C]
fishers_list=[A,B,C]
parstoplot = ['Omegam', 'Omegab', 'h', 'ns', 'sigma8', 'mnu', 'Neff', 'b1', 'b2', 'b3', 'b4', 'b5', 'b6', 'b7', 'b8', 'b9', 'b10', 'AIA', 'etaIA']
fishers_list = [fo.reshuffle(fish,parstoplot) for fish in fishers_list]



plot_options = {'fishers_list':fishers_list, 
          'colors': snscolors,
          'fish_labels': [
                    r'$\mathtt{CF/CAMB}$',
                    r'$\mathtt{CF/CLASS}$',
                    r'$\mathtt{MP/Fisher}$',
                    #'MP Fish'
                    #,'MP Fish old'
                    ],
          'plot_pars': fishers_list[1].get_param_names()[:],
          'axis_custom_factors': {'all':3},  ## Axis limits cover 3-sigma bounds of first Fisher matrix
          'plot_method': 'Gaussian',
          'file_format': '.pdf',   ##file format for all the plots
          'outpath' : './plots/',  ## directory where to store the files, if non-existent, it will be created
          'outroot':'Photo_Opt_mnu+Neff',  ## file name root for all the plots, extra names can be added individually
          'legend_title':'Photo Opt',
          'legend_title_fontsize':30,
          'dots_legend_fontsize':30,
          'xtickfontsize' :45,
          'ytickfontsize' :35,
          'xticklabsize'  :45,
          'yticklabsize'  :35,
          'ylabelfontsize':35,
          'patches_legend_fontsize':35,
          'yrang':[-12.5,12.5],
          'transform_latex_dict':transform_latex_dict,
          'figure_title':r'$\Lambda$CDM+$\sum m_\nu$+$N_\mathrm{eff}$'
          ,'xticksrotation':45
          } 
import matplotlib.pyplot as plt
plt.style.use('../../Euclid_KP_nu/plots/plot-style.txt')

fish_plotter = fpp.fisher_plotting(**plot_options)
fish_plotter.load_gaussians()
fish_plotter.compare_errors(plot_options)

# %%
A= fm.fisher_matrix(file_name='../../Euclid_KP_nu/results/cosmicfish_internal/photometric/pessimistic/CosmicFish_v1.0_wCDM+mnu+Neff_internal_camb-Pessimistic-3PT_WLGCph_fishermatrix.txt')
B= fm.fisher_matrix(file_name='../../Euclid_KP_nu/results/cosmicfish_internal/photometric/pessimistic/CosmicFish_v1.0_wCDM+mnu+Neff_internal_class-Pessimistic-3PT_WLGCph_fishermatrix.txt')
C= fm.fisher_matrix(file_name='../../Euclid_KP_nu/results/montepython_fisher/photo_pess_mnu+Neff/fisher.mat')
#C= fm.fisher_matrix(file_name='/home/sefa/fishers/photo_opt_wCDM/fisher.mat')

#fishers_list=[A,B,C]
fishers_list=[A,B,C]
parstoplot = ['Omegam', 'Omegab', 'h', 'ns', 'sigma8', 'mnu', 'Neff', 'b1', 'b2', 'b3', 'b4', 'b5', 'b6', 'b7', 'b8', 'b9', 'b10', 'AIA', 'etaIA']
fishers_list = [fo.reshuffle(fish,parstoplot) for fish in fishers_list]

for i in range(10):
    transform_latex_dict['b{}'.format(i+1)] = 'b_{{{}}}'.format(i+1)
plot_options = {'fishers_list':fishers_list, 
          'colors': snscolors,
          'fish_labels': [
                    r'$\mathtt{CF/CAMB}$',
                    r'$\mathtt{CF/CLASS}$',
                    r'$\mathtt{MP/Fisher}$',
                    #'MP Fish'
                    #,'MP Fish old'
                    ],
          'plot_pars': fishers_list[1].get_param_names()[:],
          'axis_custom_factors': {'all':3},  ## Axis limits cover 3-sigma bounds of first Fisher matrix
          'plot_method': 'Gaussian',
          'file_format': '.pdf',   ##file format for all the plots
          'outpath' : './plots/',  ## directory where to store the files, if non-existent, it will be created
          'outroot':'Photo_Pess_mnu+Neff',  ## file name root for all the plots, extra names can be added individually
          'legend_title':'Photo Pess',
          'legend_title_fontsize':30,
          'xtickfontsize' :45,
          'ytickfontsize' :35,
          'xticklabsize'  :45,
          'yticklabsize'  :35,
          'ylabelfontsize':35,
          'patches_legend_fontsize':35,
          'dots_legend_fontsize':30,
          'yrang':[-12.5,12.5],
          'transform_latex_dict':transform_latex_dict,
          'figure_title':r'$\Lambda$CDM+$\sum m_\nu$+$N_\mathrm{eff}$'
          ,'xticksrotation':45
          } 
import matplotlib.pyplot as plt
plt.style.use('../../Euclid_KP_nu/plots/plot-style.txt')

fish_plotter = fpp.fisher_plotting(**plot_options)
fish_plotter.load_gaussians()
fish_plotter.compare_errors(plot_options)

# %%
A= fm.fisher_matrix(file_name='../../Euclid_KP_nu/results/cosmicfish_internal/spectroscopic/optimistic/CosmicFish_v1.0_wCDM+mnu+Neff_internal_camb-Optimistic-3PT_GCsp_fishermatrix.txt')
B= fm.fisher_matrix(file_name='../../Euclid_KP_nu/results/cosmicfish_internal/spectroscopic/optimistic/CosmicFish_v1.0_wCDM+mnu+Neff_internal_class-Optimistic-3PT_GCsp_fishermatrix.txt')
C= fm.fisher_matrix(file_name='../../Euclid_KP_nu/results/montepython_fisher/spec_opt_mnu+Neff/fisher.mat')

fishers_list=[A,B,C]
parstoplot= ['Omegam','Omegab','h','ns','sigma8','mnu','Neff','lnbgs8_1','lnbgs8_2','lnbgs8_3','lnbgs8_4','Ps_1','Ps_2','Ps_3','Ps_4']
fishers_list = [fo.reshuffle(fish,parstoplot) for fish in fishers_list]

plot_options = {'fishers_list':fishers_list, 
          'colors': snscolors,
          'fish_labels': [
                    r'$\mathtt{CF/CAMB}$',
                    r'$\mathtt{CF/CLASS}$',
                    r'$\mathtt{MP/Fisher}$',
                    #'MP Fish'
                    #,'MP Fish old'
                    ],
          'plot_pars': fishers_list[1].get_param_names()[:],
          'axis_custom_factors': {'all':3},  ## Axis limits cover 3-sigma bounds of first Fisher matrix
          'plot_method': 'Gaussian',
          'file_format': '.pdf',   ##file format for all the plots
          'outpath' : './plots/',  ## directory where to store the files, if non-existent, it will be created
          'outroot':'Spectro_Opt_mnu+Neff',  ## file name root for all the plots, extra names can be added individually
          'legend_title':'Spectro Opt',
          'legend_title_fontsize':30,
          'xtickfontsize' :45,
          'ytickfontsize' :35,
          'xticklabsize'  :45,
          'yticklabsize'  :35,
          'ylabelfontsize':35,
          'patches_legend_fontsize':35,
          'dots_legend_fontsize':30,
          'yrang':[-12.5,12.5],
          'transform_latex_dict':transform_latex_dict,
          'figure_title':r'$\Lambda$CDM+$\sum m_\nu$+$N_\mathrm{eff}$'
          ,'xticksrotation':46
          } 
import matplotlib.pyplot as plt
plt.style.use('../../Euclid_KP_nu/plots/plot-style.txt')

fish_plotter = fpp.fisher_plotting(**plot_options)
fish_plotter.load_gaussians()
fish_plotter.compare_errors(plot_options)

# %%
A= fm.fisher_matrix(file_name='../../Euclid_KP_nu/results/cosmicfish_internal/spectroscopic/pessimistic/CosmicFish_v1.0_wCDM+mnu+Neff_internal_camb-Pessimistic-3PT_GCsp_fishermatrix.txt')
B= fm.fisher_matrix(file_name='../../Euclid_KP_nu/results/cosmicfish_internal/spectroscopic/pessimistic/CosmicFish_v1.0_wCDM+mnu+Neff_internal_class-Pessimistic-3PT_GCsp_fishermatrix.txt')
C= fm.fisher_matrix(file_name='../../Euclid_KP_nu/results/montepython_fisher/spec_pess/fisher.mat')

fishers_list=[A,B,C]
parstoplot= ['Omegam','Omegab','h','ns','sigma8','mnu','Neff','lnbgs8_1','lnbgs8_2','lnbgs8_3','lnbgs8_4','Ps_1','Ps_2','Ps_3','Ps_4']
fishers_list = [fo.reshuffle(fish,parstoplot) for fish in fishers_list]

plot_options = {'fishers_list':fishers_list, 
          'colors': snscolors,
          'fish_labels': [
                    r'$\mathtt{CF/CAMB}$',
                    r'$\mathtt{CF/CLASS}$',
                    r'$\mathtt{MP/Fisher}$',
                    #'MP Fish'
                    #,'MP Fish old'
                    ],
          'plot_pars': fishers_list[1].get_param_names()[:],
          'axis_custom_factors': {'all':3},  ## Axis limits cover 3-sigma bounds of first Fisher matrix
          'plot_method': 'Gaussian',
          'file_format': '.pdf',   ##file format for all the plots
          'outpath' : './plots/',  ## directory where to store the files, if non-existent, it will be created
          'outroot':'Spectro_Pess_mnu+Neff',  ## file name root for all the plots, extra names can be added individually
          'legend_title':'Spectro Pess',
          'legend_title_fontsize':30,
          'xtickfontsize' :45,
          'ytickfontsize' :35,
          'xticklabsize'  :45,
          'yticklabsize'  :35,
          'ylabelfontsize':35,
          'patches_legend_fontsize':35,
          'dots_legend_fontsize':30,
          'yrang':[-12.5,12.5],
          'transform_latex_dict':transform_latex_dict,
          'figure_title':r'$\Lambda$CDM+$\sum m_\nu$+$N_\mathrm{eff}$'
          ,'xticksrotation':46
          } 
import matplotlib.pyplot as plt
plt.style.use('../../Euclid_KP_nu/plots/plot-style.txt')

fish_plotter = fpp.fisher_plotting(**plot_options)
fish_plotter.load_gaussians()
fish_plotter.compare_errors(plot_options)

# %%



