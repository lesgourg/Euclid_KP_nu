import os, sys, argparse

os.chdir(os.path.dirname(os.path.realpath(__file__)))
sys.path.append('../../')

parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('--error-only',action='store_true',dest='error_only',
                    help='    Plot error comparions plots only',
                    default=False)
args = parser.parse_args()
error_only = args.error_only

from plot_master import plotter

fish_files =  [
    '../../../results/cosmicfish_internal/photometric/pessimistic/CosmicFish_v1.0_nulcdm_internal_class-Pessimistic-3PT_WLGCph_fishermatrix.txt',
    '../../../results/cosmicfish_internal/photometric/pessimistic/CosmicFish_v1.0_nulcdm_internal_camb-Pessimistic-3PT_WLGCph_fishermatrix.txt',
    '../../../results/montepython_fisher/photometric/pessimistic_HP/fisher.mat'
    #'../../../../fisher_for_public/All_Results/pessimistic/flat/EuclidISTF_GCph_WL_XC_nulcdm_flat_pessimistic.txt'
              ]

labels = [
         r'CF class INT',
         r'CF camb INT',
         r'MP'
         #r'IST:F'
        ]

cutnames=['Omegam', 'Omegab', 'ns', 'h','sigma8', 'mnu', 'Neff', 'b1', 'b2', 'b3', 'b4', 'b5', 'b6', 'b7', 'b8', 'b9', 'b10','AIA', 'etaIA']

compare_errors_dict={'ncol_legend':4, 'legend_title':'XCph pess', 'xticksrotation':45}# 'legend_title_fontsize':16}

plotter(fish_files=fish_files,labels=labels,pars=cutnames,error_only=error_only, compare_errors_dict=compare_errors_dict)
