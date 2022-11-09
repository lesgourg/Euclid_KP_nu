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
    '../../../../../TestResults/FFMuF2/FFCT_Zf6/photometricExt/optimistic/CosmicFish_v1.0_w0wa_external_camb-Optimistic-3PT_P3_WLGCph_fishermatrixZf6.txt',
    '../../../../../TestResults/FFMuF2/FFCT_Zf4/photometricExt/optimistic/CosmicFish_v1.0_w0wa_external_camb-Optimistic-3PT_P3_WLGCph_fishermatrixZf4.txt',
    '../../../../../TestResults/FFMuF2/FFCT_Zf2/photometricExt/optimistic/CosmicFish_v1.0_w0wa_external_camb-Optimistic-3PT_P3_WLGCph_fishermatrixZf2.txt',
    '../../../../../TestResults/FFMuF2/FFCT_Zf1/photometricExt/optimistic/CosmicFish_v1.0_w0wa_external_camb-Optimistic-3PT_P3_WLGCph_fishermatrix.txt'
              ]

labels = [
         r'CFcambExt F6',
         r'CFcambExt F4',
         r'CFcambExt F2',
         r'CFcambExt F1'
        ]

cutnames=['Omegam', 'Omegab', 'ns', 'h','sigma8','w0', 'wa', 'b1', 'b2', 'b3', 'b4', 'b5', 'b6', 'b7', 'b8', 'b9', 'b10','AIA', 'etaIA']

compare_errors_dict={'ncol_legend':4, 'legend_title':'XCph opt', 'xticksrotation':45}# 'legend_title_fontsize':16}

plotter(fish_files=fish_files,labels=labels,pars=cutnames,error_only=error_only, compare_errors_dict=compare_errors_dict)
