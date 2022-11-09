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
    '../../../results/cosmicfish_external/photometric/optimistic/CosmicFish_v1.0_w0wa_external_camb-Optimistic-3PT_P2_WLGCph_fishermatrix.txt',
    '../../../results/cosmicfish_external/photometric/optimistic/CosmicFish_v1.0_w0wa_external_camb-Optimistic-3PT_P3_WLGCph_fishermatrix.txt'
              ]

labels=  [r'CF_ext_camb XCph opt P2',
          r'CF_ext_camb Xcph opt P3']

cutnames=['Omegam', 'Omegab', 'ns', 'h','sigma8','w0', 'wa', 'b1', 'b2', 'b3', 'b4', 'b5', 'b6', 'b7', 'b8', 'b9', 'b10','AIA', 'etaIA']

plotter(fish_files=fish_files,labels=labels,pars=cutnames,error_only=error_only)
