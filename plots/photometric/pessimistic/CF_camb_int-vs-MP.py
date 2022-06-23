import errno
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
    '../../../results/cosmicfish_internal/photometric/pessimistic/CosmicFish_v0.9_varying_mnu__internal_camb-Pessimistic-3PT_WLGCph_fishermatrix.txt',
    '../../../results/montepython_fisher/photometric/pessimistic/fisher.mat'
              ]

labels = [r'CF_int_camb XCph pess',
          r'MP XCph pess' ]

cutnames=['Neff', 'mnu', 'Omegam', 'Omegab', 'ns', 'h','sigma8', 'b1', 'b2', 'b3', 'b4', 'b5', 'b6', 'b7', 'b8', 'b9', 'b10','AIA', 'etaIA']

plotter(fish_files=fish_files,labels=labels,pars=cutnames,error_only=error_only) 
