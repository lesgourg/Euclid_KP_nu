import os, sys, argparse

os.chdir(os.path.dirname(os.path.realpath(__file__)))
sys.path.append('../../')
from plot_master import *

parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument('--error-only',action='store_true',dest='error_only',
                    help='    Plot error comparions plots only',
                    default=False)
args = parser.parse_args()
error_only = args.error_only

fish_files =  [
    #'../../../results/muRunsLinear/optimistic_HP3/fisher3.mat',
    #'../../../results/muRunsLinear/optimistic_HP/fisher9.mat',
    '../../../results/muRunsLinear/optimistic_HP50/fisher50.mat',
    '../../../results/muRunsLinear/optimistic_HP100/fisher100.mat'
              ]

labels=  [
         #r'3',
         #r'9',
          r'50',
          r'100'
          ]

cutnames=['Omegam', 'Omegab','ns', 'h','sigma8', 'w0', 'wa',  'lnbgs8_1', 'lnbgs8_2', 'lnbgs8_3', 'lnbgs8_4', 'Ps_1', 'Ps_2', 'Ps_3', 'Ps_4']

transf_labels={'\\ln(b_g \\sigma_8)_1': r'b_{g1}',
               '\\ln(b_g \\sigma_8)_2': r'b_{g2}',
               '\\ln(b_g \\sigma_8)_3': r'b_{g3}',
                '\\ln(b_g \\sigma_8)_4': r'b_{g4}'
               }

compare_errors_dict={'ncol_legend':4, 'legend_title':'GCsp opt', 'xticksrotation':45, 'transform_latex_dict': transf_labels}# 'legend_title_fontsize':16}
plotter(fish_files=fish_files,labels=labels,pars=cutnames,error_only=error_only, compare_errors_dict=compare_errors_dict)