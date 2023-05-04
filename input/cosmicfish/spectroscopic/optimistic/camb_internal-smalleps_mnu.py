#Initializing params and modules
import sys,os
os.chdir(os.path.dirname(os.path.realpath(__file__)))
sys.path.append('../../')
from cosmicfish_internal_master import internal_runs

##################################################

## Change the keys here to change the runs accordingly
obs_opts = ['GCsp']                 # ['GCsp'], ['WLxGCph'],['WL'],['GCph']
codes_list = ['camb']                 # ['class'], ['camb']
specifications = ['Optimistic']        # ['Pessimistic-ql', 'Pessimistic', 'Optimistic']
derivs_dict={'GCsp': '4PT_FWD'}         # derivative options: ['4PT_FWD', '3PT'] 
freepars_new = {"Omegam":0.01,
                "Omegab":0.01 ,
                "h":0.01,
                "ns":0.01,
                "sigma8":0.01,
                'mnu' : 0.01,  #use small eps for testing
                'Neff' : 0.01
                }

###################################################

internal_runs(obs_opts=obs_opts,codes_list=codes_list, 
              specifications=specifications, 
              derivatives_dictionary=derivs_dict,
              freepars_eps_dict=freepars_new,
              name='smalleps_mnu'
              )
