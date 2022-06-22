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

###################################################
# Derivative options
# Change the value of the derivatives dict if any other available derivative method needs to be used
# Check readme to find what derivative methods are available in cosmicfish

# derivatives_dict = {'GCsp' : 'own' , 'WL' : '3PT','WLxGCph' : '3PT','GCph' : '3PT','WL':'3PT' }

internal_runs(obs_opts=obs_opts,codes_list=codes_list,specifications=specifications)