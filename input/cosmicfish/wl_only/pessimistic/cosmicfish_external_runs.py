#Initializing params and modules
import sys,os
os.chdir(os.path.dirname(os.path.realpath(__file__)))
sys.path.append('../../')
from cosmicfish_external_master import external_runs

##################################################

## Change the keys here to change the runs accordingly
obs_opts = ['WL']                 # ['GCsp'], ['WLxGCph'],['WL'],['GCph']
codes_list = ['class']                 # ['class'], ['camb']
specifications = ['Pessimistic']        # ['Pessimistic-ql', 'Pessimistic', 'Optimistic']

## Path to the external files for class or camb 
## {codename} placeholder will be filled with either class or camb,
#   so make sure to have the same directory name.
ext_dir= '/scratch/work/casas/Euclid_WP3/external_input/default_{codename}_euclid_WP3/'

###################################################
# Derivative options
# Change the value of the derivatives dict if any other available derivative method needs to be used
# Check readme to find what derivative methods are available in cosmicfish

# derivatives_dict = {'GCsp' : 'own' , 'WL' : '3PT','WLxGCph' : '3PT','GCph' : '3PT','WL':'3PT' }

external_runs(observables=obs_opts,codes_list=codes_list,specifications=specifications)