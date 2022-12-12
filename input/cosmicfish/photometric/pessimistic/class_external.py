#Initializing params and modules
import sys,os
os.chdir(os.path.dirname(os.path.realpath(__file__)))
sys.path.append('../../')
from cosmicfish_external_master import external_runs

##################################################

## Change the keys here to change the runs accordingly
obs_opts = ['WLxGCph']                 # ['GCsp'], ['WLxGCph'],['WL'],['GCph']
codes_list = ['class']                 # ['class'], ['camb']
specifications = ['Pessimistic']        # ['Pessimistic-ql', 'Pessimistic', 'Optimistic']
precisions = ['LP']

# Path to the external files for class or camb
# {codename} placeholder will be filled with either class or camb,
# so make sure to have the same directory name as given below.
# '../../../input_4_cast/output/{codename}_w0wa_{precision}'

###################################################
# Derivative options
# Change the value of the derivatives dict if any other available derivative method
# needs to be used
# Check readme to find what derivative methods are available in cosmicfish

# Example suppose you want to use the below derivatives, then
# dr = {'GCsp' : 'own' , 'WL' : '3PT','WLxGCph' : '3PT','GCph' : '3PT','WL':'3PT' }

# Pass the dictionary as an argument derivatives_dictionary = dr to the below given fxn


external_runs(observables=obs_opts,
codes_list=codes_list,
specifications=specifications,
precision_list=precisions)
