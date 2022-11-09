#Initializing params and modules
import os, sys
from time import time
os.chdir(os.path.dirname(os.path.realpath(__file__)))
sys.path.append('../../../cosmicfish_reloaded')

#Importing main module
from cosmicfishpie.fishermatrix import cosmicfish
from cosmicfishpie.utilities.utils import printing as upr


###################################################
# Dictionaries for translation

obs_dict = {'GCsp' : ['GCsp'], 'WLxGCph' : ['WL', 'GCph'], 'WL' : ['WL'], 'GCph' : ['GCph']}
derivatives_default = {'GCsp' : '3PT' , 'WL' : '3PT','WLxGCph' : '3PT','GCph' : '3PT' }
paths_dict = {'WL':'lensing','WLxGCph':'photometric','GCsp':'spectroscopic'}

model = 'w0waCDM'

###################################################
def internal_runs(obs_opts,codes_list,specifications,derivatives_dictionary=derivatives_default,name=''):
    derivatives_default.update(derivatives_dictionary)
    derivatives_dict = derivatives_default.copy()
    start_time = time()
    options = {'derivatives': '3PT',
            'accuracy': 6,
            'feedback': 1,
            'outroot': 'w0wa',
            'survey_name': 'Euclid',
            'cosmo_model' : model,
            'code':'class',
            'specs_dir' : './survey_specifications/',
            'class_config_yaml' : './boltzmann_yaml_files/class/default.yaml',
            'camb_config_yaml' : './boltzmann_yaml_files/camb/default.yaml',
            'results_dir': '../../results/cosmicfish_internal/'
            }


    fiducial = {"Omegam": 0.32,
                "Omegab": 0.05,
                "h":0.67,
                "ns":0.96,
                "sigma8":0.815584,
                'mnu': 0.06,
                'Neff': 3.046,
                'w0' : -1.0,
                'wa' : 0.0
                }


    ## Main block where the code loops over different probes, codes, specifications
    for code in codes_list :
        for obs in obs_opts :
            for specifs in specifications :
                freepars = {"Omegam":0.01,
                            "Omegab":0.01 ,
                            "h":0.01,
                            "ns":0.01,
                            "sigma8":0.01,
                            'w0' : 0.01,
                            'wa' : 0.01
                            }
                options.update({
                                'derivatives' : derivatives_dict[obs],
                                'survey_name': 'Euclid-ISTF-'+specifs,
                                'outroot'  : 'w0wa'+'_internal_'+code+'-'+specifs+'-'+derivatives_dict[obs]+name,
                                'code': code,
                                'results_dir': '../../results/cosmicfish_internal/'+paths_dict[obs].lower()+'/'+specifs.lower()+'/'
                            })

                print(' *****************Internal Run: ******{coden}--{obsn}--{specn}****************'.format(coden=code, obsn=obs, specn=specifs))
                cosmoFM = cosmicfish.FisherMatrix(fiducialpars=fiducial,
                                                freepars=freepars,
                                                options=options, observables=obs_dict[obs],
                                                cosmoModel=options['cosmo_model'],
                                                surveyName=options['survey_name'])
                print(options)
                cosmoFM.compute()

    end_time =  time()
    print('Total Time taken = {}', end_time-start_time)
