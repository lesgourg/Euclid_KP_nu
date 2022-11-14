#Initializing params and modules
import os, sys
from time import time
os.chdir(os.path.dirname(os.path.realpath(__file__)))
sys.path.append('../../../cosmicfish_reloaded/')

#Importing main module
from cosmicfishpie.fishermatrix import cosmicfish

###################################################
#Dictionaries for translation

obs_dict = {'GCsp' : ['GCsp'], 'WLxGCph' : ['WL', 'GCph'], 'WL' : ['WL'], 'GCph' : ['GCph']}
derivatives_default = {'GCsp' : '3PT' , 'WLxGCph' : '3PT','GCph' : '3PT','WL':'3PT' }
paths_dict = {'WL':'wl_only','WLxGCph':'photometric','GCsp':'spectroscopic'}
ext_default = '../../../input_4_cast/output/{codename}_nulcdm_WP3_{precision}'
precision_def = ['HP']

###################################################

def external_runs(observables,codes_list,specifications,precision_list=precision_def,derivatives_dictionary=derivatives_default,name='') :

    derivatives_default.update(derivatives_dictionary)
    derivatives_dict = derivatives_default.copy()
    start_time = time()
    options = {'derivatives': '3PT',
            'accuracy': 1,
            'feedback': 1,
            'outroot': 'nulcdm',
            'survey_name': 'Euclid',
            'cosmo_model' : 'LCDM',
            'code':'external',
            'specs_dir' : './survey_specifications/',
            'class_config_yaml' : './boltzmann_yaml_files/class/default.yaml',
            'camb_config_yaml' : './boltzmann_yaml_files/camb/default.yaml',
            'results_dir': '../../results/cosmicfish_external/'
            }

    fiducial = {"Omegam": 0.3145714273,
                "Omegab": 0.0491989,
                "h":0.6737,
                "ns":0.96605,
                "sigma8":0.81,
                'mnu': 0.06,
                'Neff': 3.044,
                }

    external = {'directory':'./',
                'paramnames': ['Omegam', 'Omegab', 'h', 'ns', 'sigma8', 'mnu','Neff'],
                'folder_paramnames': ['Om', 'Ob', 'h', 'ns', 's8','Mnu','Neff'],
                'k-units' : 'h/Mpc',
                'r-units' : 'Mpc',
                'eps_values': [0.00625, 0.01, 0.0125, 0.01875, 0.02, 0.025, 0.03, 0.0375, 0.05, 0.10]}
#I dont know why i need to do this twice ask Santiago
    envkey = 'OMP_NUM_THREADS'
    print("The value of {:s} is: ".format(envkey), os.environ.get(envkey))
    os.environ[envkey] = str(8)
    os.environ[envkey] = str(8)
    print("The value of {:s} is: ".format(envkey), os.environ.get(envkey))
    for precision in precision_list :
        for code in codes_list :
            for obs in observables :
                for specifs in specifications :
                    freepars = {"Omegam":0.01,
                                "Omegab":0.01 ,
                                "h":0.01,
                                "ns":0.01,
                                "sigma8":0.01,
                                "mnu" : 0.1,
                                "Neff" : 0.01
                                } 
                    options.update({ 
                                    'derivatives' : derivatives_dict[obs],
                                    'survey_name': 'Euclid-ISTF-'+specifs,
                                    'outroot'  : 'nulcdm'+'_external_'+code+'-'+specifs+'-'+derivatives_dict[obs] + '_' + precision + name,
                                    'code': 'external',
                                    'results_dir' : '../../results/cosmicfish_external/'+paths_dict[obs].lower()+'/'+specifs.lower()+'/'
                                    })
                    external.update({'directory':ext_default.format(codename=code,precision=precision)
                                        })
                    print("Reading from dir: ", external['directory']) 
                    print(' *****************External Run: ******{coden}--{obsn}--{specn}****************'.format(coden=code, obsn=obs, specn=specifs))
                    cosmoFM = cosmicfish.FisherMatrix(fiducialpars=fiducial, 
                                                    freepars=freepars,
                                                    options=options, observables=obs_dict[obs],
                                                    cosmoModel=options['cosmo_model'], 
                                                    surveyName=options['survey_name'], extfiles=external)
                    print('********\n')


                    print(options)
                    cosmoFM.compute()

    end_time =  time()

    print('Total Time taken = {}', end_time-start_time)