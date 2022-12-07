from get_sigmas import *
os.chdir(os.path.dirname(os.path.realpath(__file__)))

## Column ordering in the Table follows the pattern : Fisher1, Fisher2, .... MCMC
## Therefore, pass the Fishers in the required order inside the dictonaries paths, and names
## paths : contains the paths to the Fishers
## names : Column heading for the latex table
## filename : name of the file (do not add .tex as it will be automatically added)
## create_tables : a function that returns a pandas dataframe and it can be further used if needed
## save_table : converts the dataframe into latex table


################################################# Spec Pess ###############################################
paths = {'mcmc' : ['../../results/montepython_mcmc/nulcdm_spec_opt'], 'fisher' : ['../../results/cosmicfish_internal/spectroscopic/pessimistic/CosmicFish_v1.0_nulcdm_internal_class-Pessimistic-3PT_GCsp_fishermatrix.txt','../../results/montepython_fisher/spectroscopic/pessimistic_HP/fisher.mat']}

names = {'fisher' : ['CosmicFish','MontePython Fisher']}

filename = 'spec_pess'

df = create_tables(paths_dict=paths,names_dict=names,probe='GCsp')
save_table(df,filename=filename)
