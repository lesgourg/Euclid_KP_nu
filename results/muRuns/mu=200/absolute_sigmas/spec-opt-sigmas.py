from get_sigmas import *
os.chdir(os.path.dirname(os.path.realpath(__file__)))

## Column ordering in the Table follows the pattern : Fisher1, Fisher2, .... MCMC
## Therefore, pass the Fishers in the required order inside the dictonaries paths, and names
## paths : contains the paths to the Fishers
## names : Column heading for the latex table
## filename : name of the file (do not add .tex as it will be automatically added)
## create_tables : a function that returns a pandas dataframe and it can be further used if needed
## save_table : converts the dataframe into latex table

################################################# Spec Opt #################################################
paths = {'mcmc' : ['../../results/montepython_mcmc/w0wa_spec_opt'], 'fisher' : ['../../external_fishers/soapfish/km030/SOAPFish_C2-NL2-km030-w0wa_cdm-1p0E-02.dat','../../results/cosmicfish_internal/spectroscopic/optimistic/CosmicFish_v0.9_w0wa_internal_class-Optimistic-own_GCsp_fishermatrix.txt','../../results/montepython_fisher/spectroscopic/optimistic_HP/fisher.mat']}

names = {'fisher' : ['SoapFish','CosmicFish','MontePython Fisher']}

filename = 'spec_opt'

df = create_tables(paths_dict=paths,names_dict=names,probe='GCsp')
save_table(df,filename=filename)
