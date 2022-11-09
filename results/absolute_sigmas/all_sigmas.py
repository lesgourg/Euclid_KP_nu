from get_sigmas import *
os.chdir(os.path.dirname(os.path.realpath(__file__)))

## Column ordering in the Table follows the pattern : Fisher1, Fisher2, .... MCMC
## Therefore, pass the Fishers in the required order inside the dictonaries paths, and names
## paths : contains the paths to the Fishers
## names : Column heading for the latex table
## filename : name of the file (do not add .tex as it will be automatically added)
## create_tables : a function that returns a pandas dataframe and it can be further used if needed
## save_table : converts the dataframe into latex table


################################################ Photo Pess #####################################
paths = {'mcmc' : ['../../results/montepython_mcmc/w0wa_photo_pess'], 'fisher' : ['../../../fisher_for_public/All_Results/pessimistic/flat/EuclidISTF_GCph_WL_XC_w0wa_flat_pessimistic.txt','../../results/cosmicfish_internal/photometric/pessimistic/CosmicFish_v0.9_w0wa_internal_class-Pessimistic-3PT_WLGCph_fishermatrix.txt','../../results/montepython_fisher/photometric/pessimistic/fisher.mat']}

names = {'fisher' : ['IST:Fisher','CosmicFish','MontePython Fisher']}

filename = 'photo_pess'

df = create_tables(paths_dict=paths,names_dict=names,probe='WLxGCph')
save_table(df,filename=filename)

################################################ Photo Opt#########################################
paths = {'mcmc' : ['../../results/montepython_mcmc/w0wa_photo_opt'], 'fisher' : ['../../../fisher_for_public/All_Results/optimistic/flat/EuclidISTF_GCph_WL_XC_w0wa_flat_optimistic.txt','../../results/cosmicfish_internal/photometric/optimistic/CosmicFish_v0.9_w0wa_internal_class-Optimistic-3PT_WLGCph_fishermatrix.txt','../../results/montepython_fisher/photometric/optimistic/fisher.mat']}

names = {'fisher' : ['IST:Fisher','CosmicFish','MontePython Fisher']}

filename = 'photo_opt'

df = create_tables(paths_dict=paths,names_dict=names,probe='WLxGCph')
save_table(df,filename=filename)


################################################# Spec Pess ###############################################
paths = {'mcmc' : ['../../results/montepython_mcmc/w0wa_spec_opt'], 'fisher' : ['../../external_fishers/soapfish/km025/SOAPFish_C2-NL2-km025-w0wa-flat-full-1p0E-02-last.dat','../../results/cosmicfish_internal/spectroscopic/pessimistic/CosmicFish_v0.9_w0wa_internal_class-Pessimistic-own_GCsp_fishermatrix.txt','../../results/montepython_fisher/spectroscopic/pessimistic/fisher.mat']}

names = {'fisher' : ['SoapFish','CosmicFish','MontePython Fisher']}

filename = 'spec_pess'

df = create_tables(paths_dict=paths,names_dict=names,probe='GCsp')
save_table(df,filename=filename)


################################################# Spec Opt #################################################
paths = {'mcmc' : ['../../results/montepython_mcmc/w0wa_spec_opt'], 'fisher' : ['../../external_fishers/soapfish/km030/SOAPFish_C2-NL2-km030-w0wa-flat-full-1p0E-02-last.dat','../../results/cosmicfish_internal/spectroscopic/optimistic/CosmicFish_v0.9_w0wa_internal_class-Optimistic-own_GCsp_fishermatrix.txt','../../results/montepython_fisher/spectroscopic/optimistic/fisher.mat']}

names = {'fisher' : ['SoapFish','CosmicFish','MontePython Fisher']}

filename = 'spec_opt'

df = create_tables(paths_dict=paths,names_dict=names,probe='GCsp')
save_table(df,filename=filename)
