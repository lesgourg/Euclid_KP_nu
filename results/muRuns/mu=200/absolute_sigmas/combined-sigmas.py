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
paths = {'fisher' : ['../../results/cosmicfish_internal/combined/pess_pess/CosmicFish_v0.9_w0wa_internal_class-combined-pess-pess_fishermatrix.txt','../../results/cosmicfish_internal/combined/pess_opt/CosmicFish_v0.9_w0wa_internal_class-combined-pess-opt_fishermatrix.txt','../../results/cosmicfish_internal/combined/opt_pess/CosmicFish_v0.9_w0wa_internal_class-combined-opt-pess_fishermatrix.txt','../../results/cosmicfish_internal/combined/opt_opt/CosmicFish_v0.9_w0wa_internal_class-combined-opt-opt_fishermatrix.txt']}

names = {'fisher' : ['pess/pess','pess/opt','opt/pess','opt/opt']}

filename = 'combined'

df = create_tables(paths_dict=paths,names_dict=names,probe='combined')
save_table(df,filename=filename)
