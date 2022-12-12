import pandas as pd
from glob import glob
from dataframe import *
import os
os.chdir(os.path.dirname(os.path.realpath(__file__)))

spec_opt_files = glob('../../plots/spectroscopic/optimistic/*error*.txt')
photo_opt_files = glob('../../plots/photometric/optimistic/*error*.txt')
spec_pess_files = glob('../../plots/spectroscopic/pessimistic/*error*.txt')
photo_pess_files = glob('../../plots/photometric/pessimistic/*error*.txt')

## Choose the probes that needs to be printed
## Make sure to maintain the correpondance between the elements in below lists

probes = {'spec_opt':spec_opt_files,  'photo_opt':photo_opt_files, 'spec_pess' : spec_pess_files, 'photo_pess' : photo_pess_files }

p = 'spec_pess'
df1=pd.DataFrame(index=codes,columns=codes,dtype='float')
df1=df1.fillna(0)
write_data(df1,'cosmo_and_nuisance',probes[p])
save_table(df1,filename=p,save=True)
