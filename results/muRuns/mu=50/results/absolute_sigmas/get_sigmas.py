from glob import glob
import os
import numpy as np
import pandas as pd
import sys

THIS_FILE = os.path.dirname(os.path.realpath(__file__))
CF_PATH = '../../../cosmicfish_reloaded'
sys.path.append(os.path.join(THIS_FILE,CF_PATH))
from cosmicfishpie.analysis import fisher_matrix as cfm
from cosmicfishpie.analysis import fisher_operations as cfo

names_mptocf={
'wa_fld' : 'wa', 'w0_fld' : 'w0','Omega_b':'Omegab', 'h' : 'h','n_s' : 'ns','sigma8' : 'sigma8','Omega_m_camb' : 'Omegam',
'bias_1' : 'b1','bias_2' : 'b2', 'bias_3' : 'b3', 'bias_4' : 'b4', 'bias_5' : 'b5','bias_6' : 'b6','bias_7' : 'b7','bias_8' : 'b8','bias_9' : 'b9','bias_10' : 'b10',
'aIA' : 'AIA', 'etaIA' :'etaIA',
'lnbsigma8_0' : 'lnbgs8_1', 'lnbsigma8_1' : 'lnbgs8_2' , 'lnbsigma8_2' : 'lnbgs8_3' , 'lnbsigma8_3' : 'lnbgs8_4',
'P_shot0' : 'Ps_1', 'P_shot1' : 'Ps_2' , 'P_shot2' : 'Ps_3', 'P_shot3' : 'Ps_4'
}
labels_dict={
'wa':'$w_a$', 'w0':'$w_0$','Omegab':'$\Omega_\mathrm{b}$', 'h':'$h$','ns':'$n_s$','sigma8':'$\sigma_8$','Omegam' : '$\Omega_\mathrm{m}$',
'b1' : '$b_1$','b2' : '$b_2$', 'b3' : '$b_3$', 'b4' : '$b_4$', 'b5' : '$b_5$','b6' : '$b_6$','b7' : '$b_7$','b8' : '$b_8$','b9' : '$b_9$','b10' : '$b_{10}$',
'AIA' : '$A_{\mathrm{IA}}$', 'etaIA' :'$\eta_\mathrm{IA}$', 'betaIA' : '$\beta_\mathrm{IA}$',
'lnbgs8_1' : '$\ln(b_g \sigma_8)_1$', 'lnbgs8_2' : '$\ln(b_g \sigma_8)_2$', 'lnbgs8_3' : '$\ln(b_g \sigma_8)_3$', 'lnbgs8_4' : '$\ln(b_g \sigma_8)_4$',
'Ps_1'  :  '$P_{S1}$', 'Ps_2'  :  '$P_{S2}$','Ps_3'  :  '$P_{S3}$','Ps_4'  :  '$P_{S4}$'
}
cosmo_pars = ['Omegab','h','ns','sigma8','Omegam','w0','wa']
nuisance_wlxgcph = ['b1','b2','b3','b4','b5','b6','b7','b8','b9','b10','AIA','etaIA']
nuisance_gcsp = ['lnbgs8_1','lnbgs8_2','lnbgs8_3','lnbgs8_4','Ps_1','Ps_2','Ps_3','Ps_4']


def hinfo_to_table(folder,kind='pandas',save=False,only_sigmas = False,only_cosmo=True) :
    folder=os.path.normpath(folder)
    runname=os.path.basename(folder)
    file = glob(os.path.join(folder,'*.h_info'))[0]
    with open(file, 'r') as f :
        next(f)
        con = f.readlines()

    temp = [(i.strip()) for i in con]
    temp = [(i.strip()).split(':') for i in con]
    temp = [i for i in temp if i != ['']]
    x = np.array([np.fromstring(i[1].strip(), dtype = float ,  sep='\t') for i in temp])

    index = ['R-1','Best fit','mean','sigma','1-sigma-','1-sigma+ ','2-sigma-','2-sigma+','3-sigma-','3-sigma+','1-sigma >',
             '1-sigma <','2-sigma >','2-sigma <' ,'3-sigma >','3-sigma <','95% CL >','95% CL <']

    #paramnames = np.genfromtxt(glob(os.path.join(folder,'*.paramnames'))[0],dtype=str,usecols=0)
    with open(glob(os.path.join(folder,'*.bestfit'))[0], 'r') as f :
        paramnames = [names_mptocf[i.strip()] for i in f.readline().split('#')[1].strip().split(',')]


    data=pd.DataFrame(data=x,index=index,columns=paramnames)
    data=data.transpose()
    if only_cosmo :
        data = data.loc[cosmo_pars]

    if save :
        data.to_csv(runname + 'h_info')

    if kind == 'pandas' :
        if only_sigmas :
            return pd.DataFrame(data['sigma']).rename(columns={'sigma':'MontePython MCMC'})
    elif kind == 'numpy' :
        return x.T

def fisher_to_table(file_path,pars,name='CosmicFish') :
    cf1 = cfm.fisher_matrix(file_name=file_path)
    #cf1_r = cfo.reshuffle(cf1,names=pars)
    cf_m = cfo.marginalise(cf1,cosmo_pars)
    return pd.DataFrame(data=cf_m.get_confidence_bounds(),index=cosmo_pars,columns=[name])

def create_tables(paths_dict,names_dict,probe,only_cosmo = True) :
    if probe.lower() == 'wlxgcph' :
        pars = cosmo_pars+eval('nuisance_'+'wlxgcph')
    if probe.lower() == 'gcsp' :
        pars = cosmo_pars+eval('nuisance_'+'gcsp')
    if probe.lower() == 'combined' :
        pars = cosmo_pars+eval('nuisance_'+'wlxgcph')+eval('nuisance_'+'gcsp')
    if not only_cosmo:
        df = pd.DataFrame(index=pars)
    else :
        df = pd.DataFrame(index=cosmo_pars)
    for i,j in zip(paths_dict['fisher'],names_dict['fisher']) :
        cf1=fisher_to_table(i,pars=pars,name=j)
        df = df.join(cf1)
    try:
        for i in paths_dict['mcmc'] :
            mp=hinfo_to_table(i,only_sigmas=True)
            df = df.join(mp)
    except:
        pass
    df = df.rename(index=labels_dict)
    try:
        df = df.apply(lambda x : (np.vectorize(float)(np.vectorize(np.format_float_positional)(x,precision=3,unique=True,min_digits=4,trim='-'))).astype(str),axis=0)
    except:
        df = df.apply(lambda x : (np.vectorize(float)(np.vectorize(np.format_float_positional)(x,precision=5,unique=True,trim='-'))).astype(str),axis=0)
    return df

def save_table(df,filename):
    with open(filename + '.tex','w') as f :
        f.writelines(df.style.to_latex(column_format='|c|c|c|c|c|',position_float='centering'))
