import numpy as np
import os
from getdist import paramnames
from getdist import loadMCSamples
from getdist import plots
import seaborn as sns
import sys
os.chdir(os.path.dirname(os.path.realpath(__file__)))
sys.path.append('../../cosmicfish_reloaded/')
from cosmicfishpie.analysis import fisher_matrix as cfm
from cosmicfishpie.analysis import fisher_operations as cfo

################# Getdist default plot settings #############

plot_settings = plots.GetDistPlotSettings()
plot_settings.legend_fontsize = 20
plot_settings.legend_loc = 'upper right'

sns.set_theme(context='paper', style='ticks', palette='colorblind')

################## Parameters of different probes ################

cosmo_pars = ['Omegab','h','ns','sigma8','Omegam','mnu','Neff']
nuisance_wlxgcph = ['b1','b2','b3','b4','b5','b6','b7','b8','b9','b10','AIA','etaIA']
nuisance_gcsp = ['lnbgs8_1','lnbgs8_2','lnbgs8_3','lnbgs8_4','Ps_1','Ps_2','Ps_3','Ps_4']

############### Dictionaries for conversions of names and labels #################
names_mptocf={
'wa_fld' : 'wa', 'w0_fld' : 'w0','Omega_b':'Omegab', 'h' : 'h','n_s' : 'ns','sigma8' : 'sigma8','Omega_m_camb' : 'Omegam',
'bias_1' : 'b1','bias_2' : 'b2', 'bias_3' : 'b3', 'bias_4' : 'b4', 'bias_5' : 'b5','bias_6' : 'b6','bias_7' : 'b7','bias_8' : 'b8','bias_9' : 'b9','bias_10' : 'b10',
'aIA' : 'AIA', 'etaIA' :'etaIA',
'lnbsigma8_0' : 'lnbgs8_1', 'lnbsigma8_1' : 'lnbgs8_2' , 'lnbsigma8_2' : 'lnbgs8_3' , 'lnbsigma8_3' : 'lnbgs8_4',
'P_shot0' : 'Ps_1', 'P_shot1' : 'Ps_2' , 'P_shot2' : 'Ps_3', 'P_shot3' : 'Ps_4',
'N_eff_camb':'Neff' , 'm_nu_camb':'mnu'
}

labels_dict={
'wa':'w_a', 'w0':'w_0','Omegab':'\Omega_\mathrm{b}', 'h':'h','ns':'n_s','sigma8':'\sigma_8','Omegam' : '\Omega_\mathrm{m}',
'b1' : 'b_1','b2' : 'b_2', 'b3' : 'b_3', 'b4' : 'b_4', 'b5' : 'b_5','b6' : 'b_6','b7' : 'b_7','b8' : 'b_8','b9' : 'b_9','b10' : 'b_{10}',
'AIA' : 'A_{\mathrm{IA}}', 'etaIA' :'\eta_\mathrm{IA}', 'betaIA' : '\beta_\mathrm{IA}',
'lnbgs8_1' : '\ln(b_g \sigma_8)_1', 'lnbgs8_2' : '\ln(b_g \sigma_8)_2', 'lnbgs8_3' : '\ln(b_g \sigma_8)_3', 'lnbgs8_4' : '\ln(b_g \sigma_8)_4',
'Ps_1'  :  'P_{S1}', 'Ps_2'  :  'P_{S2}','Ps_3'  :  'P_{S3}','Ps_4'  :  'P_{S4}',
'Neff' : 'N_{eff}', 'mnu':'m_\nu'
}

mp_labels_dict={'wa_fld':'w_a', 'w0_fld':'w_0','Omega_b':'\Omega_\mathrm{b}', 'h':'h','n_s':'n_s','sigma8':'\sigma_8','Omega_m_camb' : '\Omega_\mathrm{m}',
'bias_1' : 'b_1','bias_2' : 'b_2', 'bias_3' : 'b_3', 'bias_4' : 'b_4', 'bias_5' : 'b_5','bias_6' : 'b_6','bias_7' : 'b_7','bias_8' : 'b_8','bias_9' : 'b_9','bias_10' : 'b_{10}',
'aIA' : 'a_{\mathrm{IA}}', 'etaIA' :'\eta_\mathrm{IA}',
'N_eff_camb':'N_{eff}', 'm_nu_camb':'m_\nu'
}



############################## Fisher matrix loading class ############################

class Fisher:
    def __init__(self,path):
        self.path = os.path.abspath(path)
        if 'spectroscopic' in self.path :
            self.nuisancenames = nuisance_gcsp.copy()
        elif 'photometric' in self.path :
            self.nuisancenames = nuisance_wlxgcph.copy()
        elif 'combined' in self.path :
            self.nuisancenames = nuisance_wlxgcph.copy() + nuisance_gcsp.copy()
        cf1 = cfm.fisher_matrix(file_name = self.path)
        cf1 = cfo.reshuffle(cf1,cosmo_pars+self.nuisancenames)
        self.fiducials = cf1.get_param_fiducial().copy()
        self.paramnames = cf1.get_param_names().copy()
        self.labels = [labels_dict[i] for i in self.paramnames ]
        self.fishermatrix = cf1.get_fisher_matrix().copy()
        self.cosmonames = cosmo_pars.copy()



############################## MCMC chains loading class ############################

class mcmc():
    def __init__(self,pathlist,probe=None):
        self.pathlist = pathlist
        self.samples = self.loadChains()

    def loadChainsParamnames(self,path):
        path = path + '.paramnames'
        params = np.genfromtxt(path,dtype=str,usecols=0)
        mpnames = [names_mptocf[i] for i in params]
        mpnames
        mplabels = [labels_dict[i] for i in mpnames]
        return paramnames.ParamNames(names=mpnames,labels=mplabels)

    def loadChains(self) :
        samples = loadMCSamples(self.pathlist[0])
        samples.setParamNames(self.loadChainsParamnames(self.pathlist[0]))
        for i in self.pathlist[1:] :
            samples1=loadMCSamples(i)
            samples1.setParamNames(self.loadChainsParamnames(i))
            samples = samples.getCombinedSamplesWithSamples(samples1)
            del samples1
        return samples
