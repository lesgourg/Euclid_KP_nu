#You need to have CAMB and CLASS installed in your environment
import camb as camb
from classy import Class

from copy import deepcopy
from scipy.interpolate import RectBivariateSpline, UnivariateSpline
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

plt.style.use('../../Euclid_KP_nu/plots/plot-style.txt')

#Set Accuracy Settings for Class and Camb
Class_dict_DP = {
  'k_per_decade_for_bao': 70,
  'k_per_decade_for_pk': 10,
  'l_max_g' : 20,
  'l_max_pol_g' : 15,
  'radiation_streaming_approximation' : 2,
  'radiation_streaming_trigger_tau_over_tau_k' : 240.,
  'radiation_streaming_trigger_tau_c_over_tau' : 100.,
  'tol_ncdm_synchronous' : 1.e-5,
  'ncdm_fluid_trigger_tau_over_tau_k' : 41.,
  'l_max_ncdm' : 22,
  'background_Nloga' : 6000,
  'thermo_Nz_log' : 20000,
  'thermo_Nz_lin' : 40000,
  'tol_perturbations_integration' : 1.e-6,
  'P_k_max_1/Mpc': 50,
  'output' : 'mPk,mTk',
  'non linear' : 'hmcode2020',
  'z_max_pk' : 2.0,
  'N_ncdm': 1,
  'T_cmb' : 2.7255,
  'nonlinear_min_k_max' : 80.,
  'hmcode_tol_sigma': 1.e-8,
  'k_pivot' : 0.05,
}

Class_dict_HP = {
  'k_per_decade_for_bao': 70,
  'k_per_decade_for_pk': 10,
  'l_max_g' : 20,
  'l_max_pol_g' : 15,
  'radiation_streaming_approximation' : 2,
  'radiation_streaming_trigger_tau_over_tau_k' : 240.,
  'radiation_streaming_trigger_tau_c_over_tau' : 100.,
  'tol_ncdm_synchronous' : 1.e-5,
  'ncdm_fluid_trigger_tau_over_tau_k' : 90.,
  'l_max_ncdm' : 25,
  'background_Nloga' : 6000,
  'thermo_Nz_log' : 20000,
  'thermo_Nz_lin' : 40000,
  'tol_perturbations_integration' : 1.e-6,
  'P_k_max_1/Mpc': 50,
  'output' : 'mPk,mTk',
  'non linear' : 'hmcode2020',
  'z_max_pk' : 2.0,
  'N_ncdm': 1,
  'T_cmb' : 2.7255,
  'nonlinear_min_k_max' : 80.,
  'hmcode_tol_sigma': 1.e-8,
  'k_pivot' : 0.05,
}

Class_dict_HP_lowk = {
  'k_per_decade_for_bao': 70,
  'k_per_decade_for_pk': 10,
  'l_max_g' : 20,
  'l_max_pol_g' : 15,
  'radiation_streaming_approximation' : 2,
  'radiation_streaming_trigger_tau_over_tau_k' : 240.,
  'radiation_streaming_trigger_tau_c_over_tau' : 100.,
  'tol_ncdm_synchronous' : 1.e-5,
  'ncdm_fluid_approximation': 3,
  'evolver' : 0,
  'l_max_ncdm' : 25,
  'background_Nloga' : 6000,
  'thermo_Nz_log' : 20000,
  'thermo_Nz_lin' : 40000,
  'tol_perturbations_integration' : 1.e-6,
  'P_k_max_1/Mpc': 15,
  'output' : 'mPk,mTk',
  'non linear' : 'hmcode2020',
  'z_max_pk' : 2.0,
  'N_ncdm': 1,
  'T_cmb' : 2.7255,
  'nonlinear_min_k_max' : 80.,
  'hmcode_tol_sigma': 1.e-8,
  'k_pivot' : 0.05,
}

Class_dict_UHP = {
  'k_per_decade_for_bao': 70,
  'k_per_decade_for_pk': 10,
  'l_max_g' : 20,
  'l_max_pol_g' : 15,
  'radiation_streaming_approximation' : 2,
  'radiation_streaming_trigger_tau_over_tau_k' : 240.,
  'radiation_streaming_trigger_tau_c_over_tau' : 100.,
  'tol_ncdm_synchronous' : 1.e-5,
  'l_max_ncdm' : 25,
  'ncdm_fluid_approximation': 3,
  'background_Nloga' : 6000,
  'thermo_Nz_log' : 20000,
  'thermo_Nz_lin' : 40000,
  'tol_perturbations_integration' : 1.e-6,
  'evolver' : 0,
  'P_k_max_1/Mpc': 15,
  'output' : 'mPk,mTk',
  'z_max_pk' : 2.0,
  'N_ncdm': 1,
  'T_cmb' : 2.7255,
  'k_pivot' : 0.05,
}

Camb_dict_DP={
  'kmax' : 50,
  'k_per_logint' : 50,
  'AccuracyBoost' : 2,
  'lAccuracyBoost' : 2,
  'lSampleBoost' : 1,
  'accurate_massive_neutrino_transfers' : True,
  'halofit_version' : 9,
  'dark_energy_model' : 'ppf',
  'num_nu_massive' : 1,
  'TCMB': 2.7255,
  'WantTensors' : False,
  'WantVectors' : False,
  'Reion.Reionization' : False,
  'WantCls' : False,
  'Want_CMB' : False,
  'Want_CMB_lensing': False,
  'Want_cl_2D_array' : False,
  'DoLateRadTruncation': True,
  'DoLensing': False,
}

Camb_dict_HP={
  'kmax' : 50,
  'k_per_logint' : 50,
  'AccuracyBoost' : 3,
  'lAccuracyBoost' : 3,
  'lSampleBoost' : 1,
  'accurate_massive_neutrino_transfers' : True,
  'halofit_version' : 9,
  'dark_energy_model' : 'ppf',
  'num_nu_massive' : 1,
  'TCMB': 2.7255,
  'WantTensors' : False,
  'WantVectors' : False,
  'Reion.Reionization' : False,
  'WantCls' : False,
  'Want_CMB' : False,
  'Want_CMB_lensing': False,
  'Want_cl_2D_array' : False,
  'DoLateRadTruncation': True,
  'DoLensing': False,
}


fiducial = {"Omegam": 0.3145714273,
            "Omegab": 0.0491989,
            "h":0.6737,
            "ns":0.96605,
            "As":2.1e-9,
            'mnu': 0.06,
            'Neff': 3.044,
            'w0': -1,
            'wa':  0
            }

#Basis Conversion for Class
def basechange_class(cosmopars,shareDeltaNeff=False):
    # transforms cosmopars into cosmopars that can be read by CLASS
 
    classpars = deepcopy(cosmopars)
    if 'h'      in classpars: 
        classpars['h']    = classpars.pop('h') 
        h = classpars['h']
    if 'H0'     in classpars: 
        classpars['H0']   = classpars.pop('H0')
        h = classpars['H0'] / 100.
    if 'ns'     in classpars: classpars['n_s'] = classpars.pop('ns')

    if 'w0'     in classpars: 
        classpars['w0_fld']    = classpars.pop('w0')
        classpars['Omega_Lambda'] = 0
    if 'wa'     in classpars: classpars['wa_fld'] = classpars.pop('wa')
    if 'As'     in classpars: classpars['A_s'] = classpars.pop('As')   

    Neff = classpars.pop('Neff')
    fidNeff = 3.044

    if shareDeltaNeff:
        classpars['N_ur'] = 2./3.*Neff #This version does not have the discontinuity at Nur = 1.99
        g_factor = Neff/3.
    else:
        classpars['N_ur'] = Neff - fidNeff/3.
        g_factor = fidNeff/3.
    neutrino_mass_fac  = 94.07
    if 'mnu' in classpars: 
        classpars['T_ncdm'] = (4./11.)**(1./3.) * g_factor**(1./4.)
        classpars['Omega_ncdm'] = classpars['mnu'] * g_factor**(0.75) / neutrino_mass_fac / h**2
        classpars.pop('mnu')
    elif 'Omeganu' in classpars:
        classpars['Omega_ncdm'] = classpars.pop('Omeganu')
    if 'Omegab' in classpars: classpars['Omega_b'] = classpars.pop('Omegab')
    if 'Omegam' in classpars:
         classpars['Omega_cdm'] = classpars.pop('Omegam') - classpars['Omega_b'] - classpars['Omega_ncdm']

    return classpars

#Basis Conversion for Camb
def basechange_camb(cosmopars,shareDeltaNeff=False):
    # transforms cosmopars into cosmopars that can be read by CAMB
    cambpars = deepcopy(cosmopars)
    if 'h'      in cambpars: cambpars['H0']    = cambpars.pop('h')*100
    if 'Omegab' in cambpars: cambpars['ombh2'] = cambpars.pop('Omegab')*(cambpars['H0']/100)**2
    if 'Omegak' in cambpars: cambpars['omk']   = cambpars.pop('Omegak')
    if 'w0'     in cambpars: 
        cambpars['w']     = cambpars.pop('w0')

    cambpars['share_delta_neff']=shareDeltaNeff
    fidNeff = 3.044
    if 'Neff' in cambpars:
        Neff = cambpars.pop('Neff')
        if shareDeltaNeff:
            cambpars['num_nu_massless'] = Neff - cambpars['num_nu_massive']
        else:
            cambpars['num_nu_massless'] = Neff - fidNeff/3
    
    else:
        Neff = cambpars['num_nu_massive']+cambpars['num_nu_massless']
    
    if shareDeltaNeff:
        g_factor = Neff/3
    else:
        g_factor = fidNeff/3
    
    neutrino_mass_fac  = 94.07
    h2 = (cambpars['H0']/100)**2
    if 'mnu' in cambpars: 
        Onu                = cambpars['mnu']/neutrino_mass_fac*(g_factor)** 0.75/h2
        onuh2  = Onu*h2
        cambpars['omnuh2'] = onuh2
    elif 'Omeganu' in cambpars:
        cambpars['omnuh2'] = cambpars.pop('Omeganu')*h2
        onuh2 = cambpars['omnuh2']
    elif 'omnuh2' in cambpars:
        onuh2 = cambpars['omnuh2']
    if 'Omegam' in cambpars: #TO BE GENERALIZED
        cambpars['omch2']  = cambpars.pop('Omegam')*h2-cambpars['ombh2']-onuh2

    rescaleAs = False
    if 'sigma8' in cambpars:
        insigma8 = cambpars['sigma8']
        cambpars['As'] = 2.1e-9
        cambpars.pop('sigma8')
        rescaleAs = True
    
    try:
        myparpars=camb.set_params(**cambpars) # to see which methods are being called: verbose=True
    except camb.CAMBUnknownArgumentError as argument:
        print("Remove parameter from cambparams: ", str(argument))

    if rescaleAs==True:
        cambpars['As'] = rescale_LP(cambpars,camb,insigma8)

    return cambpars

def rescale_LP(cambpars,camb,insigma8) :
    cambpars_LP = cambpars.copy()
    ini_As = 2.1e-9
    
    boost  = 3
    cambpars_LP['AccuracyBoost'] = boost
    cambpars_LP['lAccuracyBoost'] = boost
    cambpars_LP['lSampleBoost'] = boost
    cambpars_LP['kmax'] = 20
    
    pars = camb.set_params(redshifts=[0.],**cambpars_LP)
    results = camb.get_results(pars)
    
    test_sig8=np.array(results.get_sigma8())
    final_As = ini_As*(insigma8/test_sig8[-1])**2.

    return final_As


def ready_camb(cosmopars,camb_dict):
    #Input: cosmological parameters in a dict
    #Output: Interpolating funktion for Pmm and Pcb
    zmax = 5  
    zsamples= 100
    camb_zarray = (np.linspace(0. , zmax, zsamples)[::-1])

    cambcosmopars = basechange_camb(cosmopars)

    cambpars  = deepcopy(camb_dict)
    cambpars.update(cambcosmopars)

    cambclasspars=camb.set_params(**cambpars)

    cambclasspars.set_matter_power(redshifts=camb_zarray, 
                                        k_per_logint=cambpars['k_per_logint'], 
                                        kmax=cambpars['kmax'],
                                        accurate_massive_neutrino_transfers=cambpars['accurate_massive_neutrino_transfers']
                                        )
    cambres = camb.get_results(cambclasspars)
    return cambres

def ready_class(cosmopars,class_dict):
    #Input: cosmological parameters in a dict
    #Output: Interpolating funktion for Pmm and Pcb

    classcosmopars = basechange_class(cosmopars)
    classpars  = deepcopy(class_dict)
    classpars.update(classcosmopars)
    
    classres = Class()
    classres.set(classpars)
    print(classpars)
    classres.compute()
    
    return classres    

class classresults: #will get you your matterpowerspectra interpolators from class
    def __init__(self,cosmopars,class_dict_linear,class_dict_nonlinear):
        classres_linear = ready_class(cosmopars=cosmopars,class_dict=class_dict_linear)
        
        #Calculate the Matter fractions for CB Powerspectrum
        f_cdm= classres_linear.Omega0_cdm()/classres_linear.Omega_m()
        f_b  = classres_linear.Omega_b()   /classres_linear.Omega_m()
        f_cb = f_cdm+f_b
        f_nu = 1-f_cb

        ## rows are k, and columns are z
        ## interpolating function Pk_l (k,z)
        Pk_l, k, z    = classres_linear.get_pk_and_k_and_z(nonlinear=False)
        Pk_cb_l, k, z = classres_linear.get_pk_and_k_and_z(only_clustering_species=True,nonlinear=False)        
        self.Pk_l = RectBivariateSpline(z[::-1],k,(np.flip(Pk_l,axis=1)).transpose())
        self.Pk_cb_l = RectBivariateSpline(z[::-1],k,(np.flip(Pk_cb_l,axis=1)).transpose())
        self.kgrid = k
        self.zgrid = z[::-1]


        classres_nonlinear = ready_class(cosmopars=cosmopars,class_dict=class_dict_nonlinear)
        ## interpolating function Pk_nl (k,z)
        Pk_nl, k, z = classres_nonlinear.get_pk_and_k_and_z(nonlinear=True)
        self.Pk_nl = RectBivariateSpline(z[::-1],k,(np.flip(Pk_nl,axis=1)).transpose())

        tk, k ,z = classres_linear.get_transfer_and_k_and_z()
        T_cb = (f_b*tk['d_b']+f_cdm*tk['d_cdm'])/f_cb
        T_nu = tk['d_ncdm[0]']

        pm = classres_linear.get_primordial()
        pk_prim = UnivariateSpline(pm['k [1/Mpc]'],pm['P_scalar(k)'])(k)*(2.*np.pi**2)/np.power(k,3)
        
        pk_cnu  = T_nu * T_cb * pk_prim[:,None]
        pk_nunu = T_nu * T_nu * pk_prim[:,None]
        Pk_cb_nl= 1./f_cb**2 * (Pk_nl-2*pk_cnu*f_nu*f_cb-pk_nunu*f_nu*f_nu)

        self.Pk_cb_nl = RectBivariateSpline(z[::-1],k,(np.flip(Pk_cb_nl,axis=1)).transpose())
        classres_linear.struct_cleanup()
        classres_linear.empty()
        classres_nonlinear.struct_cleanup()
        classres_nonlinear.empty()

class cambresults: #will get you your matterpowerspectra interpolators from camb
    def __init__(self,cosmopars,camb_dict):
        cambres = ready_camb(cosmopars=cosmopars,camb_dict=camb_dict)
        Pk_l, self.zgrid, self.kgrid = cambres.get_matter_power_interpolator(
                                                     hubble_units=False,
                                                     k_hunit=False,
                                                     var1='delta_tot',
                                                     var2='delta_tot',
                                                     nonlinear=False,
                                                     extrap_kmax=100,
                                                     return_z_k=True)
        Pk_nl = cambres.get_matter_power_interpolator(
                                                     hubble_units=False,
                                                     k_hunit=False,
                                                     var1='delta_tot',
                                                     var2='delta_tot',
                                                     nonlinear=True,
                                                     extrap_kmax=100,
                                                     return_z_k=False)
        Pk_cb_l = cambres.get_matter_power_interpolator(
                                                     hubble_units=False,
                                                     k_hunit=False,
                                                     var1='delta_nonu',
                                                     var2='delta_nonu',
                                                     nonlinear=False,
                                                     extrap_kmax=100,
                                                     return_z_k=False)

        self.Pk_l = RectBivariateSpline(self.zgrid, self.kgrid, Pk_l.P(self.zgrid, self.kgrid))
        self.Pk_nl = RectBivariateSpline(self.zgrid, self.kgrid, Pk_nl.P(self.zgrid, self.kgrid))
        self.Pk_cb_l = RectBivariateSpline(self.zgrid, self.kgrid, Pk_cb_l.P(self.zgrid, self.kgrid))
        
        Om_m = (cambres.get_Omega('cdm', z=0)+cambres.get_Omega('baryon', z=0)+cambres.get_Omega('nu', z=0))
        
        #Calculate the Non linear cb powerspectrum using Gabrieles Approximation
        f_cdm=cambres.get_Omega('cdm',z=0)/Om_m
        f_b  =cambres.get_Omega('baryon',z=0)/Om_m
        f_cb =f_cdm+f_b
        f_nu =1-f_cb
        Pk_cross_l = cambres.get_matter_power_interpolator(
                                                     hubble_units=False,
                                                     k_hunit=False,
                                                     var1='delta_nonu',
                                                     var2='delta_nu',
                                                     nonlinear=False,
                                                     extrap_kmax=100,
                                                     return_z_k=False)
        Pk_nunu_l = cambres.get_matter_power_interpolator(
                                                     hubble_units=False,
                                                     k_hunit=False,
                                                     var1='delta_nu',
                                                     var2='delta_nu',
                                                     nonlinear=False,
                                                     extrap_kmax=100,
                                                     return_z_k=False)
        Pk_cb_nl=1/f_cb**2 * (Pk_nl.P(self.zgrid, self.kgrid)-2 * Pk_cross_l.P(self.zgrid, self.kgrid)*f_cb*f_nu - Pk_nunu_l.P(self.zgrid, self.kgrid) * f_nu**2)
        self.Pk_cb_nl = RectBivariateSpline(self.zgrid, self.kgrid, Pk_cb_nl)

#####################
#Start of the Script#
#####################
k = np.logspace(-3,1,num=200)
z = [0]

# Compare Fiducial MM Powerspectra
####################################################################
classres_DP=classresults(cosmopars=fiducial,class_dict_linear=Class_dict_DP,class_dict_nonlinear=Class_dict_DP)
classres_HP=classresults(cosmopars=fiducial,class_dict_linear=Class_dict_HP,class_dict_nonlinear=Class_dict_HP)
classres_UHP=classresults(cosmopars=fiducial,class_dict_linear=Class_dict_UHP,class_dict_nonlinear=Class_dict_HP_lowk)

Class_mm_lin_DP = np.array([classres_DP.Pk_l (zi,k) for zi in z])[:,0,:]
Class_mm_nol_DP = np.array([classres_DP.Pk_nl(zi,k) for zi in z])[:,0,:]
Class_mm_lin_HP = np.array([classres_HP.Pk_l (zi,k) for zi in z])[:,0,:]
Class_mm_nol_HP = np.array([classres_HP.Pk_nl(zi,k) for zi in z])[:,0,:]
Class_mm_lin_UHP =np.array([classres_UHP.Pk_l (zi,k) for zi in z])[:,0,:]
Class_mm_nol_UHP =np.array([classres_UHP.Pk_nl(zi,k) for zi in z])[:,0,:]

###################################################################
cambres_DP = cambresults(cosmopars=fiducial,camb_dict=Camb_dict_DP)
cambres_HP = cambresults(cosmopars=fiducial,camb_dict=Camb_dict_HP)

Camb_mm_lin_DP  = np.array([cambres_DP.Pk_l (zi,k) for zi in z]) [:,0,:]
Camb_mm_nol_DP  = np.array([cambres_DP.Pk_nl(zi,k)for zi in z])  [:,0,:]
Camb_mm_lin_HP  = np.array([cambres_HP.Pk_l (zi,k) for zi in z]) [:,0,:]
Camb_mm_nol_HP  = np.array([cambres_HP.Pk_nl(zi,k)for zi in z])  [:,0,:]

#####################################################################
fig , axs = plt.subplots(2,1,figsize=(16,9),sharex=True,sharey=True)
axs = np.array(axs).flatten()

snscolors=sns.color_palette("colorblind")
colors = iter(snscolors)


c = next(colors)
axs[0].axvspan(1e-3,0.25,hatch='/',alpha=0.2,color=c,label='spectroscopic')
c = next(colors)
axs[1].axvspan(1e-3,10,hatch='\\',alpha=0.2,color=c,label='photometric')
axs[0].axvspan(1e-3,1e-3,hatch='\\',alpha=0.2,color=c,label='photometric')



buf=[]

c = next(colors)
axs[0].plot(k,(Class_mm_lin_DP[0,:]-Camb_mm_lin_HP[0,:])/(Camb_mm_lin_HP[0,:])*100,ls='-',color=c)
axs[1].plot(k,(Class_mm_nol_DP[0,:]-Camb_mm_nol_HP[0,:])/(Camb_mm_nol_HP[0,:])*100,ls='-' ,color=c)
axs[0].scatter(buf,buf,c=c,label='CLASS:DP')

c = next(colors)
axs[0].plot(k,(Class_mm_lin_HP[0,:]-Camb_mm_lin_HP[0,:])/(Camb_mm_lin_HP[0,:])*100,ls='-',color=c)
axs[1].plot(k,(Class_mm_nol_HP[0,:]-Camb_mm_nol_HP[0,:])/(Camb_mm_nol_HP[0,:])*100,ls='-' ,color=c)
axs[0].scatter(buf,buf,c=c,label='CLASS:HP')

c = next(colors)
axs[0].plot(k,(Class_mm_lin_UHP[0,:]-Camb_mm_lin_HP[0,:])/(Camb_mm_lin_HP[0,:])*100,ls='-',color=c)
axs[1].plot(k,(Class_mm_nol_UHP[0,:]-Camb_mm_nol_HP[0,:])/(Camb_mm_nol_HP[0,:])*100,ls='-' ,color=c)
axs[0].scatter(buf,buf,c=c,label='CLASS:UHP')


axs[0].set_xscale('log')
axs[1].set_xscale('log')

axs[0].set_title('Linear Power Spectrum')
axs[1].set_title('Non-Linear Power Spectrum')

axs[0].set_ylabel('% rel. dev. vs CAMB:HP $\mathrm{P}_\mathrm{mm}$')
axs[1].set_ylabel('% rel. dev. vs CAMB:HP $\mathrm{P}_\mathrm{mm}$')

axs[1].set_xlabel('k [1/Mpc]')
axs[1].set_xlim(1e-3,10)


axs[0].legend()

plt.savefig('Fiducial_deviation_mm_appendix_A.pdf')

# # Compare Log-Derivative of CB Powerspectra wrt deriv_param
# deriv_param= 'wa'
# deriv_eps  = 0.01

# meps = deepcopy(fiducial)
# peps = deepcopy(fiducial)
# meps[deriv_param] = meps[deriv_param]*(1-deriv_eps)
# peps[deriv_param] = peps[deriv_param]*(1+deriv_eps)

# mepsclassres=classresults(cosmopars=meps)
# mepscambres = cambresults(cosmopars=meps)
# pepsclassres=classresults(cosmopars=peps)
# pepscambres = cambresults(cosmopars=peps)


# mepsClass_cb_lin = np.array([mepsclassres.Pk_cb_l (zi,k) for zi in z])[:,0,:]
# mepsClass_cb_nol = np.array([mepsclassres.Pk_cb_nl(zi,k) for zi in z])[:,0,:]
# mepsCamb_cb_lin =  np.array([mepscambres.Pk_cb_l (zi,k) for zi in z])[:,0,:]
# mepsCamb_cb_nol =  np.array([mepscambres.Pk_cb_nl(zi,k)for zi in z])[:,0,:]
# pepsClass_cb_lin = np.array([pepsclassres.Pk_cb_l (zi,k) for zi in z])[:,0,:]
# pepsClass_cb_nol = np.array([pepsclassres.Pk_cb_nl(zi,k) for zi in z])[:,0,:]
# pepsCamb_cb_lin =  np.array([pepscambres.Pk_cb_l (zi,k) for zi in z])[:,0,:]
# pepsCamb_cb_nol =  np.array([pepscambres.Pk_cb_nl(zi,k)for zi in z])[:,0,:]

# dClass_cb_lin=(pepsClass_cb_lin-mepsClass_cb_lin)/Class_cb_lin
# dClass_cb_nol=(pepsClass_cb_nol-mepsClass_cb_nol)/Class_cb_nol
# dCamb_cb_lin =(pepsCamb_cb_lin-mepsCamb_cb_lin)/Camb_cb_lin
# dCamb_cb_nol =(pepsCamb_cb_nol-mepsCamb_cb_nol)/Camb_cb_nol

# fig , axs = plt.subplots(2,2,figsize=(16,9))
# axs = np.array(axs).flatten()

# for 0 , axi in enumerate(axs):
#     axi.plot(k,(dClass_cb_lin[0,:]-dCamb_cb_lin[0,:])/(dClass_cb_lin[0,:]+dCamb_cb_lin[0,:])*200,ls='--',c='black')
#     axi.plot(k,(dClass_cb_nol[0,:]-dCamb_cb_nol[0,:])/(dClass_cb_nol[0,:]+dCamb_cb_nol[0,:])*200,ls='-',c='black')
#     axi.set_xscale('log')
#     axi.set_title('z: {}'.format(z[0]))
#     axi.set_ylabel('% rel. dev. Class vs Camb $\partial \, \mathrm{P}_\mathrm{cb}$')
#     axi.axvspan(1e-3,0.3,hatch='/',alpha=0.05,color='green',label='optimistic')
#     axi.axvspan(1e-3,0.25,hatch='\\',alpha=0.05,color='red',label='pessimistic')


# axs[2].set_xlabel('k [1/Mpc]')
# axs[3].set_xlabel('k [1/Mpc]')
# axs[0].plot(buf,buf,ls='--',c='black',label='linear')
# axs[0].plot(buf,buf,ls='-',c='black',label='non linear')
# axs[0].legend()

# plt.savefig('Derivativ_deviation_cb.png')

# mepsClass_mm_lin = np.array([mepsclassres.Pk_l (zi,k) for zi in z])[:,0,:]
# mepsClass_mm_nol = np.array([mepsclassres.Pk_nl(zi,k) for zi in z])[:,0,:]
# mepsCamb_mm_lin =  np.array([mepscambres.Pk_l (zi,k) for zi in z])[:,0,:]
# mepsCamb_mm_nol =  np.array([mepscambres.Pk_nl(zi,k)for zi in z])[:,0,:]
# pepsClass_mm_lin = np.array([pepsclassres.Pk_l (zi,k) for zi in z])[:,0,:]
# pepsClass_mm_nol = np.array([pepsclassres.Pk_nl(zi,k) for zi in z])[:,0,:]
# pepsCamb_mm_lin =  np.array([pepscambres.Pk_l (zi,k) for zi in z])[:,0,:]
# pepsCamb_mm_nol =  np.array([pepscambres.Pk_nl(zi,k)for zi in z])[:,0,:]

# dClass_mm_lin=(pepsClass_mm_lin-mepsClass_mm_lin)/Class_mm_lin
# dClass_mm_nol=(pepsClass_mm_nol-mepsClass_mm_nol)/Class_mm_nol
# dCamb_mm_lin =(pepsCamb_mm_lin-mepsCamb_mm_lin)/Camb_mm_lin
# dCamb_mm_nol =(pepsCamb_mm_nol-mepsCamb_mm_nol)/Camb_mm_nol

# fig , axs = plt.subplots(2,2,figsize=(16,9))
# axs = np.array(axs).flatten()

# for 0 , axi in enumerate(axs):
#     axi.plot(k,(dClass_mm_lin[0,:]-dCamb_mm_lin[0,:])/(dClass_mm_lin[0,:]+dCamb_mm_lin[0,:])*200,ls='--',c='black')
#     axi.plot(k,(dClass_mm_nol[0,:]-dCamb_mm_nol[0,:])/(dClass_mm_nol[0,:]+dCamb_mm_nol[0,:])*200,ls='-',c='black')
#     axi.set_xscale('log')
#     axi.set_title('z: {}'.format(z[0]))
#     axi.set_ylabel('% rel. dev. Class vs Camb $\partial \, \mathrm{P}_\mathrm{mm}$')
#     axi.axvspan(1e-3,0.3,hatch='/',alpha=0.05,color='green',label='optimistic')
#     axi.axvspan(1e-3,0.25,hatch='\\',alpha=0.05,color='red',label='pessimistic')


# axs[2].set_xlabel('k [1/Mpc]')
# axs[3].set_xlabel('k [1/Mpc]')
# axs[0].plot(buf,buf,ls='--',c='black',label='linear')
# axs[0].plot(buf,buf,ls='-',c='black',label='non linear')
# axs[0].legend()

# plt.savefig('Derivativ_deviation_mm.png')
