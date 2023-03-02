## CAMB

```
class: <CAMBparams>
 WantCls = False
 WantTransfer = True
 WantScalars = True
 WantTensors = False
 WantVectors = False
 WantDerivedParameters = True
 Want_cl_2D_array = False
 Want_CMB = False
 Want_CMB_lensing = False
 DoLensing = False
 NonLinear = NonLinear_both
 Transfer: <TransferParams>
   high_precision = True
   accurate_massive_neutrinos = True
   kmax = 50.0
   k_per_logint = 50
   PK_num_redshifts = 100
   PK_redshifts = [5.0, 4.94949494949495, 4.898989898989899, 4.848484848484849, 4.797979797979798, 4.747474747474747, 4.696969696969697, ...]
 want_zstar = False
 want_zdrag = False
 min_l = 2
 max_l = 4000
 max_l_tensor = 600
 max_eta_k = 18000.0
 max_eta_k_tensor = 1200.0
 ombh2 = 0.02232998788914099
 omch2 = 0.11980025128030217
 omk = 0.0
 omnuh2 = 0.0006451383989381787
 H0 = 67.36999999999999
 TCMB = 2.7255
 YHe = 0.24
 num_nu_massless = 2.044
 num_nu_massive = 1
 nu_mass_eigenstates = 1
 share_delta_neff = True
 nu_mass_degeneracies = [1.0153333333333332]
 nu_mass_fractions = [1.0]
 nu_mass_numbers = [1]
 InitPower: <InitialPowerLaw>
   tensor_parameterization = tensor_param_rpivot
   ns = 0.96605
   nrun = 0.0
   nrunrun = 0.0
   nt = -0.0
   ntrun = -0.0
   r = 0.0
   pivot_scalar = 0.05
   pivot_tensor = 0.05
   As = 2.0944696981908336e-09
   At = 1.0
 Recomb: <Recfast>
   min_a_evolve_Tm = 0.0011098779505118728
   RECFAST_fudge = 1.125
   RECFAST_fudge_He = 0.86
   RECFAST_Heswitch = 6
   RECFAST_Hswitch = True
   AGauss1 = -0.14
   AGauss2 = 0.079
   zGauss1 = 7.28
   zGauss2 = 6.73
   wGauss1 = 0.18
   wGauss2 = 0.33
 Reion: <TanhReionization>
   Reionization = False
   use_optical_depth = True
   redshift = 10.0
   optical_depth = 0.058
   delta_redshift = 0.5
   fraction = -1.0
   include_helium_fullreion = True
   helium_redshift = 3.5
   helium_delta_redshift = 0.4
   helium_redshiftstart = 5.5
   tau_solve_accuracy_boost = 1.0
   timestep_boost = 1.0
   max_redshift = 50.0
 DarkEnergy: <DarkEnergyPPF>
   w = -1.0
   wa = 0.0
   cs2 = 1.0
   use_tabulated_w = False
 NonLinearModel: <Halofit>
   Min_kh_nonlinear = 0.005
   halofit_version = mead2020
   HMCode_A_baryon = 3.13
   HMCode_eta_baryon = 0.603
   HMCode_logT_AGN = 7.8
 Accuracy: <AccuracyParams>
   AccuracyBoost = 3.0
   lSampleBoost = 1.0
   lAccuracyBoost = 3.0
   AccuratePolarization = True
   AccurateBB = False
   AccurateReionization = True
   TimeStepBoost = 1.0
   BackgroundTimeStepBoost = 1.0
   IntTolBoost = 1.0
   SourcekAccuracyBoost = 1.0
   IntkAccuracyBoost = 1.0
   TransferkBoost = 1.0
   NonFlatIntAccuracyBoost = 1.0
   BessIntBoost = 1.0
   LensingBoost = 1.0
   NonlinSourceBoost = 1.0
   BesselBoost = 1.0
   LimberBoost = 1.0
   SourceLimberBoost = 1.0
   KmaxBoost = 1.0
   neutrino_q_boost = 1.0
 SourceTerms: <SourceTermParams>
   limber_windows = True
   limber_phi_lmin = 100
   counts_density = True
   counts_redshift = True
   counts_lensing = False
   counts_velocity = True
   counts_radial = False
   counts_timedelay = True
   counts_ISW = True
   counts_potential = True
   counts_evolve = False
   line_phot_dipole = False
   line_phot_quadrupole = False
   line_basic = True
   line_distortions = True
   line_extra = False
   line_reionization = False
   use_21cm_mK = True
 z_outputs = []
 scalar_initial_condition = initial_adiabatic
 InitialConditionVector = []
 OutputNormalization = 1
 Alens = 1.0
 MassiveNuMethod = Nu_best
 DoLateRadTruncation = True
 Evolve_baryon_cs = False
 Evolve_delta_xe = False
 Evolve_delta_Ts = False
 Do21cm = False
 transfer_21cm_cl = False
 Log_lvalues = False
 use_cl_spline_template = True
 SourceWindows = []
 CustomSources: <CustomSources>
   num_custom_sources = 0
   c_source_func = None
   custom_source_ell_scales = []
```


## CLASS 

```
{'k_per_decade_for_bao': 50, 
'k_per_decade_for_pk': 50, 
'l_max_g': 20, 'l_max_pol_g': 15, 
'radiation_streaming_approximation': 2, 
'radiation_streaming_trigger_tau_over_tau_k': 240.0, 
'radiation_streaming_trigger_tau_c_over_tau': 100.0, 
'tol_ncdm_synchronous': 1e-05, 'l_max_ncdm': 22,
'ncdm_fluid_trigger_tau_over_tau_k': 100.0, 
'background_Nloga': 6000, 'thermo_Nz_log': 20000, 
'thermo_Nz_lin': 40000, 'tol_perturbations_integration': 1e-06,
'P_k_max_1/Mpc': 50, 'output': 'mPk,mTk', 'non linear': 'hmcode',
'z_max_pk': 5.0, 'N_ncdm': 1, 'T_cmb': 2.7255, 'nonlinear_min_k_max': 80.0, 
'hmcode_tol_sigma': 1e-08, 'k_pivot': 0.05, 
'Omega_k': 0, 'YHe': 0.2454006, 'sigma8': 0.81, 'h': 0.6737,
'N_ur': 2.029333333333333, 'T_ncdm': 0.7163687246184776, 'Omega_ncdm': 0.0014207234756587472, 
'Omega_b': 0.0491989, 'Omega_cdm': 0.26395180382434125, 'n_s': 0.96605}
```
