Go to the input_4_cast directory to run the code, and check that you have correctly checked out the branch wp3, not the main branch

    cd ../input_4_cast
    git status

----------------------------

* For producing the CLASS external files, you will use the input files

        ../Euclid_KP_nu/input/input_4_cast/default_class_WP3_LP.ini
        ../Euclid_KP_nu/input/input_4_cast/default_class_WP3_MP.ini
        ../Euclid_KP_nu/input/input_4_cast/default_class_WP3_HP.ini

    (for respectively low precision, medium precision, high precision).
    The relative paths at the beginning may need some editing if, on your machine, class is not located in the directory ../class/

    Then run with e.g.

        python3 run.py ../Euclid_KP_nu/input/input_4_cast/default_class_WP3_LP.ini

----------------------------

* For producing the CAMB external files, you will use the input file

        ../Euclid_KP_nu/input/input_4_cast/default_camb_WP3_LP.ini
        ../Euclid_KP_nu/input/input_4_cast/default_camb_WP3_MP.ini
        ../Euclid_KP_nu/input/input_4_cast/default_camb_WP3_HP.ini

    (for respectively low precision, medium precision, high precision).
    The relative paths at the beginning may need some editing if, on your machine, class is not located in the directory ../CAMB/

    Also, set number_of_threads to the number of cores in your machine

    Then run with e.g.

        python3 run.py ../Euclid_KP_nu/input/input_4_cast/default_camb_WP3_LP.ini

-----------------------------

* For comparing CLASS and CAMB files, run e.g.

    python3 compare.py --save_plots --threshold 0.1 output/default_camb_euclid_WP3_LP output/default_class_euclid_WP3_LP


    All spectra should agree very well (typically better than 0.1% for
    the linear spectrum with MP or HP) excepted when considering the
    case of CAMB with N_eff[fiducial]>3 and N_eff[minus epsilon]<3,
    for which the CAMB treatement leads to a small discontinuity. So,
    for large epsilons, it is normal that the spectra of
    Neff_mn_epsilon disagree by a greater amount, e.g. 0.3%.  Thus for
    N_eff[fiducial]=3.044 it is not recommended to rely on the CAMB
    Neff_mn_epsilon output for epsilon euqal to 0.015 or greater.

-----------------------------

When everything is done, go backs

     cd ../Euclid_KP_nu
