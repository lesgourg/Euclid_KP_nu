# before running, since this is the pessimistic case, go to the root montepython/ directory and edit:

    montepython/likelihoods/euclid_photometric/euclid_photometric.data

# check that we use only the weak lensing part, that is,

    euclid_photometric.probe = ['WL']

# and

    euclid_photometric.fiducial_file  = "euclid_wl_fiducial.dat"

#  check that the pessimistic settings are on, that is:

    euclid_photometric.lmax_WL = 1500

# We assume here that you run within the montepython/ directory, not within this one. Then the commands are:

1) to remove a possible previous fiducial model generated with different settings

    rm data/euclid_wl_fiducial.dat

2) to remove possible results of a previous attempt to run in the same directory

    rm -rf ../Euclid_KP_nu/results/montepython_fisher/wl_only/pessimistic

3) to create a new fiducial model

    python3 montepython/MontePython.py run -p ../Euclid_KP_nu/input/montepython_fisher/wl_only/pessimistic/wl_only_pess.param -o ../Euclid_KP_nu/results/montepython_fisher/wl_only/pessimistic -f 0

4) to run around this fiducial model

    python3 montepython/MontePython.py run -o ../Euclid_KP_nu/results/montepython_fisher/wl_only/pessimistic --fisher --fisher-step-it 1 --fisher-tol 10000

5) to create the .paramnames file understood by cosmicfish

    cd ../Euclid_KP_nu
    python3 input/montepython_fisher/paramnames_for_cosmicfish.py results/montepython_fisher/wl_only/pessimistic