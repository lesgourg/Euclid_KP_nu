# before running, since this is the optimistic case, go to the root montepython/ directory and edit:

    montepython/likelihoods/euclid_photometric/euclid_photometric.data

# and check that the pessimistic settings are on, that is:

    euclid_photometric.fiducial_file  = "euclid_xc_fiducial_hmcode.dat"

    euclid_photometric.probe = ['WL_GCph_XC']

    euclid_photometric.lmax_WL = 5000
    euclid_photometric.lmax_GC = 3000
    euclid_photometric.lmax_XC = 3000

    also, you should have:

    euclid_photometric.use_halofit = False

# We assume here that you run within the montepython/ directory, not within this one. Then the commands are:

1) to remove a possible previous fiducial model generated with different settings
    rm data/euclid_xc_fiducial_hmcode.dat

2) to remove possible results of a previous attempt to run in the same directory

    rm -rf ../Euclid_KP_nu/results/montepython_fisher/photometric/optimistic

3) to create a new fiducial model

    python3 montepython/MontePython.py run -p ../Euclid_KP_nu/input/montepython_fisher/photometric/optimistic/photometric_opt.param -o ../Euclid_KP_nu/results/montepython_fisher/photometric/optimistic -f 0

4) to run around this fiducial model

    python3 montepython/MontePython.py run -o ../Euclid_KP_nu/results/montepython_fisher/photometric/optimistic --fisher --fisher-step-it 1 --fisher-tol 10000

5) to create the .paramnames file understood by cosmicfish

    cd ../Euclid_KP_nu
    python3 input/montepython_fisher/paramnames_for_cosmicfish.py results/montepython_fisher/photometric/optimistic