# before running, since this is the pessimistic case, go to the root montepython/ directory and edit:

    montepython/likelihoods/euclid_spectroscopic/euclid_spectroscopic.data

# and check that the pessimistic settings are on, that is:

    euclid_spectroscopic.kmax = 0.15

# We assume here that you run within the montepython/ directory, not within this one. Then the commands are:

1) to remove a possible previous fiducial model generated with different settings

    rm data/euclid_pk_fiducial.dat

2) to remove possible results of a previous attempt to run in the same directory

    rm -rf ../Euclid_KP_nu/results/montepython_fisher/spectroscopic/quasilinear

3) to create a new fiducial model

    python3 montepython/MontePython.py run -p ../Euclid_KP_nu/input/montepython_fisher/spectroscopic/quasilinear/spectroscopic_quasi.param -o ../Euclid_KP_nu/results/montepython_fisher/spectroscopic/quasilinear -f 0

4) to run around this fiducial model

    python3 montepython/MontePython.py run -o ../Euclid_KP_nu/results/montepython_fisher/spectroscopic/quasilinear --fisher --fisher-step-it 1 --fisher-tol 10000

5) to create the .paramnames file understood by cosmicfish

    cd ../Euclid_KP_nu
    python3 input/montepython_fisher/paramnames_for_cosmicfish.py results/montepython_fisher/spectroscopic/quasilinear