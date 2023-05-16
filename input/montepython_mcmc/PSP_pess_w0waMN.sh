#!/bin/zsh
CASE=PSP_w0waMN
rm -r results/montepython_mcmc/${CASE}
cd ../montepython
cp montepython/likelihoods/euclid_spectroscopic/euclid_spectroscopic.data.pessimistic montepython/likelihoods/euclid_spectroscopic/euclid_spectroscopic.data
cp montepython/likelihoods/euclid_photometric_z/euclid_photometric_z.data.pessimistic montepython/likelihoods/euclid_photometric_z/euclid_photometric_z.data
rm data/euclid_xc_fiducial.dat
rm data/euclid_pk_fiducial.dat
rm data/fake_planck_realistic_fiducial.dat
python3 montepython/MontePython.py run -p ../Euclid_KP_nu/input/montepython_mcmc/${CASE}.param -o ../Euclid_KP_nu/results/montepython_mcmc/${CASE} -f 0
python3 montepython/MontePython.py run -o ../Euclid_KP_nu/results/montepython_mcmc/${CASE} -f 0 -N 1
python3 montepython/MontePython.py run -o ../Euclid_KP_nu/results/montepython_mcmc/${CASE} -N 100000 --update 50 --superupdate 20 --conf default.conf # add your covmat here: -c ...covmat
cd ../Euclid_KP_nu
