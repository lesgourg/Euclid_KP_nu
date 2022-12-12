PYTHON=python3

# The tests below are performed using the spectra P_L, P_NL produced by input_4_cast for class/camb with different precision settings.

# First test: comparison of fiducial P_L and P_NL
# This test is not done with the output files located in
#
# Euclid_w0wa/input/input_4_cast
#
# Indeed, for this test, we need precisely the same value of A_s, and
# thus slightly different values of sigma_8.
# So we use here some customised input files with two different sigma_8
# giving the same A_s. However, we can do short runs because we only need
# the fiducial models, no parameter variations.
# Additionally, for this test, we use 'halofit_min_k_nonlinear=3.282e-3' in CLASS settings,
# to match the L/NL transition int he two codes. We do not do it elsewhere.

echo "Run and plot comparison of fiducial spectra? (y/n)"
read answer
if [ "$answer" = "y" ] ; then

    echo "For this purpose, do you want to do short runs of input_4_cast with customised input files?  (y/n)"
    read answer
    if [ "$answer" = "y" ] ; then

        cd ../../input_4_cast
        $PYTHON run.py ../Euclid_w0wa/test_precision/class_w0wa_DP.ini
        $PYTHON run.py ../Euclid_w0wa/test_precision/class_w0wa_HP.ini
        $PYTHON run.py ../Euclid_w0wa/test_precision/camb_w0wa_P1.ini
        $PYTHON run.py ../Euclid_w0wa/test_precision/camb_w0wa_P2.ini
        $PYTHON run.py ../Euclid_w0wa/test_precision/camb_w0wa_P3.ini
        cd ../Euclid_w0wa/test_precision/

    fi

    echo "Now doing the plot..."
    $PYTHON plot_fiducial_ratios.py

fi

# Second test: comparison of first derivatives of P_L and P_NL w.r.t each parameter
# This test is done with the output files located in
#
# Euclid_w0wa/input/input_4_cast

echo "Run and plot comparison of first derivatives of spectra? (y/n)"
read answer
if [ "$answer" = "y" ] ; then

    echo "For this purpose, do you want to do long runs of input_4_cast with normal input files? Usually this would have alredy been done and you don't want to erase it...? (y/n)"
    read answer
    if [ "$answer" = "y" ] ; then

        cd ../../input_4_cast
        $PYTHON run.py ../Euclid_w0wa/input/input_4_cast/class_w0wa_DP.ini
        $PYTHON run.py ../Euclid_w0wa/input/input_4_cast/class_w0wa_HP.ini
        $PYTHON run.py ../Euclid_w0wa/input/input_4_cast/camb_w0wa_P1.ini
        $PYTHON run.py ../Euclid_w0wa/input/input_4_cast/camb_w0wa_P2.ini
        $PYTHON run.py ../Euclid_w0wa/input/input_4_cast/camb_w0wa_P3.ini
        cd ../Euclid_w0wa/test_precision/

    fi

    echo "Now doing the plot..."
    $PYTHON plot_derivative.py
fi

# Third test:
# For photo/opt and then for spectro/opt:
# - run CosmicFish with CAMBext P2 versus CAMB ext P3, and compare (un)marginalised errors;
# - run CosmicFish with and CLASSext DP versus CLASSext HP, and compare (un)marginalised errors.

PROBE=photometric
PROBE_SHORT=photo
LKL=euclid_photometric
CASE=optimistic
CASE_SHORT=opt

#########

echo "For photometric/optimistic:"
echo "Shall we run CosmicFish CLASS external DP? (y/n)"
read answer
if [ "$answer" = "y" ] ; then
    $PYTHON ../input/cosmicfish/$PROBE/$CASE/class_external_DP.py
fi

echo "For photometric/optimistic:"
echo "Shall we run CosmicFish CLASS external HP? (y/n)"
read answer
if [ "$answer" = "y" ] ; then
    $PYTHON ../input/cosmicfish/$PROBE/$CASE/class_external_HP.py
fi

echo "For photometric/optimistic:"
echo "For CosmicFish/external/CLASS DP vs HP:"
echo "Shall we recompute error ratios and redo error plots ? (y/n)"
read answer
if [ "$answer" = "y" ] ; then

    # plot error comparisons
    cd ../plots
    #
    $PYTHON $PROBE/$CASE/CF_class_ext_DP-vs-HP.py --error-only
    cd ../test_precision

fi

############

echo "For photometric/optimistic:"
echo "Shall we run CosmicFish CAMB external P2? (y/n)"
read answer
if [ "$answer" = "y" ] ; then
    $PYTHON ../input/cosmicfish/$PROBE/$CASE/camb_external_P2.py
fi

echo "For photometric/optimistic:"
echo "Shall we run CosmicFish CAMB external P3? (y/n)"
read answer
if [ "$answer" = "y" ] ; then
    $PYTHON ../input/cosmicfish/$PROBE/$CASE/camb_external_P3.py
fi

echo "For photometric/optimistic:"
echo "For CosmicFish/external/CAMB P2 vs P3:"
echo "Shall we recompute error ratios and redo error plots ? (y/n)"
read answer
if [ "$answer" = "y" ] ; then

    # plot error comparisons
    cd ../plots
    #
    $PYTHON $PROBE/$CASE/CF_camb_ext_P2-vs-P3.py --error-only
    cd ../test_precision

fi

####################

echo "For photometric/optimistic:"
echo "Shall we run MontePython Fisher CLASS DP ? (y/n)"
read answer
if [ "$answer" = "y" ] ; then

    cd ../../montepython
    rm -r ../Euclid_w0wa/results/montepython_fisher/$PROBE/${CASE}_DP
    cp montepython/likelihoods/$LKL/$LKL.data.$CASE montepython/likelihoods/$LKL/$LKL.data
    rm data/euclid_xc_fiducial.dat
    $PYTHON montepython/MontePython.py run -p ../Euclid_w0wa/input/montepython_fisher/$PROBE/$CASE/${PROBE}_${CASE_SHORT}_DP.param -o ../Euclid_w0wa/results/montepython_fisher/$PROBE/${CASE}_DP -f 0
    $PYTHON montepython/MontePython.py run -o ../Euclid_w0wa/results/montepython_fisher/$PROBE/${CASE}_DP --fisher --fisher-step-it 1 --fisher-tol 10000
    cd ../Euclid_w0wa/test_precision

fi

echo "For photometric/optimistic:"
echo "Shall we run MontePython Fisher CLASS HP ? (y/n)"
read answer
if [ "$answer" = "y" ] ; then

    cd ../../montepython
    rm -r ../Euclid_w0wa/results/montepython_fisher/$PROBE/${CASE}_HP
    cp montepython/likelihoods/$LKL/$LKL.data.$CASE montepython/likelihoods/$LKL/$LKL.data
    rm data/euclid_xc_fiducial.dat
    $PYTHON montepython/MontePython.py run -p ../Euclid_w0wa/input/montepython_fisher/$PROBE/$CASE/${PROBE}_${CASE_SHORT}_HP.param -o ../Euclid_w0wa/results/montepython_fisher/$PROBE/${CASE}_HP -f 0
    $PYTHON montepython/MontePython.py run -o ../Euclid_w0wa/results/montepython_fisher/$PROBE/${CASE}_HP --fisher --fisher-step-it 1 --fisher-tol 10000
    cd ../Euclid_w0wa/test_precision

fi

echo "For photometric/optimistic:"
echo "For MontePyhton Fisher CLASS DP vs HP:"
echo "Shall we recompute error ratios and redo error plots ? (y/n)"
read answer
if [ "$answer" = "y" ] ; then

    # plot error comparisons
    cd ../plots
    #
    $PYTHON $PROBE/$CASE/CF_MP_DP-vs-HP.py --error-only
    cd ../test_precision

fi

####################

PROBE=spectroscopic
PROBE_SHORT=spectro
LKL=euclid_spectroscopic

CASE=optimistic
CASE_SHORT=opt

echo "For spectroscopic/optimistic:"
echo "Shall we run CosmicFish CLASS external DP? (y/n)"
read answer
if [ "$answer" = "y" ] ; then
    $PYTHON ../input/cosmicfish/$PROBE/$CASE/class_external_DP.py
fi

echo "For spectroscopic/optimistic:"
echo "Shall we run CosmicFish CLASS external HP? (y/n)"
read answer
if [ "$answer" = "y" ] ; then
    $PYTHON ../input/cosmicfish/$PROBE/$CASE/class_external_HP.py
fi

echo "For spectroscopic/optimistic:"
echo "For CosmicFish/external/CLASS DP vs HP:"
echo "Shall we recompute error ratios and redo error plots ? (y/n)"
read answer
if [ "$answer" = "y" ] ; then

    # plot error comparisons
    cd ../plots
    #
    $PYTHON $PROBE/$CASE/CF_class_ext_DP-vs-HP.py --error-only
    cd ../test_precision

fi

##################

echo "For spectroscopic/optimistic:"
echo "Shall we run CosmicFish CAMB external P2? (y/n)"
read answer
if [ "$answer" = "y" ] ; then
    $PYTHON ../input/cosmicfish/$PROBE/$CASE/camb_external_P2.py
fi

echo "For spectroscopic/optimistic:"
echo "Shall we run CosmicFish CAMB external P3? (y/n)"
read answer
if [ "$answer" = "y" ] ; then
    $PYTHON ../input/cosmicfish/$PROBE/$CASE/camb_external_P3.py
fi

echo "For spectroscopic/optimistic:"
echo "For CosmicFish/external/CAMB P2 vs P3:"
echo "Shall we recompute error ratios and redo error plots ? (y/n)"
read answer
if [ "$answer" = "y" ] ; then

    # plot error comparisons
    cd ../plots
    #
    $PYTHON $PROBE/$CASE/CF_camb_ext_P2-vs-P3.py --error-only
    rm CF_camb_ext_P2-vs-P3_cosmo_and_nuisance_matrix_ratio.pdf
    cd ../test_precision

fi

####################

echo "For spectroscopic/optimistic:"
echo "Shall we run MontePython Fisher CLASS DP ? (y/n)"
read answer
if [ "$answer" = "y" ] ; then

    cd ../../montepython
    rm -r ../Euclid_w0wa/results/montepython_fisher/$PROBE/${CASE}_DP
    cp montepython/likelihoods/$LKL/$LKL.data.$CASE montepython/likelihoods/$LKL/$LKL.data
    rm data/euclid_xc_fiducial.dat
    $PYTHON montepython/MontePython.py run -p ../Euclid_w0wa/input/montepython_fisher/$PROBE/$CASE/${PROBE}_${CASE_SHORT}_DP.param -o ../Euclid_w0wa/results/montepython_fisher/$PROBE/${CASE}_DP -f 0
    $PYTHON montepython/MontePython.py run -o ../Euclid_w0wa/results/montepython_fisher/$PROBE/${CASE}_DP --fisher --fisher-step-it 1 --fisher-tol 10000
    cd ../Euclid_w0wa/test_precision

fi

echo "For photometric/optimistic:"
echo "Shall we run MontePython Fisher CLASS HP ? (y/n)"
read answer
if [ "$answer" = "y" ] ; then

    cd ../../montepython
    rm -r ../Euclid_w0wa/results/montepython_fisher/$PROBE/${CASE}_HP
    cp montepython/likelihoods/$LKL/$LKL.data.$CASE montepython/likelihoods/$LKL/$LKL.data
    rm data/euclid_xc_fiducial.dat
    $PYTHON montepython/MontePython.py run -p ../Euclid_w0wa/input/montepython_fisher/$PROBE/$CASE/${PROBE}_${CASE_SHORT}_HP.param -o ../Euclid_w0wa/results/montepython_fisher/$PROBE/${CASE}_HP -f 0
    $PYTHON montepython/MontePython.py run -o ../Euclid_w0wa/results/montepython_fisher/$PROBE/${CASE}_HP --fisher --fisher-step-it 1 --fisher-tol 10000
    cd ../Euclid_w0wa/test_precision

fi

echo "For spectroscopic/optimistic:"
echo "For MontePyhton Fisher CLASS DP vs HP:"
echo "Shall we recompute error ratios and redo error plots ? (y/n)"
read answer
if [ "$answer" = "y" ] ; then

    # plot error comparisons
    cd ../plots
    #
    $PYTHON $PROBE/$CASE/CF_MP_DP-vs-HP.py --error-only
    cd ../test_precision

fi

########################
