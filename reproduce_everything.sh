PYTHON=python3

# NOTE: the IST:Fisher matrix located in the ../fisher_for_public directory
# can only be handled after a renaming of the parameters, performed by running
# once the commands:
# python external_fishers/replace_paramnames.py ../fisher_for_public/All_Results/pessimistic/flat/
# python external_fishers/replace_paramnames.py ../fisher_for_public/All_Results/optimistic/flat/

# Select here the probe (photometric/spectroscopic)
#
PROBE=photometric
PROBE_SHORT=photo
LKL=euclid_photometric_z
#
#PROBE=spectroscopic
#PROBE_SHORT=spec
#LKL=euclid_spectroscopic

# Select here the case (pessimistic/optimistic)
#
#CASE=pessimistic
#CASE_SHORT=pess
#
CASE=optimistic
CASE_SHORT=opt

# Select here the precision
CLASS_PREC=HP
CAMB_PREC=P3

# placeholder for running optionally input_4_cast

echo "Shall we rerun the $PROBE $CASE case ? (y/n)"
read answer
if [ "$answer" = "y" ] ; then

    # run CosmicFish

    echo "Shall we rerun CosmicFish CLASS external ? (y/n)"
    read answer
    if [ "$answer" = "y" ] ; then
        $PYTHON input/cosmicfish/$PROBE/$CASE/class_external_${CLASS_PREC}.py
    fi

    echo "Shall we rerun CosmicFish CAMB external ? (y/n)"
    read answer
    if [ "$answer" = "y" ] ; then
        $PYTHON input/cosmicfish/$PROBE/$CASE/camb_external_${CAMB_PREC}.py
    fi

    echo "Shall we rerun CosmicFish CLASS internal ? (y/n)"
    read answer
    if [ "$answer" = "y" ] ; then
        $PYTHON input/cosmicfish/$PROBE/$CASE/class_internal.py
    fi

    echo "Shall we rerun CosmicFish CAMB internal ? (y/n)"
    read answer
    if [ "$answer" = "y" ] ; then
        $PYTHON input/cosmicfish/$PROBE/$CASE/camb_internal.py
    fi

    # run MP

    echo "Shall we rerun MontePython Fisher ? (y/n)"
    read answer
    if [ "$answer" = "y" ] ; then
        cd ../montepython
        cp montepython/likelihoods/$LKL/$LKL.data.$CASE montepython/likelihoods/$LKL/$LKL.data
        if [ "$PROBE" = "photometric" ] ; then
            rm data/euclid_xc_fiducial.dat
        fi
        if [ "$PROBE" = "spectroscopic" ] ; then
            rm data/euclid_pk_fiducial.dat
        fi
        rm -rf ../Euclid_w0wa/results/montepython_fisher/$PROBE/${CASE}_${CLASS_PREC}
        $PYTHON montepython/MontePython.py run -p ../Euclid_w0wa/input/montepython_fisher/$PROBE/$CASE/${PROBE}_${CASE_SHORT}_${CLASS_PREC}.param -o ../Euclid_w0wa/results/montepython_fisher/$PROBE/${CASE}_${CLASS_PREC} -f 0
        $PYTHON montepython/MontePython.py run -o ../Euclid_w0wa/results/montepython_fisher/$PROBE/${CASE}_${CLASS_PREC} --fisher --fisher-step-it 1 --fisher-tol 10000
        cd ../Euclid_w0wa
        $PYTHON input/montepython_fisher/paramnames_for_cosmicfish.py results/montepython_fisher/$PROBE/${CASE}_${CLASS_PREC}
    fi

    # TBD: run MCMC
    echo "Shall we rerun MontePython MCMC? (y/n)"
    read answer
    if [ "$answer" = "y" ] ; then

        echo "Shall we erase previous output directory with the chains? [WARNING: make a backup before saying yes] (y/n)"
        read answer
        if [ "$answer" = "y" ] ; then
            rm -r ../Euclid_w0wa/results/montepython_mcmc/w0wa_${PROBE_SHORT}_${CASE_SHORT}
            cd ../montepython
            cp montepython/likelihoods/$LKL/$LKL.data.$CASE montepython/likelihoods/$LKL/$LKL.data
            rm data/euclid_xc_fiducial.dat
            $PYTHON montepython/MontePython.py run -p ../Euclid_w0wa/input/montepython_fisher/$PROBE/$CASE/${PROBE}_$CASE_SHORT.param -o ../Euclid_w0wa/results/montepython_mcmc/w0wa_${PROBE_SHORT}_${CASE_SHORT} -f 0
            $PYTHON montepython/MontePython.py run -o ../Euclid_w0wa/results/montepython_mcmc/w0wa_${PROBE_SHORT}_${CASE_SHORT} -N 100000 --update 50 --superupdate 20 --covmat .... --conf default.conf
            cd ../Euclid_w0wa
        else
            echo "Will not erase old directory and will not rerun. Make a back up of the old chains and rerun the script."
        fi
    fi

    echo "Shall we recompute error ratios and redo error plots ? (y/n)"
    echo "(this step is necessary to update the main comparison table)"
    read answer
    if [ "$answer" = "y" ] ; then

        # plot error comparisons
        cd plots
        #
        # uncomment if you prefer old files to be removed rather than overwritten
        #rm $PROBE/$CASE/*.txt
        #rm $PROBE/$CASE*/.pdf
        #
        $PYTHON $PROBE/$CASE/CF_camb_ext-vs-MP.py --error-only
        rm $PROBE/$CASE/CF_camb_ext-vs-MP_cosmo_and_nuisance_error_comparison.pdf
        #
        $PYTHON $PROBE/$CASE/CF_camb_ext-vs-camb_int.py --error-only
        # keep error comparison plot for the paper
        #
        $PYTHON $PROBE/$CASE/CF_camb_int-vs-MP.py --error-only
        rm $PROBE/$CASE/CF_camb_int-vs-MP_cosmo_and_nuisance_error_comparison.pdf
        #
        $PYTHON $PROBE/$CASE/CF_class_ext-vs-MP.py --error-only
        rm $PROBE/$CASE/CF_class_ext-vs-MP_cosmo_and_nuisance_error_comparison.pdf
        #
        $PYTHON $PROBE/$CASE/CF_class_ext-vs-camb_ext.py --error-only
        rm $PROBE/$CASE/CF_class_ext-vs-camb_ext_cosmo_and_nuisance_error_comparison.pdf
        #
        $PYTHON $PROBE/$CASE/CF_class_ext-vs-camb_int.py --error-only
        rm $PROBE/$CASE/CF_class_ext-vs-camb_int_cosmo_and_nuisance_error_comparison.pdf
        #
        $PYTHON $PROBE/$CASE/CF_class_ext-vs-class_int.py --error-only
        # keep error comparison plot for the paper
        #
        $PYTHON $PROBE/$CASE/CF_class_int-vs-MP.py --error-only
        # keep error comparison plot for the paper
        #
        $PYTHON $PROBE/$CASE/CF_class_int-vs-camb_ext.py --error-only
        rm $PROBE/$CASE/CF_class_int-vs-camb_ext_cosmo_and_nuisance_error_comparison.pdf
        #
        $PYTHON $PROBE/$CASE/CF_class_int-vs-camb_int.py --error-only
        # keep error comparison plot for the paper
        #
        # run comparison plot between various pipelines
        $PYTHON $PROBE/$CASE/4codes-CF_class_camb-vs-ISTF-vs-MP.py --error-only
        # keep error comparison plot for the paper
        #
        cd ..
    fi

    echo "Shall we replot main the comparison table ? (y/n)"
    read answer
    if [ "$answer" = "y" ] ; then

        # run comparison table
        cd results/comparison_table
        $PYTHON table-$PROBE_SHORT-$CASE_SHORT.py
        cd ../..
    fi

    echo "Shall we replot main the contours with getdist ? (y/n)"
    read answer
    if [ "$answer" = "y" ] ; then

        # run getdist
        cd getdist
        $PYTHON ${PROBE_SHORT}_$CASE_SHORT.py
        cd ..
    fi

    echo "Shall we redo the LaTeX table with absolute errors ? (y/n)"
    read answer
    if [ "$answer" = "y" ] ; then

        if [ "$PROBE" = "photometric" ] ; then
           echo "setting the parameter names in fisher_for_public"
           python3 external_fishers/replace_paramnames.py ../fisher_for_public/All_Results/$CASE/flat/
        fi

        # run absolute error script
        cd results/absolute_sigmas
        $PYTHON $PROBE_SHORT-$CASE_SHORT-sigmas.py
        cd ../..
    fi

fi
