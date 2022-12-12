PYTHON=python3

cd ../../montepython

PROBE=photometric
PROBE_SHORT=photo
PROBE_FIDUCIAL=euclid_xc_fiducial.dat
CASE=pessimistic
CASE_SHORT=pess
PRECISION=DP

### IMPORTANT NOTE
#
# These runs only work in combination with what is currently a local modification of Montepyhton.
# In sampler.py, there is a new function get_fisher_vs_stepsize(), outputting either
# one non-diagonal coefficient or two diagonal coefficients of the Fisher matrix.
# For each of the two, you would need to run the script below to go through all probes/cases/precisions.
# However, the SP-MP paper only contains the figure for the non-diagonal coefficient.

echo "Do $PROBE_SHORT $CASE_SHORT $PRECISION [1/8]? (y/n)"
read answer
if [ "$answer" = "y" ] ; then
    rm data/$PROBE_FIDUCIAL
    rm -rf ../Euclid_w0wa/test_stepsizes/$PROBE/{$CASE}_$PRECISION
    python3 montepython/MontePython.py run -p ../Euclid_w0wa/test_stepsizes/${PROBE}_${CASE_SHORT}_$PRECISION.param -o ../Euclid_w0wa/test_stepsizes/$PROBE/${CASE}_$PRECISION -f 0
    python3 montepython/MontePython.py run -o ../Euclid_w0wa/test_stepsizes/$PROBE/${CASE}_$PRECISION --method FvS
fi

PROBE=photometric
PROBE_SHORT=photo
PROBE_FIDUCIAL=euclid_xc_fiducial.dat
CASE=pessimistic
CASE_SHORT=pess
PRECISION=HP

echo "Do $PROBE_SHORT $CASE_SHORT $PRECISION [2/8]? (y/n)"
read answer
if [ "$answer" = "y" ] ; then
    rm data/$PROBE_FIDUCIAL
    rm -rf ../Euclid_w0wa/test_stepsizes/$PROBE/{$CASE}_$PRECISION
    python3 montepython/MontePython.py run -p ../Euclid_w0wa/test_stepsizes/${PROBE}_${CASE_SHORT}_$PRECISION.param -o ../Euclid_w0wa/test_stepsizes/$PROBE/${CASE}_$PRECISION -f 0
    python3 montepython/MontePython.py run -o ../Euclid_w0wa/test_stepsizes/$PROBE/${CASE}_$PRECISION --method FvS
fi

PROBE=photometric
PROBE_SHORT=photo
PROBE_FIDUCIAL=euclid_xc_fiducial.dat
CASE=optimistic
CASE_SHORT=opt
PRECISION=DP

echo "Do $PROBE_SHORT $CASE_SHORT $PRECISION [3/8]? (y/n)"
read answer
if [ "$answer" = "y" ] ; then
    rm data/$PROBE_FIDUCIAL
    rm -rf ../Euclid_w0wa/test_stepsizes/$PROBE/{$CASE}_$PRECISION
    python3 montepython/MontePython.py run -p ../Euclid_w0wa/test_stepsizes/${PROBE}_${CASE_SHORT}_$PRECISION.param -o ../Euclid_w0wa/test_stepsizes/$PROBE/${CASE}_$PRECISION -f 0
    python3 montepython/MontePython.py run -o ../Euclid_w0wa/test_stepsizes/$PROBE/${CASE}_$PRECISION --method FvS
fi

PROBE=photometric
PROBE_SHORT=photo
PROBE_FIDUCIAL=euclid_xc_fiducial.dat
CASE=optimistic
CASE_SHORT=opt
PRECISION=HP

echo "Do $PROBE_SHORT $CASE_SHORT $PRECISION [4/8]? (y/n)"
read answer
if [ "$answer" = "y" ] ; then
    rm data/$PROBE_FIDUCIAL
    rm -rf ../Euclid_w0wa/test_stepsizes/$PROBE/{$CASE}_$PRECISION
    python3 montepython/MontePython.py run -p ../Euclid_w0wa/test_stepsizes/${PROBE}_${CASE_SHORT}_$PRECISION.param -o ../Euclid_w0wa/test_stepsizes/$PROBE/${CASE}_$PRECISION -f 0
    python3 montepython/MontePython.py run -o ../Euclid_w0wa/test_stepsizes/$PROBE/${CASE}_$PRECISION --method FvS
fi

PROBE=spectroscopic
PROBE_SHORT=spectro
PROBE_FIDUCIAL=euclid_pk_fiducial.dat
CASE=pessimistic
CASE_SHORT=pess
PRECISION=DP

echo "Do $PROBE_SHORT $CASE_SHORT $PRECISION [5/8]? (y/n)"
read answer
if [ "$answer" = "y" ] ; then
    rm data/$PROBE_FIDUCIAL
    rm -rf ../Euclid_w0wa/test_stepsizes/$PROBE/{$CASE}_$PRECISION
    python3 montepython/MontePython.py run -p ../Euclid_w0wa/test_stepsizes/${PROBE}_${CASE_SHORT}_$PRECISION.param -o ../Euclid_w0wa/test_stepsizes/$PROBE/${CASE}_$PRECISION -f 0
    python3 montepython/MontePython.py run -o ../Euclid_w0wa/test_stepsizes/$PROBE/${CASE}_$PRECISION --method FvS
fi

PROBE=spectroscopic
PROBE_SHORT=spectro
PROBE_FIDUCIAL=euclid_pk_fiducial.dat
CASE=pessimistic
CASE_SHORT=pess
PRECISION=HP

echo "Do $PROBE_SHORT $CASE_SHORT $PRECISION [6/8]? (y/n)"
read answer
if [ "$answer" = "y" ] ; then
    rm data/$PROBE_FIDUCIAL
    rm -rf ../Euclid_w0wa/test_stepsizes/$PROBE/{$CASE}_$PRECISION
    python3 montepython/MontePython.py run -p ../Euclid_w0wa/test_stepsizes/${PROBE}_${CASE_SHORT}_$PRECISION.param -o ../Euclid_w0wa/test_stepsizes/$PROBE/${CASE}_$PRECISION -f 0
    python3 montepython/MontePython.py run -o ../Euclid_w0wa/test_stepsizes/$PROBE/${CASE}_$PRECISION --method FvS
fi

PROBE=spectroscopic
PROBE_SHORT=spectro
PROBE_FIDUCIAL=euclid_pk_fiducial.dat
CASE=optimistic
CASE_SHORT=opt
PRECISION=DP

echo "Do $PROBE_SHORT $CASE_SHORT $PRECISION [7/8]? (y/n)"
read answer
if [ "$answer" = "y" ] ; then
    rm data/$PROBE_FIDUCIAL
    rm -rf ../Euclid_w0wa/test_stepsizes/$PROBE/{$CASE}_$PRECISION
    python3 montepython/MontePython.py run -p ../Euclid_w0wa/test_stepsizes/${PROBE}_${CASE_SHORT}_$PRECISION.param -o ../Euclid_w0wa/test_stepsizes/$PROBE/${CASE}_$PRECISION -f 0
    python3 montepython/MontePython.py run -o ../Euclid_w0wa/test_stepsizes/$PROBE/${CASE}_$PRECISION --method FvS
fi

PROBE=spectroscopic
PROBE_SHORT=spectro
PROBE_FIDUCIAL=euclid_pk_fiducial.dat
CASE=optimistic
CASE_SHORT=opt
PRECISION=HP

echo "Do $PROBE_SHORT $CASE_SHORT $PRECISION [8/8]? (y/n)"
read answer
if [ "$answer" = "y" ] ; then
    rm data/$PROBE_FIDUCIAL
    rm -rf ../Euclid_w0wa/test_stepsizes/$PROBE/{$CASE}_$PRECISION
    python3 montepython/MontePython.py run -p ../Euclid_w0wa/test_stepsizes/${PROBE}_${CASE_SHORT}_$PRECISION.param -o ../Euclid_w0wa/test_stepsizes/$PROBE/${CASE}_$PRECISION -f 0
    python3 montepython/MontePython.py run -o ../Euclid_w0wa/test_stepsizes/$PROBE/${CASE}_$PRECISION --method FvS
fi

cd ../Euclid_w0wa/test_stepsizes
