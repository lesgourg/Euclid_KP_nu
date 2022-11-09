Go to the input_4_cast directory to run the code, and check that you have correctly checked out to the branch wp3, not the main branch

    cd ../input_4_cast
    git status

----------------------------

* For producing the CLASS external files, you will use the input files

        ../Euclid_w0wa/input/input_4_cast/class_w0wa_LP.ini
        ../Euclid_w0wa/input/input_4_cast/class_w0wa_MP.ini
        ../Euclid_w0wa/input/input_4_cast/class_w0wa_HP.ini

    (for respectively low precision, medium precision, high precision).
    The relative paths at the beginning may need some editing if, on your machine, class is not located in the directory ../class/

    Then run with e.g.

        python3 run.py ../Euclid_w0wa/input/input_4_cast/class_w0wa_LP.ini

----------------------------

* For producing the CAMB external files, you will use the input file

        ../Euclid_w0wa/input/input_4_cast/camb_w0wa_LP.ini
        ../Euclid_w0wa/input/input_4_cast/camb_w0wa_MP.ini
        ../Euclid_w0wa/input/input_4_cast/camb_w0wa_HP.ini

    (for respectively low precision, medium precision, high precision).
    The relative paths at the beginning may need some editing if, on your machine, class is not located in the directory ../CAMB/

    Also, set number_of_threads to the number of cores in your machine

    Then run with e.g.

        python3 run.py ../Euclid_w0wa/input/input_4_cast/camb_w0wa_HP.ini

-----------------------------

* For comparing CLASS and CAMB files, run e.g.

    python3 compare.py --save_plots --threshold 0.01 output/camb_w0wa_HP output/class_w0wa_HP

    All spectra should agree very well (typically better than 0.1% for
    the linear spectrum with MP or HP)

-----------------------------

When everything is done, go back

     cd ../Euclid_w0wa