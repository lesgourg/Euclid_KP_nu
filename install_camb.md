For this task you can clone the master branch of camb in a directory CAMB/ on
the same level as this one. At some point it should become possible to use the public repository:

    cd ..
    git clone --recursive https://github.com/cmbant/CAMB

However, so far, a small bug makes input_4_cast incompatible with this version, and we use instead

    cd ..
    git clone --recursive https://github.com/santiagocasas/CAMB

then:

    cd CAMB
    python3 setup.py install [--user]
    cd ../../Euclid_KP_nu