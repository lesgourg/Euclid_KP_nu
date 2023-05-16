For this task you can clone the master branch of camb in a directory CAMB/ on
the same level as this one. It should be equivalent to use the public repository:

    cd ..
    git clone --recursive https://github.com/cmbant/CAMB

or the version of Santiago (which differs only for things relevant for input_4_cast, not used in KP_nu):

    cd ..
    git clone --recursive https://github.com/santiagocasas/CAMB

then:

    cd CAMB
    python3 setup.py install [--user]
    cd ../../Euclid_KP_nu