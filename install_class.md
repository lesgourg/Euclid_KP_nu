For this task you can clone the master branch of class in a directory class/ on the same level as this one, either from the public repository:

    cd ..
    git clone https://github.com/lesgourg/class_public.git class/

or the private one:

    cd ..
    git clone https://github.com/lesgourg/class.git class/

then:

    cd class
    make clean
    PYTHON=python3 make -j
    cd ../Euclid_KP_nu