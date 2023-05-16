For this task you can clone the class_nonlinear branch of class in a directory class/ on the same level as this one from the private repository:

    cd ..
    git clone https://github.com/lesgourg/class.git class/

then:

    cd class
    git checkout class_nonlinear
    make clean
    PYTHON=python3 make -j
    cd ../Euclid_KP_nu