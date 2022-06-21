For this project, at the moment, you need to clone montepython from the private repository (in a directory montepython at the same level as this one) and checkout the branch euclid_cl_xc:

    cd ..
    git clone https://github.com/lesgourg/montepython montepython
    cd montepython
    git co euclid_cl_xc

To set up the path to class you need to do:

    cp default.conf.template default.conf

the edit default.conf and set the path "root" to the directory that contains class/, montepython/, Euclid_KP_nu/, etc.

Return here with

    cd ../Euclid_KP_nu