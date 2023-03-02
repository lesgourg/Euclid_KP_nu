# Euclid_KPnu

This directory contains everything to reproduce forecasts for the KP_nu project of WP3.

To reproduce the results, the first the first step is to clone here adequate versions of the various codes, following instructions int he respective .md files

Then, the files in

    input/

show how to run the code (input files, scripts), while the files in

    results/

show the output that we obtained and that you should be able to reproduce.

    plots/

contains comparison plots.

To create the default conda environnement, execute once:

    conda env create -f environment-pip.yml

and then, to activate it, at the beginning of each new session,

    conda activate w0waMP
