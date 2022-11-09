# Euclid_w0wa

This directory contains everything to reproduce forecasts for LCDM+w0+wa.

To reproduce the results, the first the first step is to clone here adequate versions of the various codes, following instructions int he respective .md files

Then, the files in

    input/

show how to run the code (input files, scripts), while the files in

    results/

show the output that we obtained and that you should be able to reproduce.

    plots/

contains comparison plots.

Final table comparing all codes in google sheets : [Google sheets link](https://docs.google.com/spreadsheets/d/1DcKexttMQ_1JmajdQiL0LbSD3m1FXQKU/edit?usp=sharing&ouid=116025791120069747063&rtpof=true&sd=true)

To create the default conda environnement, execute once:

    conda env create -f environment-pip.yml

and then, to activate it, at the beginning of each new session,

    conda activate w0waMP