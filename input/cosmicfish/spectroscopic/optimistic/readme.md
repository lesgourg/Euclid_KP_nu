To run CosmiFish in internal or external mode, you don't need to go to
the CosmicFish executable directory, you can run from here. The input
file is a python script, e.g. class_internal.py. You can edit it to
specify which case you want to run. This file only contains the
settings that differ from run to run. Those settings that are supposed
to be common to all the runs of this repository are in a common input
file located in

    results/cosmifish/internal_master.py
    results/cosmifish/external_master.py


You can run direcly with e.g.

    python3 class_external.py

Results are automatically written in the results/directory, with a naming of the output directory and files choosen automatically by CosmicFish as a function of the input.
