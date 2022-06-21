# Plots

 * This directory contains all the code comparison plots along with the plotting routines.

 * The specific versions used for each plot are in the sub-directories, e.g. photometric/optimistic/...

 * Each such sub-directory contains the python plotting file along with the plots, for example,
    `./photometric/optimistic/CF_camb_int-vs-MP.py` is the python file that needs to be executed to run the plotting routine for the comparison between CosmicFish in internal mode with Camb and MontePython.

# Plotting many code comparisons in parallel 

 Many such plotting files can be launched in parallel by using the ampersand `&` symbol. This can be automated using the python script   `plot_parallel.py`. For example, to run all plotting routines involving for the probe photometric-optimistic the following command has to be used.

    python plot_parallel.py -d ./photometric/optimistic/
        
# Arguments of the python file 

 *  `-d` or `--dir` :  Specify the probe by giving the path to the folder containing the plotting scripts of that probe, if not passed the routine searches for all possible probe combinations
 
 *  `-c1` or `--code1` : To specify the or all comparisons involving a given code combination like `camb_ext` then use the argument , which specifies the first code.

 *  `-c2` or `--code2` : Same as above but specifies the second code for comparison

 *  The code combinations can be specified using the following names :  `camb, camb_ext, camb_int, class, class_ext, class_int, MP`.
    Just `camb`, `class` requests the code to plot both internal and external modes
 
 *  The help message can be triggered using the `-h` or `--help`