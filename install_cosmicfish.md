For this project you should use (for the moment) the branch classtest of cosmicfish. You should clone it in a directory cosmicfish_reloaded/ on the same level as this one

    cd ..
    git clone https://gitlab.com/matmartinelli/cosmicfish_reloaded.git
    cd cosmicfish_reloaded
    git checkout classtest

You can install the packages that are required using pip, or you can create a conda environment. For installing the packages use,

    pip install -r requirements.txt

Or to create a conda environment use,

    conda env create -f environment.yaml

Then return to this directory:

    cd ../Euclid_KP_nu

(Remark: In a given scripts example.py, inside the dictionary `options`, the absolute or relative path (relative to example.py) to 'survey specifications directory must be provided with the key `specs_dir` and similarly for the boltzmann yaml files the paths must be specified in the key `class_config_yaml` and `camb_config_yaml`)
