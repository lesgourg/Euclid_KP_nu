import sys
import os
from glob import glob
import argparse
import yaml
FILEPATH = os.path.abspath(os.path.dirname(os.path.realpath(__file__)))

# Function to extract data from the log.param file
def paramdatafxn(outdir):
    logfile = os.path.join(outdir, 'log.param')
    namelist = []
    with open(logfile) as f:
        lines = f.readlines()
    paramdata = dict()
    for line in lines:
        li = line.strip()
        if not li.startswith('#'):
            li = li.split('=')
            if len(li) > 1:
                name = li[0]
                data = li[1]
                data = data.strip()
                if data.startswith("["):
                    data = eval(data)
                    if len(data) > 1:
                        name = name.split("[")[1].strip()[
                            :-1].strip().replace('\'', '').replace('\"', '')
                        if data[3] > 0:
                            paramdata[name] = data
                            namelist.append(name)
    return paramdata, namelist

# The below fxn is not used anywhere, it can be used at some point if needed.
# It creates a dictionary { MP parameter names : MP param Latex names }


def MPparanamnamesfxn(outdir):
    paramnamespath = max(
        glob(os.path.join(outdir, '[!fish]*'+'.paramnames')), key=os.path.getctime)
    paramnames = dict()
    with open(paramnamespath) as f:
        lines2 = f.readlines()
    for line in lines2:
        li = line.strip()
        if not li.startswith('#'):
            name = li.split("\t")[0].strip()
            latexname = li.split("\t")[1].strip()
            paramnames[name] = latexname
    return paramnames

# Main file that is used to create the fisher.paramnames file from the log.param file


def write_files(outdir, out_filename):
    paramdata, namelist = paramdatafxn(outdir)
    with open(os.path.join(FILEPATH,'mp_to_cf.yaml'), 'r') as yamlfile:
        paramnames = yaml.load(yamlfile, Loader=yaml.FullLoader)
    newparamfile = os.path.join(outdir, out_filename)
    print('Reading log.param from the directory : ', outdir)
    print('Writing a new file named : ', newparamfile)
    print('Success ! ')
    with open(newparamfile, 'w') as f:
        f.write(
            '#Paramnames file for CosmicFish, for a Fisher that was produced by MP\n')
        for i in namelist:
            f.write(str(paramnames[i][0]) + '    ' + str(paramnames[i]
                    [1]) + '    ' + str(paramdata[i][0]) + '\n')


parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)
parser.add_argument(
    dest='outdir', help="Specify the path to the output directory of MP")
args = parser.parse_args()
mpoutdir = args.outdir


write_files(outdir=mpoutdir, out_filename='fisher.paramnames')
