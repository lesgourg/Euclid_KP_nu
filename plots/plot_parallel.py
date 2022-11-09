from glob import glob
import os
import argparse
os.chdir(os.path.dirname(os.path.realpath(__file__)))

parser = argparse.ArgumentParser(formatter_class=argparse.RawTextHelpFormatter)

parser.add_argument('-d', '--dir', metavar='\b',
                    help="    Specify the path to the directory of the probe",
                    default=os.path.join('.', '*', '*'))
parser.add_argument('-c1', '--code1', metavar='\b',
                    help="    Name of the first code", default='')
parser.add_argument('-c2', '--code2', metavar='\b',
                    help="    Name of the second code"
                    "\n\nCode combination names that can be passed for -c1"
                    " and -c2 :"
                    "\n1) It can be MP or class or class_ext or class_int or \
                    camb or camb_int or camb_ext "
                    "\n2) Just class or camb corresponds to both external and \
                    internal modes\n ", default='')
parser.add_argument('--error-only', action='store_true', dest='error_only',
                    help='    Plot error comparions plots only',
                    default=False)

args = parser.parse_args()
c1 = args.code1
c2 = args.code2
error_only = args.error_only


def checker(c1, c2):
    if c1 == '':
        if c2 == '':
            wildcard = '*.py'
        else:
            wildcard = '*' + c2 + '*.py'

    elif c2 == '':
        if c1 == '':
            wildcard = '*.py'
        else:
            wildcard = '*' + c1 + '*.py'

    else:
        wildcard = '*' + c1 + '*' + c2 + '*.py'
    return wildcard


def mainfxn(filepath):
    print('File exists !')
    var = ''
    for file in glob(filepath):
        var += 'python '
        var += file
        if error_only:
            var += ' --error-only'
        else : pass
        if file == glob(filepath)[-1] : pass
        else :
            var += ' & '

    print(var)
    choice = input('\n Do you want to execute the above command ? (y/n)')
    if choice == 'y' or choice == 'Y':
        os.system(var)


wildcard1 = checker(c1, c2)
filepath1 = os.path.join(args.dir, wildcard1)
wildcard2 = checker(c2, c1)
filepath2 = os.path.join(args.dir, wildcard2)

if glob(filepath1):
    mainfxn(filepath1)

elif glob(filepath2):
    mainfxn(filepath2)
else:
    print('No matches found')
