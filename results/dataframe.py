import numpy as np
from matplotlib import cm,colors
from matplotlib.colors import ListedColormap
from itertools import combinations as cb


types = ['All Pars', 'Cosmo with fixed Nuisance', 'Cosmo with marginalized Nuisance', 'Nuisance with fixed Nuisance'] 
codes = ['Class EXT','Class INT','MP','Camb EXT' , 'Camb INT']
types_dict = {'All Pars':'cosmo_and_nuisance', 'Cosmo with fixed Nuisance' : 'cosmo_fix_nuisance', 
'Cosmo with marginalized Nuisance' : 'cosmo_marg_nuisance', 'Nuisance with fixed Nuisance' :  'nuisance_fix_cosmo'  }
name_to_cell = {'Class EXT':'class_ext','Class INT' : 'class_int', 'MP':'MP',
'Camb EXT' : 'camb_ext','Camb INT' : 'camb_int'}
all_combinations = list(cb(codes,2))

colormap = cm.get_cmap('RdYlGn')(np.linspace(0,1,500))[::-1]
gr = colormap[0:120]
ora = colormap[370:400]
rd = colormap[400:500]
nmap = np.concatenate((gr,ora,rd),axis=0)

greenmap = ListedColormap(gr)
orangemap = ListedColormap(ora)
redmap = ListedColormap(rd)
newmap = ListedColormap(nmap)
norm_green = colors.Normalize(vmin=0, vmax=10, clip=True)
mapper_green = cm.ScalarMappable(norm=norm_green, cmap=greenmap)
norm_orange = colors.Normalize(vmin=10, vmax=20, clip=True)
mapper_orange = cm.ScalarMappable(norm=norm_orange, cmap=orangemap)
norm_red = colors.Normalize(vmin=30, vmax=70, clip=True)
mapper_red = cm.ScalarMappable(norm=norm_red, cmap=redmap)


def _hexcolor(value,mapper):

    r,g,b,a = mapper.to_rgba(value)
    r = int(r*255)
    g = int(g*255)
    b = int(b*255)

    return "#%s%s%s" % tuple([hex(c)[2:].rjust(2, "0") for c in (r, g, b)])
 
def color_map(v):

    try :
        if v >0 and v <= 10.0 :
            color = _hexcolor(v,mapper_green)
            return f"background-color: " + color
        elif v>10.0 and v <= 30.0 :
            color=_hexcolor(v,mapper_orange)
            return f"background-color: " + color
        elif v>30.0 and v<99.0 :
            color=_hexcolor(v,mapper_red)
            return f"background-color: " + color
        else :
            return None
    except :
        return None


def write_data(df,category,files) :
   for cell1,cell2 in all_combinations :
    name1 = name_to_cell[cell1]
    name2 = name_to_cell[cell2]
    for file in files :
        if name1 in file and name2 in file and category in file :
            contents = np.loadtxt(file)
            max_val = np.amax(np.abs(contents))
            df[cell1][cell2] = max_val
            df[cell2][cell1] = max_val