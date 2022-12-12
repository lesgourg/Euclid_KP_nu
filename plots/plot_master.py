import os, sys
from re import A
import inspect
from pathlib import Path
import seaborn as sns

snscolors=sns.color_palette("Set1")

# dict must be of the type {'code' : 'cf', 'specs' : 'pessimistic', 'mode' : 'internal'
# , 'boltzmann' : 'class'}
def fisher_path(dict) :
        os.chdir(os.path.dirname(os.path.realpath(__file__)))
        derivatives_default = {'gcsp' : '3PT' , 'wl' : '3PT','wlxgcph' : '3PT','gcph' : '3PT'}
        probe_path = {'wl':'wl_only','wlxgcph':'photometric','gcsp':'spectroscopic'}
        resultsdir = '../results'

        if dict['code'].lower() == 'montepython' or dict['code'].lower() == 'mp' :
                probe_dir = probe_path[dict['probe'].lower()]
                specs_dir = dict['specs'].lower()
                path = os.path.join(resultsdir,'montepython_fisher',probe_dir,specs_dir)
                path = os.path.join(path,'fisher.mat')

        elif 'cosmicfish' in dict['code'].lower() or dict['code'].lower() == 'cf' :
                probe_filename = {'wlxgcph':'WLGCph','gcsp':'GCsp','wl':'WL'}
                probe_dir = probe_path[dict['probe'].lower()]
                specs_dir = dict['specs'].lower()
                mode = dict['mode'].lower()
                if mode == 'external' : precision = dict['precision'].upper() + '_'
                else : precision = ''
                path = os.path.join(resultsdir,probe_dir,specs_dir)
                code = dict['boltzmann'].lower()
                specs = (dict['specs'].lower()).capitalize()
                if 'derivative' in dict :
                        dr = dict['derivative']
                else :
                        dr = derivatives_default[dict['probe'].lower()]

                filename = 'CosmicFish_v1.0_w0wa_'+mode+'_'+code+'-'+ specs +'-'+ dr +'_' + \
                          precision + probe_filename[dict['probe'].lower()] +'_'+'fishermatrix.txt'
                path = os.path.join(resultsdir,'cosmicfish'+'_'+mode,probe_dir,specs_dir,filename)
        return os.path.realpath(path)


def plotter(fish_files, labels, pars, outpath='automatic',
            script_name='automatic', error_only=False, compare_errors_dict=dict()):

    fish_files = [os.path.abspath(i) for i in fish_files] ## This is evaluated at old CWD
    os.chdir(os.path.dirname(os.path.realpath(__file__))) ## CWD changes
    sys.path.append('../../cosmicfish_reloaded/')
    from cosmicfishpie.analysis import fisher_plot_analysis as cfa
    from cosmicfishpie.analysis import fisher_plotting as cfp
    from cosmicfishpie.analysis import fisher_matrix as cfm
    print('\n\n')
    print('----> Fishers are located at ----> \n')
    for fi in fish_files:
        print(fi,'\n')
    #print(fish_files[1],'\n')
    if outpath == 'automatic' :
        outpath = inspect.stack()[1][1]
        outpath = Path(outpath).parent
        outpath = str(outpath) +'/'
        print('Output path not provided ! defaulting to ',outpath)

    else : pass
    if script_name == 'automatic':
        script_name = inspect.stack()[1][1]
        script_name = os.path.splitext(os.path.basename(os.path.realpath(script_name)))[0]
        script_name = str(script_name)
        print(script_name)
    else : pass

    cosmo_names = ['Omegam', 'Omegab', 'ns', 'h','sigma8','mnu','Neff']
    nuisance_names = list( set(pars) - set(cosmo_names)  )
    fgroup=cfa.CosmicFish_FisherAnalysis()
    for fii in fish_files:
        ftesta = cfm.fisher_matrix(file_name=fii)
        orignames=ftesta.get_param_names()
        print(orignames)
        ftemp=cfm.fisher_matrix(file_name=fii)
        print(ftemp.name)
        fgroup.add_fisher_matrix(ftemp)

#     #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     COSMO + Nuisance      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    cutnames=pars

    fgroupRe = fgroup.reshuffle(params=cutnames, update_names=False)


    pessions = {'fishers_group': fgroupRe,
            'colors': snscolors,
            'dpi':100,
            'fish_labels': labels,
            'plot_pars': pars,
            'plot_method': 'Gaussian',
            'axis_custom_factors' : {'all':4},
            'outpath': outpath,
            'outroot': script_name + '_' + 'cosmo_and_nuisance',
            'param_labels' : pars
            }
    compare_errors_dict_opts = {'save_error':True}
    compare_errors_dict_opts.update(compare_errors_dict)
    fish_plotter = cfp.fisher_plotting(**pessions)
    if error_only :
        fish_plotter.compare_errors(compare_errors_dict_opts)
    else :
        fish_plotter.load_gaussians()
        fish_plotter.plot_fisher(filled=True)
        fish_plotter.compare_errors(compare_errors_dict_opts)
        #fish_plotter.matrix_ratio()

    return

#     #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    COSMO marginalizing Nuisance   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cutnames=pars

    fgroupRe = fgroup.reshuffle(params=cutnames, update_names=False)


    pessions = {'fishers_group': fgroupRe,
            'colors': snscolors,
            'dpi':100,
            'fish_labels': labels,
            'plot_pars': cosmo_names,
            'plot_method': 'Gaussian',
            'axis_custom_factors' : {'all':4},
            'outpath': outpath,
            'outroot': script_name + '_' + 'cosmo_marg_nuisance'
            }

    fish_plotter = cfp.fisher_plotting(**pessions)
    if error_only :
        fish_plotter.compare_errors(compare_errors_dict_opts)
    else :
        fish_plotter.load_gaussians()
        fish_plotter.plot_fisher(filled=True)
        fish_plotter.compare_errors(compare_errors_dict_opts)



#     #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    COSMO fixing Nuisance   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cutnames=cosmo_names

    fgroupRe = fgroup.reshuffle(params=cutnames, update_names=False)


    pessions = {'fishers_group': fgroupRe,
            'colors': snscolors,
            'dpi':100,
            'fish_labels': labels,
            'plot_pars': cosmo_names,
            'plot_method': 'Gaussian',
            'axis_custom_factors' : {'all':4},
            'outpath': outpath,
            'outroot': script_name+ '_' + 'cosmo_fix_nuisance'
            }

    fish_plotter = cfp.fisher_plotting(**pessions)
    if error_only :
        fish_plotter.compare_errors(compare_errors_dict_opts)
    else :
        fish_plotter.load_gaussians()
        fish_plotter.plot_fisher(filled=True)
        fish_plotter.compare_errors(compare_errors_dict_opts)

    #%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     Nuisance fixing cosmo     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cutnames=nuisance_names

    fgroupRe = fgroup.reshuffle(params=cutnames, update_names=False)


    pessions = {'fishers_group': fgroupRe,
            'colors': snscolors,
            'dpi':100,
            'fish_labels': labels,
            'plot_pars': nuisance_names,
            'plot_method': 'Gaussian',
            'axis_custom_factors' : {'all':4},
            'outpath': outpath,
            'outroot':script_name+ '_' +'nuisance_fix_cosmo'
            }

    fish_plotter = cfp.fisher_plotting(**pessions)
    if error_only :
        fish_plotter.compare_errors(compare_errors_dict_opts)
    else :
        fish_plotter.load_gaussians()
        fish_plotter.plot_fisher(filled=True)
        fish_plotter.compare_errors(compare_errors_dict_opts)
