
import numpy as np

from pylightcurve.errors import *
from pylightcurve.analysis.distributions import two_d_distribution

import matplotlib.cm as cm
import matplotlib.gridspec as gridspec
from matplotlib.figure import Figure
from matplotlib.backend_bases import FigureCanvasBase


# mcmc plots


def plot_mcmc_fitting(fitting_object, output_file):

    if not fitting_object.mcmc_run_complete:
        raise PyLCProcessError('MCMC not completed')

    results = fitting_object.results

    f1 = 3
    f2 = 1

    fbottom = 0.12
    ftop = 0.05
    fleft = 0.1
    fright = 0.05

    funit = 1.0
    frow = f1 + f2
    fcol = round(frow * 1.5)
    frow_height = (1.0 - fbottom - ftop) / frow

    fwidth = funit * fcol / (1 - fright - fleft)
    fheight = funit * frow / (1 - fbottom - ftop)

    lebels_right = fleft / 4

    fsmain = 2.0 * fheight
    fsbig = 1.5 * fsmain

    fig = Figure(figsize=(fwidth, fheight))
    canvas = FigureCanvasBase(fig)
    try:
        gs = gridspec.GridSpec(frow, fcol, fig, fleft, fbottom, 1.0 - fright, 1.0 - ftop, 0.0, 0.0)
    except TypeError:
        gs = gridspec.GridSpec(frow, fcol, fleft, fbottom, 1.0 - fright, 1.0 - ftop, 0.0, 0.0)

    # raw
    ax1 = fig.add_subplot(gs[0:f1, 0:])

    ax1.errorbar(results['input_series'][results['input_series_x']],
                 results['input_series'][results['input_series_y']],
                 results['input_series'][results['input_series_y_unc']],
                 fmt='o', color='k', ms=2)
    ax1.plot(results['input_series'][results['input_series_x']],
             results['output_series']['model'], 'r-')

    fig.text(lebels_right, fbottom + (f2 + f1 / 2) * frow_height, results['input_series_y_print'],
             fontsize=fsbig, va='center', ha='center', rotation='vertical')

    ax1.tick_params(labelbottom=False, labelsize=fsmain)

    # residuals
    ax2 = fig.add_subplot(gs[f1:f1 + f2, 0:])
    ax2.errorbar(results['input_series'][results['input_series_x']],
                 results['output_series']['residuals'],
                 results['input_series'][results['input_series_y_unc']],
                 fmt='o', color='k', ms=2)
    ax2.plot(results['input_series'][results['input_series_x']],
             results['output_series']['residuals'] * 0, 'r-')

    ax2.set_ylim(- 8 * results['statistics']['res_std'], 8 * results['statistics']['res_std'])

    ax2.set_xlabel(results['input_series_x_print'], fontsize=fsbig)
    fig.text(lebels_right, fbottom + (f2 / 2) * frow_height, 'residuals', fontsize=fsbig, va='center', ha='center',
             rotation='vertical')

    ax2.tick_params(labelsize=fsmain)

    fig.savefig(output_file, transparent=False)
    del fig


def plot_mcmc_corner(fitting_object, export_file):

    if not fitting_object.mcmc_run_complete:
        raise PyLCProcessError('MCMC not completed')

    names = []
    results = []
    print_results = []
    errors1 = []
    print_errors1 = []
    errors2 = []
    print_errors2 = []
    errors = []
    traces = []
    traces_bins = []
    traces_counts = []

    for i in fitting_object.names:
        if fitting_object.results['parameters'][i]['initial']:
            names.append(fitting_object.results['parameters'][i]['print_name'])
            results.append(fitting_object.results['parameters'][i]['value'])
            print_results.append(fitting_object.results['parameters'][i]['print_value'])
            errors1.append(fitting_object.results['parameters'][i]['m_error'])
            print_errors1.append(fitting_object.results['parameters'][i]['print_m_error'])
            errors2.append(fitting_object.results['parameters'][i]['p_error'])
            print_errors2.append(fitting_object.results['parameters'][i]['print_p_error'])
            errors.append(0.5 * (fitting_object.results['parameters'][i]['m_error'] +
                                 fitting_object.results['parameters'][i]['p_error']))
            traces.append(fitting_object.results['parameters'][i]['trace'])
            traces_bins.append(fitting_object.results['parameters'][i]['trace_bins'])
            traces_counts.append(fitting_object.results['parameters'][i]['trace_counts'])

    correlation = fitting_object.results['statistics']['corr_matrix']


    all_var = len(traces)

    fig = Figure(figsize=(2.5 * all_var + 0.5, 2.5 * all_var + 0.5))
    canvas = FigureCanvasBase(fig)
    cmap = cm.get_cmap('brg')
    try:
        gs = gridspec.GridSpec(all_var, all_var, fig, 0.5 / (2.5 * all_var + 0.5), 0.85 / (2.5 * all_var + 0.5),
                               1 - 0.5 / (2.5 * all_var + 0.5), 1 - 0.15 / (2.5 * all_var + 0.5), 0.0, 0.0)
    except TypeError:
        gs = gridspec.GridSpec(all_var, all_var, 0.5 / (2.5 * all_var + 0.5), 0.85 / (2.5 * all_var + 0.5),
                               1 - 0.5 / (2.5 * all_var + 0.5), 1 - 0.15 / (2.5 * all_var + 0.5), 0.0, 0.0)

    fontsize = 20

    for var in range(len(names)):

        ax = fig.add_subplot(gs[var, var], facecolor='w')

        ax.step(traces_bins[var], traces_counts[var], color='k', where='mid')

        ax.axvline(results[var], c='k')
        ax.axvline(results[var] - errors1[var], c='k', ls='--', lw=0.5)
        ax.axvline(results[var] + errors2[var], c='k', ls='--', lw=0.5)

        ax.set_xticks([])
        ax.set_yticks([])
        ax.tick_params(left=False, right=False, top=False, bottom=False, labelbottom=False, labelleft=False)

        ax.set_xlabel('{0}{3}{0}\n{0}{4}_{1}-{5}{2}^{1}+{6}{2}{0}'.format(
            '$', '{', '}', names[var], print_results[var], print_errors1[var], print_errors2[var]), fontsize=fontsize)

        ax.set_xlim(results[var] - 6 * errors[var], results[var] + 6 * errors[var])
        ax.set_ylim(0, ax.get_ylim()[1])

        for j in range(var + 1, all_var):

            ax2 = fig.add_subplot(gs[var, j], facecolor='w')

            binsx, binsy, final = two_d_distribution(traces[j], traces[var])
            ax2.imshow(np.where(final > 0, np.log(np.where(final > 0, final, 1)), 0),
                       extent=(np.min(binsx), np.max(binsx), np.min(binsy), np.max(binsy)),
                       cmap=cm.Greys, origin='lower', aspect='auto')

            ax2.set_xticks([])
            ax2.set_yticks([])
            ax2.tick_params(bottom=False, left=False, right=False, top=False, labelbottom=False,
                            labelleft=False, labelright=False, labeltop=False)

            ax2.set_xlim(results[j] - 6 * errors[j], results[j] + 6 * errors[j])
            ax2.set_ylim(results[var] - 6 * errors[var], results[var] + 6 * errors[var])
            text_x = ax2.get_xlim()[1] - 0.05 * (ax2.get_xlim()[1] - ax2.get_xlim()[0])
            text_y = ax2.get_ylim()[1] - 0.05 * (ax2.get_ylim()[1] - ax2.get_ylim()[0])
            ax2.text(text_x, text_y, '{0}{1}{2}'.format(r'$', str(round(correlation[var][j], 2)), '$'),
                     color=cmap(abs(correlation[var][j]) / 2.),
                     fontsize=fontsize, ha='right', va='top')

    fig.savefig(export_file, transparent=False)
    del fig


def plot_mcmc_traces(fitting_object, export_file):

    if not fitting_object.mcmc_run_complete:
        raise PyLCProcessError('MCMC not completed')

    fig = Figure(figsize=(7, 2.5 * len(fitting_object.fitted_parameters)), tight_layout=False)
    canvas = FigureCanvasBase(fig)
    cmap = cm.get_cmap('brg')

    for var_num, var in enumerate(fitting_object.fitted_parameters):

        ax = fig.add_subplot(len(fitting_object.fitted_parameters), 1, var_num + 1)
        ax.plot(fitting_object.results['parameters'][var]['trace'], 'k-', lw=0.1)
        ax.axhline(fitting_object.results['parameters'][var]['value'], c='r')
        ax.axhline(fitting_object.results['parameters'][var]['value'] - fitting_object.results['parameters'][var]['m_error'],
                   ls='--', c='r', lw=0.5)
        ax.axhline(fitting_object.results['parameters'][var]['value'] + fitting_object.results['parameters'][var]['m_error'],
                   ls='--', c='r', lw=0.5)

        ax.set_yticks([])
        ax.tick_params(left=False, right=False, labelleft=False)

        ax.set_ylabel('{0} = {1}\n{2}\n{3}'.format(r'${0}$'.format(fitting_object.results['parameters'][var]['print_name']),
                                                  r'${0}$'.format(fitting_object.results['parameters'][var]['print_value']),
                                                  r'$-{0}$'.format(fitting_object.results['parameters'][var]['print_m_error']),
                                                  r'$+{0}$'.format(fitting_object.results['parameters'][var]['print_p_error'])),
                      fontsize=15)

        if var_num != len(fitting_object.fitted_parameters) - 1:
            ax.tick_params(labelbottom='off')
        else:
            ax.set_xlabel(r'$\mathrm{iteration}$', fontsize=13)

        ax.set_ylim(fitting_object.results['parameters'][var]['value'] - 6 * fitting_object.results['parameters'][var]['m_error'],
                    fitting_object.results['parameters'][var]['value'] + 6 * fitting_object.results['parameters'][var]['m_error'])

    fig.subplots_adjust(hspace=0, wspace=0)
    fig.savefig(export_file, transparent=False)
    del fig


def plot_transit_fitting_models(results, output_file):

    for observation in results['observations']:

        data = results['observations'][observation]

        f1 = 3
        f2 = 3
        f3 = 1

        fbottom = 0.08
        ftop = 0.05
        fleft = 0.1
        fright = 0.05

        funit = 1.0
        frow = f1 + f2 + f3
        fcol = round(frow * 1.5)
        frow_height = (1.0 - fbottom - ftop) / frow

        lebels_right = fleft / 4

        fsmain = 10
        fsbig = 15
        fig = Figure(figsize=(funit * fcol / (1 - fright - fleft), funit * frow / (1 - fbottom - ftop)))
        canvas = FigureCanvasBase(fig)
        try:
            gs = gridspec.GridSpec(frow, fcol, fig, fleft, fbottom, 1.0 - fright, 1.0 - ftop, 0.0, 0.0)
        except TypeError:
            gs = gridspec.GridSpec(frow, fcol, fleft, fbottom, 1.0 - fright, 1.0 - ftop, 0.0, 0.0)

        fig.text(lebels_right, 1 - ftop / 2, 'relative flux', fontsize=fsbig, va='center', ha='left')

        # raw
        ax1 = fig.add_subplot(gs[0:f1, 0:])

        ax1.plot(data['time'], data['flux'], 'ko', ms=2)
        ax1.plot(data['time'], data['model'], 'r-', lw=1)

        fig.text(lebels_right, fbottom + (f3 + f2 + f1 / 2) * frow_height, 'raw', fontsize=fsbig, va='center',
                 ha='center', rotation='vertical')

        data_ymin = (min(data['flux']) - 3 * np.std(data['residuals']))

        data_ymax = (max(data['flux']) + 2 * np.std(data['residuals']))

        ax1.set_yticks(ax1.get_yticks()[np.where(ax1.get_yticks() > data_ymin)])

        ymin, ymax = data_ymax - 1.05 * (data_ymax - data_ymin), data_ymax

        ax1.set_ylim(ymin, ymax)

        ax1.tick_params(labelbottom=False, labelsize=fsmain)

        # de-trended
        ax2 = fig.add_subplot(gs[f1:f1+f2, 0:])

        ax2.plot(data['time'], data['detrended_flux'], 'ko', ms=2)
        ax2.plot(data['time'], data['detrended_model'], 'r-', lw=1)

        fig.text(lebels_right, fbottom + (f3 + f2 / 2) * frow_height, 'de-trended', fontsize=fsbig, va='center',
                 ha='center', rotation='vertical')

        data_ymin = (min(data['detrended_flux']) - 3 * np.std(data['detrended_residuals']))

        data_ymax = (max(data['detrended_flux']) + 2 * np.std(data['detrended_residuals']))

        ax2.set_yticks(ax2.get_yticks()[np.where(ax2.get_yticks() > data_ymin)])

        ymin, ymax = data_ymax - 1.05 * (data_ymax - data_ymin), data_ymax

        ax2.set_ylim(ymin, ymax)

        ax2.tick_params(labelbottom=False, labelsize=fsmain)

        # residuals
        ax3 = fig.add_subplot(gs[f1 + f2:f1 + f2 + f3, 0:])
        ax3.plot(data['time'], data['detrended_residuals'], 'ko', ms=2)
        ax3.plot(data['time'], np.zeros_like(data['time']), 'r-', lw=1)

        ax3.set_ylim(- 8 * np.std(data['detrended_residuals']),
                     8 * np.std(data['detrended_residuals']))

        ax3.set_xlabel('t [BJD_$\mathrm{TDB}$]', fontsize=fsbig)
        fig.text(lebels_right, fbottom + (f3 / 2) * frow_height, 'residuals', fontsize=fsbig, va='center', ha='center',
                 rotation='vertical')

        ax3.tick_params(labelsize=fsmain)

        ax1.set_title('{0} - Dataset {1} - {2}'.format(data['target'], observation + 1, data['date']))

        fig.savefig('.'.join(output_file.split('.')[:-1]) + '_dataset_{0}.'.format(observation + 1) + output_file.split('.')[-1])

