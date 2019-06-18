from __future__ import absolute_import
from __future__ import division
from __future__ import print_function

from ._1databases import *

from .tools_counter import *
from .analysis_distribution_one_d import *


class EmceeFitting:

    def __init__(self, input_data, model, initials, limits1, limits2, walkers, iterations, burn_in,
                 names=None, print_names=None, counter=True, counter_window=False, strech_prior=1000.0):

        self.input_data = (np.array(input_data[0]), np.array(input_data[1]))

        self.model = model

        self.initials = np.array(initials)

        self.limits1 = np.array(limits1)
        for i, j in enumerate(limits1):
            if not j:
                limits1[i] = np.nan

        self.limits2 = np.array(limits2)
        for i, j in enumerate(limits2):
            if not j:
                limits2[i] = np.nan

        self.walkers = walkers

        self.iterations = (iterations / walkers) * walkers

        self.burn_in = burn_in

        if names is None:
            self.names = range(1, len(initials) + 1)
        else:
            self.names = names

        if print_names is None:
            self.print_names = map(self.names, str)
        else:
            self.print_names = print_names

        self.counter = counter

        if self.counter is not False:
            if counter is True:
                self.counter = 'MCMC'
            elif isinstance(counter, str):
                self.counter = counter

        self.counter_window = counter_window

        if self.counter_window is not False:
            if counter_window is True:
                self.counter_window = 'MCMC'
            elif isinstance(counter_window, str):
                self.counter_window = counter_window

        self.strech_prior = strech_prior

        self.results = {'input_series': {},
                        'parameters': {},
                        'parameters_final': [],
                        'output_series': {},
                        'statistics': {}}

        self.fitted_parameters = []

        self.mcmc_run_complete = False

    def run_mcmc(self):

        fitted_parameters_indices = np.where(~np.isnan(self.limits1 * self.limits2))[0]

        dimensions = len(fitted_parameters_indices)

        internal_limits1 = self.limits1[fitted_parameters_indices]
        internal_limits2 = self.limits2[fitted_parameters_indices]
        internal_initials = self.initials[fitted_parameters_indices]

        walkers_initial_positions = np.random.uniform(
            (internal_initials -
             (internal_initials - internal_limits1) / self.strech_prior)[:, None] * np.ones(self.walkers),
            (internal_initials +
             (internal_limits2 - internal_initials) / self.strech_prior)[:, None] * np.ones(self.walkers))
        walkers_initial_positions = np.swapaxes(walkers_initial_positions, 0, 1)

        def internal_model(theta):
            parameters = self.initials
            parameters[fitted_parameters_indices] = theta
            return self.model(*parameters)

        def likelihood(theta, data_y, data_y_error):
            if ((internal_limits1 < theta) & (theta < internal_limits2)).all():
                chi = (data_y - internal_model(theta)) / data_y_error
                return -0.5 * (np.sum(chi * chi) + np.sum(np.log(2.0 * np.pi * (data_y_error * data_y_error))))
            else:
                return -np.inf

        def prior(theta):
            if ((internal_limits1 < theta) & (theta < internal_limits2)).all():
                return 0.0
            return -np.inf

        if self.counter or self.counter_window:
            counter = Counter(self.counter, self.counter_window, self.iterations + self.walkers, show_every=100)

            def probability(theta, data_y, data_y_error):
                counter.update()
                return prior(theta) + likelihood(theta, data_y, data_y_error)
        else:
            def probability(theta, data_y, data_y_error):
                return prior(theta) + likelihood(theta, data_y, data_y_error)

        sampler = emcee.EnsembleSampler(self.walkers, dimensions, probability, args=self.input_data)
        sampler.run_mcmc(walkers_initial_positions, int(self.iterations) // int(self.walkers))
        mcmc_results = sampler.flatchain

        self.results['input_series']['value'] = self.input_data[0]
        self.results['input_series']['error'] = self.input_data[1]

        for var in range(len(self.names)):

            if np.isnan(self.limits1[var]):

                variable = {'name': self.names[var], 'print_name': self.print_names[var],
                            'initial': None, 'min_allowed': None, 'max_allowed': None,
                            'trace': None, 'trace_bins': None, 'trace_counts': None,
                            'value': self.initials[var], 'm_error': None, 'p_error': None,
                            'print_value': self.initials[var], 'print_m_error': '-', 'print_p_error': '-'}

            else:
                trace = mcmc_results[:, np.where(fitted_parameters_indices == var)[0][0]]
                trace = trace.reshape(int(self.walkers), int(self.iterations) // int(self.walkers))
                trace = (np.swapaxes(trace, 0, 1).flatten())[self.burn_in:]

                bins, counts, value, m_error, p_error, print_value, print_m_error, print_p_error = \
                    one_d_distribution(trace, confidence_interval=0.68)

                variable = {'name': self.names[var], 'print_name': self.print_names[var],
                            'initial': self.initials[var],
                            'min_allowed': self.limits1[var], 'max_allowed': self.limits2[var],
                            'trace': trace, 'trace_bins': bins, 'trace_counts': counts,
                            'value': value, 'm_error': m_error, 'p_error': p_error,
                            'print_value': print_value, 'print_m_error': print_m_error, 'print_p_error': print_p_error}

                self.fitted_parameters.append(self.names[var])

            self.results['parameters'][self.names[var]] = variable
            self.results['parameters_final'].append(variable['value'])

        self.results['output_series']['model'] = self.model(*self.results['parameters_final'])
        self.results['output_series']['residuals'] = self.input_data[0] - self.results['output_series']['model']

        def correlation(x, y):
            n = len(x)
            mx = np.mean(x)
            sx = np.std(x)
            my = np.mean(y)
            sy = np.std(y)
            return np.round(np.sum((x - mx) * (y - my)) / ((n - 1) * sx * sy), 2)

        correlation_matrix = np.zeros((len(self.fitted_parameters), len(self.fitted_parameters)))

        for ii in range(len(self.fitted_parameters)):
            for jj in range(ii, len(self.fitted_parameters)):
                correlation_matrix[ii][jj] = \
                    correlation(self.results['parameters'][self.fitted_parameters[ii]]['trace'],
                                self.results['parameters'][self.fitted_parameters[jj]]['trace'])

        res_autocorr = np.correlate(self.results['output_series']['residuals'],
                                    self.results['output_series']['residuals'], mode='full')
        res_autocorr = res_autocorr[res_autocorr.size // 2:] / res_autocorr[res_autocorr.size // 2:][0]

        self.results['statistics']['res_autocorr'] = res_autocorr
        self.results['statistics']['res_std'] = np.std(self.results['output_series']['residuals'])
        self.results['statistics']['corr_matrix'] = correlation_matrix
        self.results['statistics']['corr_variables'] = ','.join(self.fitted_parameters)

        self.mcmc_run_complete = True

    def save_all(self, export_file):

        if not self.mcmc_run_complete:
            raise RuntimeError('MCMC not completed')

        pickle.dump(self.results, open(export_file, 'wb'))

    def save_results(self, export_file):

        if not self.mcmc_run_complete:
            raise RuntimeError('MCMC not completed')

        w = open(export_file, 'w')

        w.write('# variable\tresult\tuncertainty\n')

        for i in self.names:
            w.write('{0}\t{1}\t-{2} +{3}\n'.format(self.results['parameters'][i]['name'],
                                                   self.results['parameters'][i]['print_value'],
                                                   self.results['parameters'][i]['print_m_error'],
                                                   self.results['parameters'][i]['print_p_error']))

        w.close()

    def plot_corner(self, export_file):

        def correlation(x, y):
            n = len(x)
            mx = np.mean(x)
            sx = np.std(x)
            my = np.mean(y)
            sy = np.std(y)
            return np.round(np.sum((x - mx) * (y - my)) / ((n - 1) * sx * sy), 2)

        def td_distribution(datax, datay):

            datax = np.array(datax)
            median = np.median(datax)
            med = np.sqrt(np.median((datax - median) ** 2))
            xstep = med / 5.0
            xmin = min(datax)
            xmax = max(datax)
            x_size = int(round((xmax - xmin) / xstep)) + 1
            datax = np.int_((datax - xmin) / xstep)
            datay = np.array(datay)
            median = np.median(datay)
            med = np.sqrt(np.median((datay - median) ** 2))
            ystep = med / 5.0
            ymin = min(datay)
            ymax = max(datay)
            y_size = int(round((ymax - ymin) / ystep)) + 1
            datay = np.int_((datay - ymin) / ystep)

            yx_size = x_size * y_size
            yx = datay * x_size + datax

            yx = np.bincount(yx)
            yx = np.insert(yx, len(yx), np.zeros(yx_size - len(yx)))

            xx, yy = np.meshgrid(xmin + np.arange(x_size) * xstep, ymin + np.arange(y_size) * ystep)

            final = np.reshape(yx, (y_size, x_size))
            plt.imshow(np.where(final > 0, np.log(np.where(final > 0, final, 1)), 0),
                       extent=(np.min(xx), np.max(xx), np.min(yy), np.max(yy)),
                       cmap=plt.cm.Greys, origin='lower', aspect='auto')

        if not self.mcmc_run_complete:
            raise RuntimeError('MCMC not completed')

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

        for i in self.names:
            if self.results['parameters'][i]['initial']:
                names.append(self.results['parameters'][i]['print_name'])
                results.append(self.results['parameters'][i]['value'])
                print_results.append(self.results['parameters'][i]['print_value'])
                errors1.append(self.results['parameters'][i]['m_error'])
                print_errors1.append(self.results['parameters'][i]['print_m_error'])
                errors2.append(self.results['parameters'][i]['p_error'])
                print_errors2.append(self.results['parameters'][i]['print_p_error'])
                errors.append(0.5 * (self.results['parameters'][i]['m_error'] +
                                     self.results['parameters'][i]['p_error']))
                traces.append(self.results['parameters'][i]['trace'])
                traces_bins.append(self.results['parameters'][i]['trace_bins'])
                traces_counts.append(self.results['parameters'][i]['trace_counts'])

        all_var = len(traces)
        fig = plt.figure(figsize=(2.5 * all_var, 2.5 * all_var))
        fig.set_tight_layout(False)
        cmap = matplotlib.cm.get_cmap('brg')

        for var in range(len(names)):

            try:
                plt.subplot(all_var, all_var, all_var * var + var + 1, facecolor='w')
            except AttributeError:
                plt.subplot(all_var, all_var, all_var * var + var + 1, axisbg='w')

            plt.step(traces_bins[var], traces_counts[var], color='k', where='mid')

            plt.axvline(results[var], c='k')
            plt.axvline(results[var] - errors1[var], c='k', ls='--', lw=0.5)
            plt.axvline(results[var] + errors2[var], c='k', ls='--', lw=0.5)

            plt.xticks(plt.xticks()[0], np.ones_like(plt.yticks()[0]))
            plt.yticks(plt.yticks()[0], np.ones_like(plt.yticks()[0]))
            plt.tick_params(left='off', right='off', top='off', bottom='off', labelbottom='off', labelleft='off')

            plt.xlabel('{0}\n{1}\n{2}\n{3}'.format(r'${0}$'.format(names[var]), r'${0}$'.format(print_results[var]),
                                                   r'$-{0}$'.format(print_errors1[var]),
                                                   r'$+{0}$'.format(print_errors2[var])), fontsize=20)

            plt.xlim(results[var] - 6 * errors[var], results[var] + 6 * errors[var])
            plt.ylim(0, plt.ylim()[1])

            for j in range(var + 1, all_var):

                plt.subplot(all_var, all_var, all_var * var + 1 + j)
                td_distribution(traces[j], traces[var])

                plt.yticks(plt.yticks()[0], np.arange(len(plt.yticks()[0])))
                plt.xticks(plt.xticks()[0], np.arange(len(plt.xticks()[0])))
                plt.tick_params(bottom='off', left='off', right='off', top='off', labelbottom='off',
                                labelleft='off', labelright='off', labeltop='off')

                plt.xlim(results[j] - 6 * errors[j], results[j] + 6 * errors[j])
                plt.ylim(results[var] - 6 * errors[var], results[var] + 6 * errors[var])
                text_x = plt.xlim()[1] - 0.05 * (plt.xlim()[1] - plt.xlim()[0])
                text_y = plt.ylim()[1] - 0.05 * (plt.ylim()[1] - plt.ylim()[0])
                plt.text(text_x, text_y, '{0}{1}{2}'.format(r'$', str(correlation(traces[j], traces[var])), '$'),
                         color=cmap(abs(correlation(traces[j], traces[var])) / 2.),
                         fontsize=20, ha='right', va='top')

        plt.subplots_adjust(hspace=0, wspace=0)
        plt.savefig(export_file, bbox_inches='tight', transparent=False)
        plt.close('all')

    def plot_traces(self, export_file):

        if not self.mcmc_run_complete:
            raise RuntimeError('MCMC not completed')

        plt.figure(figsize=(7, 2.5 * len(self.fitted_parameters)))

        for var_num, var in enumerate(self.fitted_parameters):

            plt.subplot(len(self.fitted_parameters), 1, var_num + 1)
            plt.plot(self.results['parameters'][var]['trace'], 'k-', lw=0.1)
            plt.axhline(self.results['parameters'][var]['value'], c='r')
            plt.axhline(self.results['parameters'][var]['value'] - self.results['parameters'][var]['m_error'],
                        ls='--', c='r', lw=0.5)
            plt.axhline(self.results['parameters'][var]['value'] + self.results['parameters'][var]['m_error'],
                        ls='--', c='r', lw=0.5)

            plt.yticks(plt.yticks()[0], np.ones_like(plt.yticks()[0]))
            plt.tick_params(left='off', right='off', labelleft='off')

            plt.ylabel('{0}\n{1}\n{2}\n{3}'.format(r'${0}$'.format(self.results['parameters'][var]['print_name']),
                                                   r'${0}$'.format(self.results['parameters'][var]['print_value']),
                                                   r'$-{0}$'.format(self.results['parameters'][var]['print_m_error']),
                                                   r'$+{0}$'.format(self.results['parameters'][var]['print_p_error'])),
                       fontsize=15)

            if var_num != len(self.fitted_parameters) - 1:
                plt.tick_params(labelbottom='off')
            else:
                plt.xlabel(r'$\mathrm{iteration}$', fontsize=15)

            plt.ylim(self.results['parameters'][var]['value'] - 6 * self.results['parameters'][var]['m_error'],
                     self.results['parameters'][var]['value'] + 6 * self.results['parameters'][var]['m_error'])

        plt.subplots_adjust(hspace=0, wspace=0)
        plt.savefig(export_file, bbox_inches='tight', transparent=True)
        plt.close('all')
