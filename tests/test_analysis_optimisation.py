
import os
import pytest
import numpy as np
import pylightcurve as plc


def line(x_array, model_slope, model_constant, model_test):
    return model_slope * x_array + model_constant


def test_fitting():

    __location__ = os.path.abspath(os.path.dirname(__file__))
    try:
        os.mkdir(os.path.join(__location__, 'test_emcee'))
    except:
        pass

    data_x = np.arange(0, 10, 0.1)
    data_z = np.ones_like(data_x) * 0.1
    data_y = line(data_x, 1.3, 2.3, 0.0) + np.random.normal(0, data_z)

    test_fit = plc.Fitting(data_x, data_y, data_z, line,
                           [1.5, 1.0, 0.0], [-5.0, -5.0, np.nan], [5.0, 5.0, np.nan])

    with pytest.raises(plc.PyLCProcessError):
        test_fit.save_all(os.path.join(__location__, 'emcee_test_data_base.pickle'))
    with pytest.raises(plc.PyLCProcessError):
        test_fit.save_results(os.path.join(__location__, 'emcee_test_results.txt'))
    with pytest.raises(plc.PyLCProcessError):
        test_fit.plot_corner(os.path.join(__location__, 'emcee_test_correlations.pdf'))
    with pytest.raises(plc.PyLCProcessError):
        test_fit.plot_traces(os.path.join(__location__, 'emcee_test_traces.pdf'))
    with pytest.raises(plc.PyLCProcessError):
        test_fit.plot_fitting(os.path.join(__location__, 'emcee_test_fitting.pdf'))

    assert test_fit._probability(np.array([-1.0, 0.5])) == -np.inf
    test_fit._pass()
    test_fit.run(verbose=True)
    assert round(test_fit.results['parameters']['p0']['value'], 1) in [1.2, 1.3, 1.4]
    assert round(test_fit.results['parameters']['p1']['value'], 1) in [2.2, 2.3, 2.4]
    test_fit.save_all(os.path.join(__location__, 'test_emcee', 'emcee_test_data_base.pickle'))
    test_fit.save_results(os.path.join(__location__, 'test_emcee', 'emcee_test_results.txt'))
    test_fit.plot_fitting(os.path.join(__location__, 'test_emcee', 'emcee_test_fitting.pdf'))
    test_fit.plot_corner(os.path.join(__location__, 'test_emcee', 'emcee_test_correlations.pdf'))
    test_fit.plot_traces(os.path.join(__location__, 'test_emcee', 'emcee_test_traces.pdf'))

    data_y_with_outliers = np.ones_like(data_y) * data_y
    data_y_with_outliers[-1] += 10

    test_fit = plc.Fitting(
        data_x, data_y_with_outliers, data_z, line, [1.5, 1.0, 0.0], [-5.0, -5.0, np.nan], [5.0, 5.0, np.nan],
        data_x_name='time', data_y_name='flux',
        data_x_print_name='Time', data_y_print_name='Relative Flux',
        parameters_names=['s', 'c', 't'],
        parameters_print_names=['slope', 'constant', 'test'],
        filter_outliers=True, scale_uncertainties=True,
        optimiser='curve_fit'
    )

    test_fit.run()

    plc.Fitting(data_x, data_y, data_z, line,
                [1.5, 1.0, 0.0], [-5.0, -5.0, np.nan], [5.0, 5.0, np.nan],
                walkers=10, iterations=1000, burn_in=100)

    plc.Fitting(data_x, data_y, data_z, line,
                [1.5, 1.0, 0.0], [-5.0, -5.0, np.nan], [5.0, 5.0, np.nan],
                counter='Test')

    plc.EmceeFitting(data_x, data_y, data_z, line,
                     [1.5, 1.0, 0.0], [-5.0, -5.0, np.nan], [5.0, 5.0, np.nan],
                     10, 100000, 10000)

    plc.EmceeFitting(data_x, data_y, data_z, line,
                     [1.5, 1.0, 0.0], [-5.0, -5.0, np.nan], [5.0, 5.0, np.nan],
                     walkers=10, iterations=100000, burn_in=10000)

    with pytest.raises(plc.PyLCInputError):
        plc.Fitting(data_x, data_y, data_z, line, [1.5, 1.0, 0.0], [-5.0, -5.0, np.nan], [5.0, 5.0, np.nan],
                    optimiser='xaxaxa')

    with pytest.raises(plc.PyLCInputError):
        plc.Fitting(data_x, data_y, data_z, line,
                    [1.5, 1.0, 0.0], [-5.0, -5.0, np.nan], [5.0, 5.0, np.nan],
                    counter=5)

    with pytest.raises(plc.PyLCInputError):
        plc.Fitting(data_x, data_y, data_z, line, [1.5, 1.0, 0.0], [-5.0, -5.0, np.nan], [5.0, 5.0, np.nan],
                    walkers=10, iterations=1000, burn_in=1000000)



