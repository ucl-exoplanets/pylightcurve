import os
import pytest
import numpy as np
import pylightcurve as plc


def line(x_array, model_slope, model_constant, model_test):
    return model_slope * x_array + model_constant


def test_emcee():

    __location__ = os.path.abspath(os.path.dirname(__file__))

    data_x = np.arange(0, 10, 0.1)
    data_z = np.ones_like(data_x) * 0.1
    data_y = line(data_x, 1.3, 2.3, 0.0) + np.random.normal(0, data_z)

    def test_function():
        pass

    test_fit = plc.EmceeFitting(data_x, data_y, data_z, line, [1.5, 1.0, 0.0], [-5.0, -5.0, np.nan], [5.0, 5.0, np.nan],
                                200, 100000, 10000,
                                data_x_name='time', data_y_name='flux', data_x_print_name='Time',
                                data_y_print_name='Relative Flux',
                                parameters_names=['s', 'c', 't'], parameters_print_names=['slope', 'constant', 'test'],
                                function_to_call=test_function)

    test_fit.run_mcmc()
    assert round(test_fit.results['parameters']['s']['value'], 1) == 1.3
    assert round(test_fit.results['parameters']['c']['value'], 1) == 2.3

    test_fit = plc.EmceeFitting(data_x, data_y, data_z, line, [1.5, 1.0, 0.0], [-5.0, -5.0, np.nan], [5.0, 5.0, np.nan],
                                200, 100000, 10000,
                                data_x_name='time', data_y_name='flux', data_x_print_name='Time',
                                data_y_print_name='Relative Flux',
                                parameters_names=['s', 'c', 't'], parameters_print_names=['slope', 'constant', 'test'],
                                function_to_call=test_function)

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

    test_fit.run_mcmc()

    try:
        os.mkdir(os.path.join(__location__, 'test_emcee'))
    except:
        pass

    test_fit.save_all(os.path.join(__location__, 'test_emcee', 'emcee_test_data_base.pickle'))
    test_fit.save_results(os.path.join(__location__, 'test_emcee', 'emcee_test_results.txt'))
    test_fit.plot_fitting(os.path.join(__location__, 'test_emcee', 'emcee_test_fitting.pdf'))
    test_fit.plot_corner(os.path.join(__location__, 'test_emcee', 'emcee_test_correlations.pdf'))
    test_fit.plot_traces(os.path.join(__location__, 'test_emcee', 'emcee_test_traces.pdf'))


def test_values_to_print():

    assert plc.values_to_print(3.1234567, 0.12456, 0.4345) == ('3.12', '0.12', '0.43')
    assert plc.values_to_print(3.1234567, 0.4345, 0.12456) == ('3.12', '0.43', '0.12')
    assert plc.values_to_print(3.1234567, 0.456, 0.657) == ('3.12', '0.46', '0.66')
    assert plc.values_to_print(3.1234567, 0.0, 0.00) == ('3.12', '0.00', '0.00')
    assert plc.values_to_print(3.1234567, 3.3, 2.65789) == ('3.12', '3.30', '2.66')
    assert plc.values_to_print(3.1234567, 3.0, 2.0) == ('3.12', '3.00', '2.00')
