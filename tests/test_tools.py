
import pytest
import numpy as np
import pylightcurve as plc
import os
import glob
from astropy.io import fits as pf


def test_files():

    __location__ = os.path.abspath(os.path.dirname(__file__))

    xx = plc.open_dict(os.path.join(__location__, 'test1.pickle'))

    xx = {'a': 1}
    plc.save_dict(xx, os.path.join(__location__, 'test.pickle'))
    xx = plc.open_dict(os.path.join(__location__, 'test.pickle'))
    os.remove(os.path.join(__location__, 'test.pickle'))
    yy = plc.copy_dict(xx)
    assert id(yy) != id(xx)
    del xx, yy

    xx = {'a': 1}
    plc.save_yaml(xx, os.path.join(__location__, 'test.yaml'))
    xx = plc.open_yaml(os.path.join(__location__, 'test.yaml'))
    os.remove(os.path.join(__location__, 'test.yaml'))
    yy = plc.copy_yaml(xx)
    assert id(yy) != id(xx)
    del xx, yy

    for file in glob.glob(os.path.join(__location__, 'test?.fits')):

        xx = plc.open_fits(file)
        yy = pf.open(file)

        if os.path.isfile(os.path.join(__location__, 'test.fits')):
            os.remove(os.path.join(__location__, 'test.fits'))

        plc.save_fits(xx, os.path.join(__location__, 'test.fits'))

        os.remove(os.path.join(__location__, 'test.fits'))

        plc.save_fits(yy, os.path.join(__location__, 'test.fits'))

        os.remove(os.path.join(__location__, 'test.fits'))

    xx = pf.HDUList([pf.PrimaryHDU(), pf.ImageHDU(data=np.array([[1, 2], [1, 2]]), name='SCI'),
                     pf.ImageHDU(data=np.array([[1, 2], [1, 2]]), name='ERR')])

    assert plc.fits_sci(xx)[0] == 1
    assert plc.fits_err(xx)[0] == 2
    assert plc.fits_sci_err(xx)[0][0] == 1
    assert plc.fits_sci_err(xx)[0][1] == 2
    #
    # xx[1].header.set('A?χ', 3)
    # plc.save_fits(xx, os.path.join(__location__, 'test.fits'))
    # xx = plc.open_fits(os.path.join(__location__, 'test.fits'))
    # os.remove(os.path.join(__location__, 'test.fits'))
    # yy = plc.copy_fits(xx)
    # assert id(yy) != id(xx)
    # del xx, yy

    assert plc.open_dict_online('a') is False

    plc.download('https://raw.githubusercontent.com/ucl-exoplanets/pylightcurve/master/README.md',
                 os.path.join(__location__, 'test.txt'))
    os.remove(os.path.join(__location__, 'test.txt'))

    plc.download('https://raw.githubusercontent.com/ucl-exoplanets/pylightcurve/master/README.mx',
                 os.path.join(__location__, 'test.txt'))

    plc.zip_directory(os.path.join(__location__, 'test_emcee'), os.path.join(__location__, 'test_emcee.zip'))


def test_counter():

    xx = 1000

    counter1 = plc.Counter('test1', xx)
    counter2 = plc.Counter('test2', xx, ignore=10)
    counter3 = plc.Counter('test3', xx, show_every=10)
    counter4 = plc.Counter('test4', xx, increment=10)

    for i in range(xx):
        counter1.update()
        counter2.update()
        counter3.update()
        counter4.update(message='test')


def test_find_stars():
    __location__ = os.path.abspath(os.path.dirname(__file__))

    fits = plc.open_fits(os.path.join(__location__, 'test1.fits'))

    single = plc.find_single_star(fits[1].data, 685, 435)
    assert len(single) == 8

    single = plc.find_single_star(fits[1].data, 685, 435, maxfev=2)
    assert single is None