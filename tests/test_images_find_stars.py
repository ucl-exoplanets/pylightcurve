
import pytest
import pylightcurve as plc
import os


def test_find_stars():
    __location__ = os.path.abspath(os.path.dirname(__file__))

    fits = plc.open_fits(os.path.join(__location__, 'out_2016_08_03_19_08_10_qatar1b-001.fit'))

    single = plc.find_single_star(fits[1].data, 685, 435)
    assert len(single) == 8

    single = plc.find_single_star(fits[1].data, 0, 0)
    assert single is None
