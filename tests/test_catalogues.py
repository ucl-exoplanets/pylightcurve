
import pytest
import pylightcurve as plc
from pylightcurve.catalogues.simbad import fix_simbad_coordinates


def test_ecc():

    planet1 = plc.get_planet('hatp7b')
    planet2 = plc.get_planet('HAT-P-7b')
    assert planet1.name == 'HAT-P-7b'
    assert planet2.name == 'HAT-P-7b'
    assert planet1.sma_over_rs == planet2.sma_over_rs

    planet = plc.get_planet('wasp77b')
    assert planet.name == 'WASP-77Ab'

    with pytest.raises(plc.PyLCInputError):
        plc.get_planet('aaaaaaaaaa')

    assert len(plc.get_all_planets()) > 0

    planet = plc.locate_planet(plc.Degrees(330.795), plc.Degrees(18.884))
    assert planet.name == 'HD209458b'

    planet = plc.locate_planet(330.795, 18.884, 0.1)
    assert planet.name == 'HD209458b'

    with pytest.raises(plc.PyLCLibraryError):
        plc.locate_planet(plc.Degrees(330.795), plc.Degrees(17.884))

    assert plc.get_system('HD189987976') == []
    assert plc.get_system('M42') == []
    assert len(plc.get_system('HD209458b')) == 1

    assert len(plc.locate_system(330.795, 18.884, 0.1)) == 1


def test_simbad():

    target = plc.simbad_search_by_name('TYC 1949-1705-1')

    assert target.name == 'TYC 1949-1705-1'

    target = plc.simbad_search_by_name('XO-1')

    assert target.name == 'BD+28 2507'

    target = plc.simbad_search_by_coordinates(target.ra, target.dec)

    assert target.name == 'BD+28 2507b'

    assert fix_simbad_coordinates('16') == '16:00:00'
    assert fix_simbad_coordinates('16:00') == '16:00:00'

    assert plc.simbad_search_by_name('BDXXXXXXXXXXXXXXXXX') is None
    assert plc.simbad_search_by_coordinates(plc.Hours(1), plc.Degrees(0), radius=plc.Degrees(0.001)) is None
