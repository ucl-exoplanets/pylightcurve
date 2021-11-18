
import pytest
import pylightcurve as plc


def test_catalogues():

    planet1 = plc.get_planet('hatp7b')
    planet2 = plc.get_planet('HAT-P-7b')
    assert planet1.name == 'HAT-P-7b'
    assert planet2.name == 'HAT-P-7b'
    assert planet1.sma_over_rs == planet2.sma_over_rs

    planet = plc.get_planet('wasp77b')
    assert planet.name == 'WASP-77Ab'

    with pytest.raises(plc.PyLCInputError):
        planet = plc.get_planet('aaaaaaaaaa')

    assert len(plc.get_all_planets()) > 0

    planet = plc.locate_planet(plc.Degrees(330.795), plc.Degrees(18.884))
    assert planet.name == 'HD209458b'

    planet = plc.locate_planet(330.795, 18.884)
    assert planet.name == 'HD209458b'

    with pytest.raises(plc.PyLCLibraryError):
        planet = plc.locate_planet(plc.Degrees(330.795), plc.Degrees(17.884))
