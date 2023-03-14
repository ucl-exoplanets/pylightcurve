__all__ = ['simbad_search_by_name', 'simbad_search_by_coordinates']

import time
import requests
import warnings

from astroquery.simbad import Simbad

from ..databases import plc_data
from pylightcurve.spacetime.targets import FixedTarget
from pylightcurve.spacetime.angles import Degrees, Hours, _request_angle

Simbad = Simbad()
Simbad.add_votable_fields('ids')


def _flat_name(name):

    flat_name_list = [
        [' ', ''],
        ['-', ''],
        ['cancri', 'cnc'],
        ['hatp10', 'wasp11'],
        ['wasp40', 'hatp27'],
        ['wasp51', 'hatp30'],
        ['wasp86', 'kelt12'],
        ['kelt22', 'wasp173'],
    ]

    name = name.lower()

    for char in flat_name_list:
        name = name.replace(char[0], char[1])

    return name


def fix_simbad_name(input):

    try:
        input = input.decode("utf-8")
    except:
        pass

    while '  ' in input:
        input = input.replace('  ', ' ')

    return input


def fix_simbad_coordinates(input):

    input = str(input).replace(' ', ':')

    if len(input.split(':')) == 1:
        input += ':00:00'
    elif len(input.split(':')) == 2:
        input += ':00'

    return input


def simbad_search_by_name(star_name, max_trials=10):

    try:
        star = plc_data.ecc()['stars'][plc_data.ecc()['flats'][_flat_name(star_name)]]

        return FixedTarget(
            Hours(star['ra']),
            Degrees(star['dec']),
            star['simbad_id'],
            [star['simbad_id']]
        )

    except KeyError:
        pass

    star_name = star_name.replace('HD-', 'HD ')
    star_name = star_name.replace('HR-', 'HR ')
    star_name = star_name.replace('2MASS-', '2MASS ')
    star_name = star_name.replace('GJ-', 'GJ ')
    star_name = star_name.replace('HIP-', 'HIP ')
    star_name = star_name.replace('LHS-', 'LHS ')
    star_name = star_name.replace('TYC-', 'TYC ')
    star_name = star_name.replace('KIC-', 'KIC ')
    star_name = star_name.replace('ADS-', 'ADS ')
    star_name = star_name.replace('GSC-', 'GSC ')
    star_name = star_name.replace('UCAC4-', 'UCAC4 ')

    if 'BD' in star_name:
        star_name = star_name[:3] + star_name[3:].replace('-', ' ')

    result = None
    connected = False
    trials = 0
    while not connected and trials < max_trials:
        try:
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore")
                result = Simbad.query_object(star_name)
                connected = True
                time.sleep(0.1)
        except requests.ConnectionError:
            trials += 1

    if result and len(result) > 0:

        return FixedTarget(
            Hours(fix_simbad_coordinates(result[0]['RA'])),
            Degrees(fix_simbad_coordinates(result[0]['DEC'])),
            fix_simbad_name(result[0]['MAIN_ID']),
            [fix_simbad_name(ff) for ff in result[0]['IDS'].split('|')]
        )
    else:
        print('No matches found or connection is down.')
        return None


def simbad_search_by_coordinates(ra, dec, radius=Degrees(0, 1, 0), max_trials=10):

    target = FixedTarget(ra, dec)
    _request_angle(radius)

    result = None
    connected = False
    trials = 0
    while not connected and trials < max_trials:
        try:
            with warnings.catch_warnings():
                warnings.filterwarnings("ignore")
                result = Simbad.query_region("{0}d {1}d".format(target.ra.deg(), target.dec.deg()),
                                             radius='{0}d'.format(radius.deg()))
                connected = True
                time.sleep(0.1)
        except requests.ConnectionError:
            trials += 1

    if result and len(result) > 0:
        return FixedTarget(
            Hours(fix_simbad_coordinates(result[0]['RA'])),
            Degrees(fix_simbad_coordinates(result[0]['DEC'])),
            fix_simbad_name(result[0]['MAIN_ID']),
            [fix_simbad_name(ff) for ff in result[0]['IDS'].split('|')]
        )
    else:
        print('No matches found or connection is down.')
        return None
