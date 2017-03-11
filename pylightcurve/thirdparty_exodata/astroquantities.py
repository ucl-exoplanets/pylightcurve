""" Temp module until astro units are added to quantities
"""
from quantities import *

L_s = solar_luminosity = UnitQuantity(
    'solar_luminosity',
    3.839*(10**26)*W,
    symbol='L_s',
)
L_s.latex_symbol = 'L_\odot'

R_s = solar_radius = UnitLength(
    'solar_radius',
    6.995 * (10**8) * m,
    aliases=['solar_radii'],
    symbol='R_s',
)
R_s.latex_symbol = 'R_\odot'

R_e = earth_radius = UnitLength(
    'earth_radius',
    6.371 * (10**6) * m,
    aliases=['earth_radii'],
    symbol='R_e',
)
R_e.latex_symbol = 'R_\oplus'

R_j = jupiter_radius = UnitLength(
    'jupiter_radius',
    6.9911 * (10**7) * m,
    aliases=['jupiter_radii'],
    symbol='R_j',
)
R_j.latex_symbol = 'R_J'

M_s = solar_mass = UnitMass(
    'solar_mass', 1.99*(10**30)*kg,
    aliases=['solar_masses'],
    symbol='M_s',
)
M_s.latex_symbol = 'M_\odot'

M_e = earth_mass = UnitMass(
    'earth_mass', 5.97219*(10**24)*kg,
    aliases=['earth_masses'],
    symbol='M_e',
)
M_e.latex_symbol = 'M_\oplus'

M_j = jupiter_mass = UnitMass(
    'jupiter_mass', 1.8986*(10**27)*kg,
    aliases=['jupiter_masses'],
    symbol='M_j',
)
M_j.latex_symbol = 'M_J'

Gyear = giga_year = UnitTime(
    'giga_year', 10**9*year,
    symbol='Gyr',
)

JulianDay = julian_day = JD = UnitTime(
    'julian_day', day,
    symbol='JD',
)
""" Note while quantities doesnt directly support units with an offset in most cases Julian Days are treated like days.
It is useful then to know if your working in Julian Days, MJD, BJD etc"""

ModifiedJulianDate = modified_julian_day = MJD = UnitTime(
    'modified_julian_day', day,
    symbol='MJD',
)

# Compound Units
gcm3 = CompoundUnit('g /cm**3')
gcm3.latex_symbol = 'g/cm^3'

kgm3 = CompoundUnit('kg /m**3')
kgm3.latex_symbol = 'kg/m^3'

ms2 = CompoundUnit('m/s**2')
ms2.latex_symbol = 'm/s^2'