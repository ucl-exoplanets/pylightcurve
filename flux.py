__all__ = ['flux_drop']

import numpy as np

pi = np.pi

# coefficients from https://pomax.github.io/bezierinfo/legendre-gauss.html

gauss30 = [
    [0.1028526528935588, -0.0514718425553177],
    [0.1028526528935588, 0.0514718425553177],
    [0.1017623897484055, -0.1538699136085835],
    [0.1017623897484055, 0.1538699136085835],
    [0.0995934205867953, -0.2546369261678899],
    [0.0995934205867953, 0.2546369261678899],
    [0.0963687371746443, -0.3527047255308781],
    [0.0963687371746443, 0.3527047255308781],
    [0.0921225222377861, -0.4470337695380892],
    [0.0921225222377861, 0.4470337695380892],
    [0.0868997872010830, -0.5366241481420199],
    [0.0868997872010830, 0.5366241481420199],
    [0.0807558952294202, -0.6205261829892429],
    [0.0807558952294202, 0.6205261829892429],
    [0.0737559747377052, -0.6978504947933158],
    [0.0737559747377052, 0.6978504947933158],
    [0.0659742298821805, -0.7677774321048262],
    [0.0659742298821805, 0.7677774321048262],
    [0.0574931562176191, -0.8295657623827684],
    [0.0574931562176191, 0.8295657623827684],
    [0.0484026728305941, -0.8825605357920527],
    [0.0484026728305941, 0.8825605357920527],
    [0.0387991925696271, -0.9262000474292743],
    [0.0387991925696271, 0.9262000474292743],
    [0.0287847078833234, -0.9600218649683075],
    [0.0287847078833234, 0.9600218649683075],
    [0.0184664683110910, -0.9836681232797472],
    [0.0184664683110910, 0.9836681232797472],
    [0.0079681924961666, -0.9968934840746495],
    [0.0079681924961666, 0.9968934840746495]
]

gauss60 = [
    [0.0519078776312206, -0.0259597723012478],
    [0.0519078776312206, 0.0259597723012478],
    [0.0517679431749102, -0.0778093339495366],
    [0.0517679431749102, 0.0778093339495366],
    [0.0514884515009809, -0.1294491353969450],
    [0.0514884515009809, 0.1294491353969450],
    [0.0510701560698556, -0.1807399648734254],
    [0.0510701560698556, 0.1807399648734254],
    [0.0505141845325094, -0.2315435513760293],
    [0.0505141845325094, 0.2315435513760293],
    [0.0498220356905502, -0.2817229374232617],
    [0.0498220356905502, 0.2817229374232617],
    [0.0489955754557568, -0.3311428482684482],
    [0.0489955754557568, 0.3311428482684482],
    [0.0480370318199712, -0.3796700565767980],
    [0.0480370318199712, 0.3796700565767980],
    [0.0469489888489122, -0.4271737415830784],
    [0.0469489888489122, 0.4271737415830784],
    [0.0457343797161145, -0.4735258417617071],
    [0.0457343797161145, 0.4735258417617071],
    [0.0443964787957871, -0.5186014000585697],
    [0.0443964787957871, 0.5186014000585697],
    [0.0429388928359356, -0.5622789007539445],
    [0.0429388928359356, 0.5622789007539445],
    [0.0413655512355848, -0.6044405970485104],
    [0.0413655512355848, 0.6044405970485104],
    [0.0396806954523808, -0.6449728284894770],
    [0.0396806954523808, 0.6449728284894770],
    [0.0378888675692434, -0.6837663273813555],
    [0.0378888675692434, 0.6837663273813555],
    [0.0359948980510845, -0.7207165133557304],
    [0.0359948980510845, 0.7207165133557304],
    [0.0340038927249464, -0.7557237753065856],
    [0.0340038927249464, 0.7557237753065856],
    [0.0319212190192963, -0.7886937399322641],
    [0.0319212190192963, 0.7886937399322641],
    [0.0297524915007889, -0.8195375261621458],
    [0.0297524915007889, 0.8195375261621458],
    [0.0275035567499248, -0.8481719847859296],
    [0.0275035567499248, 0.8481719847859296],
    [0.0251804776215212, -0.8745199226468983],
    [0.0251804776215212, 0.8745199226468983],
    [0.0227895169439978, -0.8985103108100460],
    [0.0227895169439978, 0.8985103108100460],
    [0.0203371207294573, -0.9200784761776275],
    [0.0203371207294573, 0.9200784761776275],
    [0.0178299010142077, -0.9391662761164232],
    [0.0178299010142077, 0.9391662761164232],
    [0.0152746185967848, -0.9557222558399961],
    [0.0152746185967848, 0.9557222558399961],
    [0.0126781664768160, -0.9697017887650528],
    [0.0126781664768160, 0.9697017887650528],
    [0.0100475571822880, -0.9810672017525982],
    [0.0100475571822880, 0.9810672017525982],
    [0.0073899311633455, -0.9897878952222218],
    [0.0073899311633455, 0.9897878952222218],
    [0.0047127299269536, -0.9958405251188381],
    [0.0047127299269536, 0.9958405251188381],
    [0.0020268119688738, -0.9992101232274361],
    [0.0020268119688738, 0.9992101232274361],
]

gausstab = [np.swapaxes(gauss30, 0, 1), np.swapaxes(gauss60, 0, 1)]


# integral definitions for claret method


def integral_r_claret(limb_darkening_coefficients, r):
    a1, a2, a3, a4 = limb_darkening_coefficients
    mu44 = 1.0 - r * r
    mu24 = np.sqrt(mu44)
    mu14 = np.sqrt(mu24)
    return - (2.0 * (1.0 - a1 - a2 - a3 - a4) / 4) * mu44 \
           - (2.0 * a1 / 5) * mu44 * mu14 \
           - (2.0 * a2 / 6) * mu44 * mu24 \
           - (2.0 * a3 / 7) * mu44 * mu24 * mu14 \
           - (2.0 * a4 / 8) * mu44 * mu44


def num_claret(r, limb_darkening_coefficients, rprs, z):
    a1, a2, a3, a4 = limb_darkening_coefficients
    rsq = r * r
    mu44 = 1.0 - rsq
    mu24 = np.sqrt(mu44)
    mu14 = np.sqrt(mu24)
    return ((1.0 - a1 - a2 - a3 - a4) + a1 * mu14 + a2 * mu24 + a3 * mu24 * mu14 + a4 * mu44) \
        * r * np.arccos(np.minimum((-rprs ** 2 + z * z + rsq) / (2.0 * z * r), 1.0))


def integral_r_f_claret(limb_darkening_coefficients, rprs, z, r1, r2, precision=0):
    x1 = (r2 - r1) / 2
    x2 = (r2 + r1) / 2
    return x1 * np.sum(gausstab[precision][0][:, None] *
                       num_claret(x1[None, :] * gausstab[precision][1][:, None] + x2[None, :],
                                  limb_darkening_coefficients, rprs, z), 0)


# integral definitions for linear method


def integral_r_linear(limb_darkening_coefficients, r):
    a1 = limb_darkening_coefficients
    musq = 1 - r * r
    return (-1.0 / 6) * musq * (3.0 + a1 * (-3.0 + 2.0 * np.sqrt(musq)))


def num_linear(r, limb_darkening_coefficients, rprs, z):
    a1 = limb_darkening_coefficients
    rsq = r * r
    return (1.0 - a1 * (1.0 - np.sqrt(1.0 - rsq))) \
        * r * np.arccos(np.minimum((-rprs ** 2 + z * z + rsq) / (2.0 * z * r), 1.0))


def integral_r_f_linear(limb_darkening_coefficients, rprs, z, r1, r2, precision=0):
    x1 = (r2 - r1) / 2
    x2 = (r2 + r1) / 2
    return x1 * np.sum(gausstab[precision][0][:, None] *
                       num_linear(x1[None, :] * gausstab[precision][1][:, None] + x2[None, :],
                                  limb_darkening_coefficients, rprs, z), 0)


# integral definitions for quadratic method


def integral_r_quad(limb_darkening_coefficients, r):
    a1, a2 = limb_darkening_coefficients
    musq = 1 - r * r
    mu = np.sqrt(musq)
    return (1.0 / 12) * (-4.0 * (a1 + 2.0 * a2) * mu * musq + 6.0 * (-1 + a1 + a2) * musq + 3.0 * a2 * musq * musq)


def num_quad(r, limb_darkening_coefficients, rprs, z):
    a1, a2 = limb_darkening_coefficients
    rsq = r * r
    cc = 1.0 - np.sqrt(1.0 - rsq)
    return (1.0 - a1 * cc - a2 * cc * cc) \
        * r * np.arccos(np.minimum((-rprs ** 2 + z * z + rsq) / (2.0 * z * r), 1.0))


def integral_r_f_quad(limb_darkening_coefficients, rprs, z, r1, r2, precision=0):
    x1 = (r2 - r1) / 2
    x2 = (r2 + r1) / 2
    return x1 * np.sum(gausstab[precision][0][:, None] *
                       num_quad(x1[None, :] * gausstab[precision][1][:, None] + x2[None, :],
                                limb_darkening_coefficients, rprs, z), 0)


# integral definitions for square root method


def integral_r_sqrt(limb_darkening_coefficients, r):
    a1, a2 = limb_darkening_coefficients
    musq = 1 - r * r
    mu = np.sqrt(musq)
    return ((-2.0 / 5) * a2 * np.sqrt(mu) - (1.0 / 3) * a1 * mu + (1.0 / 2) * (-1 + a1 + a2)) * musq


def num_sqrt(r, limb_darkening_coefficients, rprs, z):
    a1, a2 = limb_darkening_coefficients
    rsq = r * r
    mu = np.sqrt(1.0 - rsq)
    return (1.0 - a1 * (1 - mu) - a2 * (1.0 - np.sqrt(mu))) \
        * r * np.arccos(np.minimum((-rprs ** 2 + z * z + rsq) / (2.0 * z * r), 1.0))


def integral_r_f_sqrt(limb_darkening_coefficients, rprs, z, r1, r2, precision=0):
    x1 = (r2 - r1) / 2
    x2 = (r2 + r1) / 2
    return x1 * np.sum(gausstab[precision][0][:, None] *
                       num_sqrt(x1[None, :] * gausstab[precision][1][:, None] + x2[None, :],
                                limb_darkening_coefficients, rprs, z), 0)


# dictionaries containing the different methods,
# if you define a new method, include the functions in the dictionary as well

integral_r = {
    'claret': integral_r_claret,
    'linear': integral_r_linear,
    'quad': integral_r_quad,
    'sqrt': integral_r_sqrt
}

integral_r_f = {
    'claret': integral_r_f_claret,
    'linear': integral_r_f_linear,
    'quad': integral_r_f_quad,
    'sqrt': integral_r_f_sqrt,
}


def integral_centred(limb_darkening_coefficients, rprs, ww1, ww2, method='claret'):
    return (integral_r[method](limb_darkening_coefficients, rprs)
            - integral_r[method](limb_darkening_coefficients, 0.0)) * np.abs(ww2 - ww1)


def integral_plus_core(limb_darkening_coefficients, rprs, z, ww1, ww2, method='claret', precision=0):
    if len(z) == 0:
        return z
    rr1 = z * np.cos(ww1) + np.sqrt(np.maximum(rprs ** 2 - (z * np.sin(ww1)) ** 2, 0))
    rr1 = np.clip(rr1, 0, 1)
    rr2 = z * np.cos(ww2) + np.sqrt(np.maximum(rprs ** 2 - (z * np.sin(ww2)) ** 2, 0))
    rr2 = np.clip(rr2, 0, 1)
    w1 = np.minimum(ww1, ww2)
    r1 = np.minimum(rr1, rr2)
    w2 = np.maximum(ww1, ww2)
    r2 = np.maximum(rr1, rr2)
    parta = integral_r[method](limb_darkening_coefficients, 0.0) * (w1 - w2)
    partb = integral_r[method](limb_darkening_coefficients, r1) * w2
    partc = integral_r[method](limb_darkening_coefficients, r2) * (-w1)
    partd = integral_r_f[method](limb_darkening_coefficients, rprs, z, r1, r2, precision=precision)
    return parta + partb + partc + partd


# def integral_plus(limb_darkening_coefficients, rprs, z, ww1, ww2):
#     split = np.where(ww1 * ww2 < 0)
#     if len(split[0]) == 0:
#         return integral_plus_core(limb_darkening_coefficients, rprs, z, np.abs(ww1), np.abs(ww2))
#     else:
#         w1 = np.minimum(ww1, ww2)
#         w2 = np.maximum(ww1, ww2)
#         intplus = integral_plus_core(
#             limb_darkening_coefficients, rprs, z, np.where(w1 * w2 < 0, 0, np.abs(w1)), np.abs(w2))
#         intplus[split] = \
#             intplus[split] + integral_plus_core(limb_darkening_coefficients, rprs, z[split], 0, np.abs(w1[split]))
#         return intplus


def integral_minus_core(limb_darkening_coefficients, rprs, z, ww1, ww2, method='claret', precision=0):
    if len(z) == 0:
        return z
    rr1 = z * np.cos(ww1) - np.sqrt(np.maximum(rprs ** 2 - (z * np.sin(ww1)) ** 2, 0))
    rr1 = np.clip(rr1, 0, 1)
    rr2 = z * np.cos(ww2) - np.sqrt(np.maximum(rprs ** 2 - (z * np.sin(ww2)) ** 2, 0))
    rr2 = np.clip(rr2, 0, 1)
    w1 = np.minimum(ww1, ww2)
    r1 = np.minimum(rr1, rr2)
    w2 = np.maximum(ww1, ww2)
    r2 = np.maximum(rr1, rr2)
    parta = integral_r[method](limb_darkening_coefficients, 0.0) * (w1 - w2)
    partb = integral_r[method](limb_darkening_coefficients, r1) * (-w1)
    partc = integral_r[method](limb_darkening_coefficients, r2) * w2
    partd = integral_r_f[method](limb_darkening_coefficients, rprs, z, r1, r2, precision=precision)
    return parta + partb + partc - partd


# def integral_minus(limb_darkening_coefficients, rprs, z, ww1, ww2):
#     split = np.where(ww1 * ww2 < 0)
#     if len(split[0]) == 0:
#         return integral_minus_core(limb_darkening_coefficients, rprs, z, np.abs(ww1), np.abs(ww2))
#     else:
#         w1 = np.minimum(ww1, ww2)
#         w2 = np.maximum(ww1, ww2)
#         intmins = integral_minus_core(
#             limb_darkening_coefficients, rprs, z, np.where(w1 * w2 < 0, 0, np.abs(w1)), np.abs(w2))
#         intmins[split] = \
#             intmins[split] + integral_minus_core(limb_darkening_coefficients, rprs, z[split], 0, np.abs(w1[split]))
#         return intmins


def flux_drop(limb_darkening_coefficients, rp_over_rs, z_over_rs, method='claret', precision=0):

    if len(z_over_rs) == 0:
        return np.array([])

    # cases
    zsq = z_over_rs * z_over_rs
    sum_z_rprs = z_over_rs + rp_over_rs
    dif_z_rprs = rp_over_rs - z_over_rs
    sqr_dif_z_rprs = zsq - rp_over_rs ** 2
    case0 = np.where((z_over_rs == 0) & (rp_over_rs <= 1))
    case1 = np.where((z_over_rs < rp_over_rs) & (sum_z_rprs <= 1))
    casea = np.where((z_over_rs < rp_over_rs) & (sum_z_rprs > 1) & (dif_z_rprs < 1))
    caseb = np.where((z_over_rs < rp_over_rs) & (sum_z_rprs > 1) & (dif_z_rprs > 1))
    case2 = np.where((z_over_rs == rp_over_rs) & (sum_z_rprs <= 1))
    casec = np.where((z_over_rs == rp_over_rs) & (sum_z_rprs > 1))
    case3 = np.where((z_over_rs > rp_over_rs) & (sum_z_rprs < 1))
    case4 = np.where((z_over_rs > rp_over_rs) & (sum_z_rprs == 1))
    case5 = np.where((z_over_rs > rp_over_rs) & (sum_z_rprs > 1) & (sqr_dif_z_rprs < 1))
    case6 = np.where((z_over_rs > rp_over_rs) & (sum_z_rprs > 1) & (sqr_dif_z_rprs == 1))
    case7 = np.where((z_over_rs > rp_over_rs) & (sum_z_rprs > 1) & (sqr_dif_z_rprs > 1) & (-1 < dif_z_rprs))
    plus_case = np.concatenate((case1[0], case2[0], case3[0], case4[0], case5[0], casea[0], casec[0]))
    minus_case = np.concatenate((case3[0], case4[0], case5[0], case6[0], case7[0]))
    star_case = np.concatenate((case5[0], case6[0], case7[0], casea[0], casec[0]))

    # cross points
    ph = np.arccos(np.clip((1.0 - rp_over_rs ** 2 + zsq) / (2.0 * z_over_rs), -1, 1))
    theta_1 = np.zeros(len(z_over_rs))
    ph_case = np.concatenate((case5[0], casea[0], casec[0]))
    theta_1[ph_case] = ph[ph_case]
    theta_2 = np.arcsin(np.minimum(rp_over_rs / z_over_rs, 1))
    theta_2[case1] = pi
    theta_2[case2] = pi / 2.0
    theta_2[casea] = pi
    theta_2[casec] = pi / 2.0
    theta_2[case7] = ph[case7]

    # flux_upper
    plusflux = np.zeros(len(z_over_rs))
    plusflux[plus_case] = integral_plus_core(limb_darkening_coefficients, rp_over_rs,
                                             z_over_rs[plus_case], theta_1[plus_case], theta_2[plus_case],
                                             method=method, precision=precision)
    if len(case0[0]) > 0:
        plusflux[case0] = integral_centred(limb_darkening_coefficients, rp_over_rs, 0.0, pi, method=method)
    if len(caseb[0]) > 0:
        plusflux[caseb] = integral_centred(limb_darkening_coefficients, 1, 0.0, pi, method=method)

    # flux_lower
    minsflux = np.zeros(len(z_over_rs))
    minsflux[minus_case] = integral_minus_core(limb_darkening_coefficients, rp_over_rs,
                                               z_over_rs[minus_case], 0.0, theta_2[minus_case],
                                               method=method, precision=precision)

    # flux_star
    starflux = np.zeros(len(z_over_rs))
    starflux[star_case] = integral_centred(limb_darkening_coefficients, 1, 0.0, ph[star_case], method=method)

    # flux_total
    total_flux = integral_centred(limb_darkening_coefficients, 1, 0.0, 2.0 * pi, method=method)

    return 1 - (2.0 / total_flux) * (plusflux + starflux - minsflux)
