import numpy as np

pi = np.pi

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
gausstab = np.swapaxes(gauss30, 0, 1)


def integral_r(a1, a2, a3, a4, r):
    musq = 1 - r * r
    mu = np.sqrt(musq)
    return ((-2.0 / 5) * a2 * np.sqrt(mu) - (1.0 / 3) * a1 * mu + (1.0 / 2) * (-1 + a1 + a2)) * musq


def num(r, a1, a2, a3, a4, rprs, z):
    rsq = r * r
    mu = np.sqrt(1.0 - rsq)
    return (1.0 - a1 * (1 - mu) - a2 * (1.0 - np.sqrt(mu))) \
        * r * np.arccos(np.minimum((-rprs ** 2 + z * z + rsq) / (2.0 * z * r), 1.0))


def integral_r_f(a1, a2, a3, a4, rprs, z, r1, r2):
    x1 = (r2 - r1) / 2
    x2 = (r2 + r1) / 2
    return x1 * np.sum(gausstab[0][:, None] * num(x1[None, :] * gausstab[1][:, None] + x2[None, :],
                                                  a1, a2, a3, a4, rprs, z), 0)


def integral_centred(a1, a2, a3, a4, rprs, ww1, ww2):
    return (integral_r(a1, a2, a3, a4, rprs) - integral_r(a1, a2, a3, a4, 0.0)) * np.abs(ww2 - ww1)


def integral_plus_core(a1, a2, a3, a4, rprs, z, ww1, ww2):
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
    parta = integral_r(a1, a2, a3, a4, 0.0) * (w1 - w2)
    partb = integral_r(a1, a2, a3, a4, r1) * w2
    partc = integral_r(a1, a2, a3, a4, r2) * (-w1)
    partd = integral_r_f(a1, a2, a3, a4, rprs, z, r1, r2)
    return parta + partb + partc + partd


def integral_plus(a1, a2, a3, a4, rprs, z, ww1, ww2):
    split = np.where(ww1 * ww2 < 0)
    if len(split[0]) == 0:
        return integral_plus_core(a1, a2, a3, a4, rprs, z, np.abs(ww1), np.abs(ww2))
    else:
        w1 = np.minimum(ww1, ww2)
        w2 = np.maximum(ww1, ww2)
        intplus = integral_plus_core(a1, a2, a3, a4, rprs, z, np.where(w1 * w2 < 0, 0, np.abs(w1)), np.abs(w2))
        intplus[split] = intplus[split] + integral_plus_core(a1, a2, a3, a4, rprs, z[split], 0, np.abs(w1[split]))
        return intplus


def integral_minus_core(a1, a2, a3, a4, rprs, z, ww1, ww2):
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
    parta = integral_r(a1, a2, a3, a4, 0.0) * (w1 - w2)
    partb = integral_r(a1, a2, a3, a4, r1) * (-w1)
    partc = integral_r(a1, a2, a3, a4, r2) * w2
    partd = integral_r_f(a1, a2, a3, a4, rprs, z, r1, r2)
    return parta + partb + partc - partd


def integral_minus(a1, a2, a3, a4, rprs, z, ww1, ww2):
    split = np.where(ww1 * ww2 < 0)
    if len(split[0]) == 0:
        return integral_minus_core(a1, a2, a3, a4, rprs, z, np.abs(ww1), np.abs(ww2))
    else:
        w1 = np.minimum(ww1, ww2)
        w2 = np.maximum(ww1, ww2)
        intmins = integral_minus_core(a1, a2, a3, a4, rprs, z, np.where(w1 * w2 < 0, 0, np.abs(w1)), np.abs(w2))
        intmins[split] = intmins[split] + integral_minus_core(a1, a2, a3, a4, rprs, z[split], 0, np.abs(w1[split]))
        return intmins


def flux_drop(limb_darkening_coefficients, rp_over_rs, z_over_rs):

    if len(z_over_rs) == 0:
        return np.array([])

    a1, a2, a3, a4 = limb_darkening_coefficients

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
    plusflux[plus_case] = integral_plus_core(a1, a2, a3, a4, rp_over_rs,
                                             z_over_rs[plus_case], theta_1[plus_case], theta_2[plus_case])
    if len(case0[0]) > 0:
        plusflux[case0] = integral_centred(a1, a2, a3, a4, rp_over_rs, 0.0, pi)
    if len(caseb[0]) > 0:
        plusflux[caseb] = integral_centred(a1, a2, a3, a4, 1, 0.0, pi)

    # flux_lower
    minsflux = np.zeros(len(z_over_rs))
    minsflux[minus_case] = integral_minus_core(a1, a2, a3, a4, rp_over_rs,
                                               z_over_rs[minus_case], 0.0, theta_2[minus_case])

    # flux_star
    starflux = np.zeros(len(z_over_rs))
    starflux[star_case] = integral_centred(a1, a2, a3, a4, 1, 0.0, ph[star_case])

    # flux_total
    total_flux = integral_centred(a1, a2, a3, a4, 1, 0.0, 2.0 * pi)

    return 1 - (2.0 / total_flux) * (plusflux + starflux - minsflux)
