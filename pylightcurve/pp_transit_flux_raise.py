__all__ = ['pp_transit_flux_raise']

import numpy as np

from transit_flux_drop import *


def get_circle_edge_points(z, r):

    point1_theta = np.ones_like(z) * (-np.pi)
    point1_rho = np.ones_like(z) * (r - z)

    test = np.where(z >= r)
    point1_theta[test] = - np.arcsin(r / z[test])
    point1_rho[test] = np.sqrt(z[test] * z[test] - r * r)

    point2_theta = - np.ones_like(point1_theta) * point1_theta
    point2_rho = np.ones_like(point1_rho) * point1_rho

    return point1_theta, point1_rho, point2_theta, point2_rho

z_test = np.arange(0, 1.5, 0.00001)
import matplotlib.pyplot as plt
for i in get_circle_edge_points(z_test, 0.1):
    plt.plot(z_test, i, 'o', ms=3)

plt.show()


# def solve_two_circles(theta_1, rho_1, r_1, theta_2, rho_2, r2):

#
#
# def pp_transit_flux_raise(method, limb_darkening_coefficients,
#                           rp_over_rs_1, z_over_rs_1, rp_over_rs_2, z_over_rs_2, z_over_rs_12, precision=3):
#
#     if len(z_over_rs_1) == 0:
#         return np.array([])
#
#     if (z_over_rs_1 >= 1 + rp_over_rs_1 or z_over_rs_2 >= 1 + rp_over_rs_2
#        or z_over_rs_12 >= rp_over_rs_1 + rp_over_rs_2):
#         return 0.0 * np.ones_like(z_over_rs_1)
#
#     if z_over_rs_12 <= rp_over_rs_1 - rp_over_rs_2:
#         return transit_flux_drop(method, limb_darkening_coefficients, rp_over_rs_2, z_over_rs_2, precision)
#
#
#
#
# def DFB(a1,a2,a3,a4,rp_over_rs_1,z_over_rs_1,rp_over_rs_2,z_over_rs_2,z_over_rs_12):
#     ###
#     if ( z_over_rs_1 >= 1 + rp_over_rs_1 or z_over_rs_2 >= 1 + rp_over_rs_2 or z_over_rs_12 >= rp_over_rs_1 + rp_over_rs_2):
#         return 0.0
#     if z_over_rs_12 <= rp_over_rs_1 - rp_over_rs_2:
#         return transit_flux_drop(a1,a2,a3,a4,rp_over_rs_2,z_over_rs_2)
#     ###
#     F0 = 2*pi*( INTR(a1,a2,a3,a4,1.0) - INTR(a1,a2,a3,a4,0.0) )
#     ###
#     starscons = [[-pi,1.0],[pi,1.0]]
#     plpls1cons = []
#     plmns1cons = []
#     plpls2cons = []
#     plmns2cons = []
#     ##
#     if z_over_rs_1 < rp_over_rs_1:
#         plpls1cons.append([-pi,rp_over_rs_1-z_over_rs_1])
#         plpls1cons.append([pi,rp_over_rs_1-z_over_rs_1])
#     elif z_over_rs_1 == rp_over_rs_1:
#         plpls1cons.append([-pi/2,0.0])
#         plpls1cons.append([pi/2,0.0])
#     else:
#         th1 = np.arcsin(rp_over_rs_1/z_over_rs_1)
#         ro1 = np.sqrt( z_over_rs_1**2 - rp_over_rs_1**2 )
#         plpls1cons.append([-th1,ro1])
#         plpls1cons.append([th1,ro1])
#         plmns1cons.append([-th1,ro1])
#         plmns1cons.append([th1,ro1])
#     ##
#     if z_over_rs_2 < rp_over_rs_2:
#         plpls2cons.append([-pi,rp_over_rs_2-z_over_rs_2])
#         plpls2cons.append([pi,rp_over_rs_2-z_over_rs_2])
#     elif z_over_rs_2 == rp_over_rs_2:
#         plpls2cons.append([-pi/2,0.0])
#         plpls2cons.append([pi/2,0.0])
#     else:
#         th2 = np.arcsin(rp_over_rs_2/z_over_rs_2)
#         ro2 = np.sqrt( z_over_rs_2**2 - rp_over_rs_2**2 )
#         plpls2cons.append([-th2,ro2])
#         plpls2cons.append([th2,ro2])
#         plmns2cons.append([-th2,ro2])
#         plmns2cons.append([th2,ro2])
#     ###
#     if ( z_over_rs_1 == 0 or z_over_rs_2==0 or z_over_rs_12 + z_over_rs_2 == z_over_rs_1 or z_over_rs_1 + z_over_rs_12 == z_over_rs_2):
#         ph = 0
#     else:
#         ph = np.arccos((z_over_rs_12**2 - z_over_rs_1**2 - z_over_rs_2**2)/(-2.0*z_over_rs_1*z_over_rs_2))
#     if z_over_rs_12 + z_over_rs_2 == z_over_rs_1:
#         th = 0
#     elif z_over_rs_1 + z_over_rs_12 == z_over_rs_2:
#         th = pi
#     else:
#         th = np.arccos(( z_over_rs_2**2 - z_over_rs_1**2 - z_over_rs_12**2 )/( -2.0*z_over_rs_1*z_over_rs_12 ))
#     om = np.arccos(( rp_over_rs_2**2 - rp_over_rs_1**2 - z_over_rs_12**2)/( -2.0*rp_over_rs_1*z_over_rs_12 ))
#     x1 = z_over_rs_1 - rp_over_rs_1*np.cos( th + om )
#     y1 = rp_over_rs_1*np.sin( th + om )
#     x2 = z_over_rs_1 - rp_over_rs_1*np.cos( th - om )
#     y2 = rp_over_rs_1*np.sin( th - om )
#     rA, thA = polar(x1,y1)
#     rB, thB = polar(x2,y2)
#     ##
#     if z_over_rs_1 <= rp_over_rs_1:
#         plpls1cons.append([thA,rA])
#         plpls1cons.append([thB,rB])
#     else:
#         if   rA < np.sqrt( z_over_rs_1**2 - rp_over_rs_1**2 ):
#             plmns1cons.append([thA,rA])
#         if   rA > np.sqrt( z_over_rs_1**2 - rp_over_rs_1**2 ):
#             plpls1cons.append([thA,rA])
#         if   rB < np.sqrt( z_over_rs_1**2 - rp_over_rs_1**2 ):
#             plmns1cons.append([thB,rB])
#         if   rB > np.sqrt( z_over_rs_1**2 - rp_over_rs_1**2 ):
#             plpls1cons.append([thB,rB])
#     ##
#     if z_over_rs_2 <= rp_over_rs_2:
#         plpls2cons.append([pilim(thA-ph),rA])
#         plpls2cons.append([pilim(thB-ph),rB])
#     else:
#         if   rA < np.sqrt( z_over_rs_2**2 - rp_over_rs_2**2 ):
#             plmns2cons.append([pilim(thA-ph),rA])
#         if   rA > np.sqrt( z_over_rs_2**2 - rp_over_rs_2**2 ):
#             plpls2cons.append([pilim(thA-ph),rA])
#         if   rB < np.sqrt( z_over_rs_2**2 - rp_over_rs_2**2 ):
#             plmns2cons.append([pilim(thB-ph),rB])
#         if   rB > np.sqrt( z_over_rs_2**2 - rp_over_rs_2**2 ):
#             plpls2cons.append([pilim(thB-ph),rB])
#     ###
#     if z_over_rs_1 == 1.0 - rp_over_rs_1:
#         plpls1cons.append([0.0,1.0])
#     if abs(1.0 - rp_over_rs_1) < z_over_rs_1 < 1.0 + rp_over_rs_1:
#         th1 = np.arccos( ( 1.0 - rp_over_rs_1**2 + z_over_rs_1**2 )/( 2.0*z_over_rs_1 ) )
#         starscons.append([-th1,1.0])
#         starscons.append([th1,1.0])
#         if z_over_rs_1 <= rp_over_rs_1:
#             plpls1cons.append([-th1,1.0])
#             plpls1cons.append([th1,1.0])
#         else:
#             if   np.sqrt( z_over_rs_1**2 - rp_over_rs_1**2 ) > 1.0:
#                 plmns1cons.append([-th1,1.0])
#                 plmns1cons.append([th1,1.0])
#             elif np.sqrt( z_over_rs_1**2 - rp_over_rs_1**2 ) < 1.0:
#                 plpls1cons.append([-th1,1.0])
#                 plpls1cons.append([th1,1.0])
#     ##
#     if z_over_rs_2 == 1.0 - rp_over_rs_2:
#         plpls2cons.append([0.0,1.0])
#     if abs(1.0 - rp_over_rs_2) < z_over_rs_2 < 1.0 + rp_over_rs_2:
#         th2 = np.arccos( ( 1.0 - rp_over_rs_2**2 + z_over_rs_2**2 )/( 2.0*z_over_rs_2 ) )
#         starscons.append([pilim(ph-th2),1.0])
#         starscons.append([pilim(ph+th2),1.0])
#         if z_over_rs_2 <= rp_over_rs_2:
#             plpls2cons.append([-th2,1.0])
#             plpls2cons.append([th2,1.0])
#         else:
#             if   np.sqrt( z_over_rs_2**2 - rp_over_rs_2**2 ) > 1.0:
#                 plmns2cons.append([-th2,1.0])
#                 plmns2cons.append([th2,1.0])
#             elif np.sqrt( z_over_rs_2**2 - rp_over_rs_2**2 ) < 1.0:
#                 plpls2cons.append([-th2,1.0])
#                 plpls2cons.append([th2,1.0])
#     ###
#     starscons.sort()
#     plpls1cons.sort()
#     plmns1cons.sort()
#     plpls2cons.sort()
#     plmns2cons.sort()
#     for i in range(len(starscons)-1):
#         if starscons[i][0]*starscons[i+1][0] < 0:
#             if starscons[i][0] + starscons[i+1][0] != 0:
#                 starscons.append([0.0,1.0])
#                 starscons.sort()
#                 break
#     for i in range(len(plpls1cons)-1):
#         if plpls1cons[i][0]*plpls1cons[i+1][0] < 0:
#             if plpls1cons[i][0] + plpls1cons[i+1][0] != 0:
#                 plpls1cons.append([0.0,rp_over_rs_1+z_over_rs_1])
#                 plpls1cons.sort()
#                 break
#     if len(plmns1cons) > 0:
#         for i in range(len(plmns1cons)-1):
#             if plmns1cons[i][0]*plmns1cons[i+1][0] < 0:
#                 if plmns1cons[i][0] + plmns1cons[i+1][0] != 0:
#                     plmns1cons.append([0.0,z_over_rs_1-rp_over_rs_1])
#                     plmns1cons.sort()
#                     break
#     for i in range(len(plpls2cons)-1):
#         if plpls2cons[i][0]*plpls2cons[i+1][0] < 0:
#             if plpls2cons[i][0] + plpls2cons[i+1][0] != 0:
#                 plpls2cons.append([0.0,rp_over_rs_2+z_over_rs_2])
#                 plpls2cons.sort()
#                 break
#     if len(plmns2cons) > 0:
#         for i in range(len(plmns2cons)-1):
#             if plmns2cons[i][0]*plmns2cons[i+1][0] < 0:
#                 if plmns2cons[i][0] + plmns2cons[i+1][0] != 0:
#                     plmns2cons.append([0.0,z_over_rs_2-rp_over_rs_2])
#                     plmns2cons.sort()
#                     break
#     ###
#     starsflux = 0.0
#     for i in range(len(starscons)-1):
#         testu = (2.0/3)*( starscons[i+1][0] + 0.5*starscons[i][0] )
#         testr = 1.0
#         if polardist(testr,testu,z_over_rs_1,0) < rp_over_rs_1**2:
#             if polardist(testr,testu,z_over_rs_2,ph) < rp_over_rs_2**2:
#                 starsflux = starsflux + INT(a1,a2,a3,a4,rp_over_rs_1,z_over_rs_1,'stars',starscons[i],starscons[i+1])
#     plpls1flux = 0.0
#     for i in range(len(plpls1cons)-1):
#         testu = (2.0/3)*( plpls1cons[i+1][0] + 0.5*plpls1cons[i][0] )
#         testr = z_over_rs_1*np.cos(testu) + np.sqrt( rp_over_rs_1**2 - (z_over_rs_1*np.sin(testu))**2 )
#         if  testr < 1.0:
#             if polardist(testr,testu,z_over_rs_2,ph) < rp_over_rs_2**2:
#                 plpls1flux = plpls1flux + INT(a1,a2,a3,a4,rp_over_rs_1,z_over_rs_1,'plpls',plpls1cons[i],plpls1cons[i+1])
#     plmns1flux = 0.0
#     if len(plmns1cons) > 0:
#         for i in range(len(plmns1cons)-1):
#             testu = (2.0/3)*( plmns1cons[i+1][0] + 0.5*plmns1cons[i][0] )
#             testr = z_over_rs_1*np.cos(testu) - np.sqrt( rp_over_rs_1**2 - (z_over_rs_1*np.sin(testu))**2 )
#             if testr < 1.0:
#                 if polardist(testr,testu,z_over_rs_2,ph) < rp_over_rs_2**2:
#                     plmns1flux = plmns1flux + INT(a1,a2,a3,a4,rp_over_rs_1,z_over_rs_1,'plmns',plmns1cons[i],plmns1cons[i+1])
#     plpls2flux = 0.0
#     for i in range(len(plpls2cons)-1):
#         testu = (2.0/3)*( plpls2cons[i+1][0] + 0.5*plpls2cons[i][0] )
#         testr = z_over_rs_2*np.cos(testu) + np.sqrt( rp_over_rs_2**2 - (z_over_rs_2*np.sin(testu))**2 )
#         if  testr < 1.0:
#             if polardist(testr,testu+ph,z_over_rs_1,0) < rp_over_rs_1**2:
#                 plpls2flux = plpls2flux + INT(a1,a2,a3,a4,rp_over_rs_2,z_over_rs_2,'plpls',plpls2cons[i],plpls2cons[i+1])
#     plmns2flux = 0.0
#     if len(plmns2cons) > 0:
#         for i in range(len(plmns2cons)-1):
#             testu = (2.0/3)*( plmns2cons[i+1][0] + 0.5*plmns2cons[i][0] )
#             testr = z_over_rs_2*np.cos(testu) - np.sqrt( rp_over_rs_2**2 - (z_over_rs_2*np.sin(testu))**2 )
#             if testr < 1.0:
#                 if polardist(testr,testu+ph,z_over_rs_1,0) < rp_over_rs_1**2:
#                     plmns2flux = plmns2flux + INT(a1,a2,a3,a4,rp_over_rs_2,z_over_rs_2,'plmns',plmns2cons[i],plmns2cons[i+1])
#     ###
#     return (starsflux + plpls1flux - plmns1flux + plpls2flux - plmns2flux)/F0
# ### THE FLUX-RISE FUNCTION
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#

