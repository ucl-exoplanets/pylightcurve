import glob
import os
import math

import numpy as np
import matplotlib.pyplot as plt

from tasks import *

pi = np.pi

gauss30 = [
	[0.1028526528935588,	-0.0514718425553177],
	[0.1028526528935588,	0.0514718425553177],
	[0.1017623897484055,	-0.1538699136085835],
	[0.1017623897484055,	0.1538699136085835],
	[0.0995934205867953,	-0.2546369261678899],
	[0.0995934205867953,	0.2546369261678899],
	[0.0963687371746443,	-0.3527047255308781],
	[0.0963687371746443,	0.3527047255308781],
	[0.0921225222377861,	-0.4470337695380892],
	[0.0921225222377861,	0.4470337695380892],
	[0.0868997872010830,	-0.5366241481420199],
	[0.0868997872010830,	0.5366241481420199],
	[0.0807558952294202,	-0.6205261829892429],
	[0.0807558952294202,	0.6205261829892429],
	[0.0737559747377052,	-0.6978504947933158],
	[0.0737559747377052,	0.6978504947933158],
	[0.0659742298821805,	-0.7677774321048262],
	[0.0659742298821805,	0.7677774321048262],
	[0.0574931562176191,	-0.8295657623827684],
	[0.0574931562176191,	0.8295657623827684],
	[0.0484026728305941,	-0.8825605357920527],
	[0.0484026728305941,	0.8825605357920527],
	[0.0387991925696271,	-0.9262000474292743],
	[0.0387991925696271,	0.9262000474292743],
	[0.0287847078833234,	-0.9600218649683075],
	[0.0287847078833234,	0.9600218649683075],
	[0.0184664683110910,	-0.9836681232797472],
	[0.0184664683110910,	0.9836681232797472],
	[0.0079681924961666,	-0.9968934840746495],
	[0.0079681924961666,	0.9968934840746495]
]
gauss30 = np.swapaxes(gauss30,0,1)


def INTR(a1, a2, a3, a4, r):
    a0 = 1.0 - a1 - a2 - a3 - a4
    RR = (1.0 - r ** 2) ** (1.0 / 4)
    AA4 = - (2.0 / 4) * a0 * (RR ** 4.0)
    AA5 = - (2.0 / 5) * a1 * (RR ** 5.0)
    AA6 = - (2.0 / 6) * a2 * (RR ** 6.0)
    AA7 = - (2.0 / 7) * a3 * (RR ** 7.0)
    AA8 = - (2.0 / 8) * a4 * (RR ** 8.0)
    return AA4 + AA5 + AA6 + AA7 + AA8


def num(r,a1,a2,a3,a4,p,z):
	arccos	=	( -p**2 + z**2 +r**2 )/( 2.0*z*r )
	arccos	=	np.where(arccos>1,1,arccos)
	return ( 1.0 - a1*(1.0 - (1.0 - r**2)**(1.0/4)) - a2*(1.0 - (1.0 - r**2)**(1.0/2)) - a3*(1.0 - (1.0 - r**2)**(3.0/4)) - a4*(1.0 - (1.0 - r**2)) )*r*np.arccos( arccos )


def intr(a1,a2,a3,a4,p,z,r1,r2):
	gauss	=	gauss30
	r1		=	np.clip(r1, 0, 1)
	r2		=	np.clip(r2, 0, 1)
	x1		=	(r2-r1)/2
	x2		=	(r2+r1)/2
	x1, l	=	np.meshgrid(x1,gauss[1])
	x2, w	=	np.meshgrid(x2,gauss[0])
	result	=	np.sum( w*num(x1*l+x2,a1,a2,a3,a4,p,z), 0 )
	return ((r2-r1)/2)*result


def INTCENT(a1, a2, a3, a4, p, ww1, ww2):
    w1 = np.minimum(ww1, ww2)
    w2 = np.maximum(ww1, ww2)
    FINAL = ( INTR(a1, a2, a3, a4, p) - INTR(a1, a2, a3, a4, 0.0) ) * ( w2 - w1 )
    return FINAL


def INTPLUS0(a1, a2, a3, a4, p, z, ww1, ww2):
    if len(z) == 0:
        return []
    rr1 = z*np.cos(ww1) + np.sqrt(np.maximum(p**2-(z*np.sin(ww1))**2,0))
    rr1 = np.clip(rr1,0,1)
    rr2 = z*np.cos(ww2) + np.sqrt(np.maximum(p**2-(z*np.sin(ww2))**2,0))
    rr2 = np.clip(rr2,0,1)
    w1 = np.minimum(ww1, ww2)
    r1 = np.minimum(rr1, rr2)
    w2 = np.maximum(ww1, ww2)
    r2 = np.maximum(rr1, rr2)
    PARTA = INTR(a1, a2, a3, a4, 0.0) * ( w1 - w2 )
    PARTB = INTR(a1, a2, a3, a4, r1) * (  w2 )
    PARTC = INTR(a1, a2, a3, a4, r2) * ( -w1 )
    PARTD = intr(a1, a2, a3, a4, p, z, r1, r2)
    FINAL = PARTA + PARTB + PARTC + PARTD
    return FINAL


def INTPLUS(a1, a2, a3, a4, p, z, ww1, ww2):
	split = np.where(ww1*ww2<0)
	if len(split[0]) == 0:
		intplus = INTPLUS0(a1, a2, a3, a4, p, z, np.abs(ww1), np.abs(ww2))
		return intplus
	else:
		w1 = np.minimum(ww1, ww2)
		w2 = np.maximum(ww1, ww2)
		intplus = INTPLUS0(a1, a2, a3, a4, p, z, np.where(w1*w2<0,0,np.abs(w1)), np.abs(w2))
		intplus[split] = intplus[split] + INTPLUS0(a1, a2, a3, a4, p, z[split], 0, np.abs(w1[split]) )
		return intplus	


def INTMINS0(a1, a2, a3, a4, p, z, ww1, ww2):
    if len(z) == 0:
        return []
    rr1 = z*np.cos(ww1) - np.sqrt(np.maximum(p**2-(z*np.sin(ww1))**2,0))
    rr1 = np.clip(rr1,0,1)
    rr2 = z*np.cos(ww2) - np.sqrt(np.maximum(p**2-(z*np.sin(ww2))**2,0))
    rr2 = np.clip(rr2,0,1)
    w1 = np.minimum(ww1, ww2)
    r1 = np.minimum(rr1, rr2)
    w2 = np.maximum(ww1, ww2)
    r2 = np.maximum(rr1, rr2)
    PARTA = INTR(a1, a2, a3, a4, 0.0) * ( w1 - w2 )
    PARTB = INTR(a1, a2, a3, a4, r1) * ( -w1 )
    PARTC = INTR(a1, a2, a3, a4, r2) * (  w2 )
    PARTD = intr(a1, a2, a3, a4, p, z, r1, r2)
    FINAL = PARTA + PARTB + PARTC - PARTD
    return FINAL


def INTMINS(a1, a2, a3, a4, p, z, ww1, ww2):
	split = np.where(ww1*ww2<0)
	if len(split[0]) == 0:
		intmins =  INTMINS0(a1, a2, a3, a4, p, z, np.abs(ww1), np.abs(ww2))
		return intmins
	else:
		w1 = np.minimum(ww1, ww2)
		w2 = np.maximum(ww1, ww2)
		intmins = INTMINS0(a1, a2, a3, a4, p, z, np.where(w1*w2<0,0,np.abs(w1)), np.abs(w2))
		intmins[split] = intmins[split] + INTMINS0(a1, a2, a3, a4, p, z[split], 0, np.abs(w1[split]) )
		return intmins	




def single_model( ldcoeffs, r_ratio, xyz, tt):
    if len(tt) == 0:
        return np.array([])
    a1, a2, a3, a4 = ldcoeffs
    p = r_ratio
    ## projected distance
    fx, fy, fz = xyz
    z = np.sqrt(fy ** 2 + fz ** 2)
    ## cases
    case0 = np.where((fx > 0) & (z == 0) & (p <= 1))
    case1 = np.where((fx > 0) & (z < p) & (z + p <= 1))
    caseA = np.where((fx > 0) & (z < p) & (z + p > 1) & (p - z < 1))
    caseB = np.where((fx > 0) & (z < p) & (z + p > 1) & (p - z > 1))
    case2 = np.where((fx > 0) & (z == p) & (z + p <= 1))
    caseC = np.where((fx > 0) & (z == p) & (z + p > 1))
    case3 = np.where((fx > 0) & (z > p) & (z + p < 1))
    case4 = np.where((fx > 0) & (z > p) & (z + p == 1))
    case5 = np.where((fx > 0) & (z > p) & (z + p > 1) & (z ** 2 - p ** 2 < 1))
    case6 = np.where((fx > 0) & (z > p) & (z + p > 1) & (z ** 2 - p ** 2 == 1))
    case7 = np.where((fx > 0) & (z > p) & (z + p > 1) & (z ** 2 - p ** 2 > 1) & (z < 1 + p))
    pcase = np.concatenate((case1[0],case2[0],case3[0],case4[0],case5[0],caseA[0],caseC[0]))
    mcase = np.concatenate((case3[0],case4[0],case5[0],case6[0],case7[0]))
    scase = np.concatenate((case5[0],case6[0],case7[0],caseA[0],caseC[0]))
    ## cross points
    th = np.arcsin(np.where(p / z > 1.0, 1.0, p / z))
    arccos = np.clip((1.0 - p ** 2 + z ** 2) / (2.0 * z),-1,1)
    ph = np.arccos( arccos )
    ## flux_upper
    plusflux		=	np.zeros(len(z))
    theta_1			=	np.zeros(len(z))
    theta_1[case5]	=	ph[case5]
    theta_1[caseA]	=	ph[caseA]
    theta_1[caseC]	=	ph[caseC]
    theta_2			=	np.full_like(th,th)
    theta_2[case1]	=	pi
    theta_2[case2]	=	pi/2.0
    theta_2[caseA]	=	pi
    theta_2[caseC]	=	pi/2.0
    plusflux[pcase]	=	INTPLUS(a1, a2, a3, a4, p, z[pcase], theta_1[pcase], theta_2[pcase] )
    if len(case0[0]) > 0:
        plusflux[case0]	=	INTCENT(a1, a2, a3, a4, p, 0.0, pi )
    if len(caseB[0]) > 0:
        plusflux[caseB]	=	INTCENT(a1, a2, a3, a4, 1, 0.0, pi )
    ## flux_lower
    minsflux		=	np.zeros(len(z))
    theta_2			=	np.full_like(th,th)
    theta_2[case7]	=	ph[case7]
    minsflux[mcase]	=	INTMINS(a1, a2, a3, a4, p, z[mcase], 0.0, theta_2[mcase] )
    ## flux_star
    starflux		=	np.zeros(len(z))
    starflux[scase]	=	INTCENT(a1, a2, a3, a4, 1, 0.0, ph[scase] )
    ## flux_final
    finalflux		=	2.0 * ( plusflux + starflux - minsflux )
    # return
    F0 = INTCENT(a1, a2, a3, a4, 1, 0.0, 2.0 * pi)
    return 1 - finalflux/F0 



def polar(x,y):
	radius = np.sqrt(x**2+y**2)
	if x>0:
		if y>=0:
			angle =  np.arctan(y/x)
		else:
			angle = 2*pi+np.arctan(y/x)
	elif x<0:
		angle = np.arctan(y/x)+pi
	else:
		if y>=0:
			angle = pi/2.0
		else:
			angle = -pi/2.0
 	if angle > pi:
 		angle = angle -2.0*pi
	return (radius,angle)



def pilim(angle):
	if angle > pi:
		return angle - 2.0*pi
	elif angle < -pi:
		return angle + 2.0*pi
	else:
		return angle


# def polarnew(xx,yy):
# 	radius = np.sqrt(xx**2+yy**2)
# 	angles = []
# 	for (x,y) in np.swapaxes([xx,yy],0,1):
# 		if x>0:
# 			if y>=0:
# 				angle =  np.arctan(y/x)
# 			else:
# 				angle = 2*pi+np.arctan(y/x)
# 		elif x<0:
# 			angle = np.arctan(y/x)+pi
# 		else:
# 			if y>=0:
# 				angle = pi/2.0
# 			else:
# 				angle = -pi/2.0
# 		if angle > pi:
# 			angle = angle -2.0*pi
# 		angles.append(angle)
# 	return (radius,np.array(angles))
# 
# 
# def limitsnew(p1,p2,z1,z2,z12):
# 	star = []
# 	pls1 = []
# 	mns1 = []
# 	pls2 = []
# 	mns2 = []
# 	star.append(np.ones_like(z1)*pi)
# 	star.append(np.ones_like(z1)*(-pi))
# 	th1 = np.arcsin(p1/z1)
# 	pls1.append(np.where(z1<p1,pi,th1))
# 	pls1.append(np.where(z1<p1,-pi,-th1))
# 	mns1.append(np.where(z1<p1,0.0,th1))
# 	mns1.append(np.where(z1<p1,0.0,-th1))
# 	th2 = np.arcsin(p2/z2)
# 	pls2.append(np.where(z2<p2,pi,th2))
# 	pls2.append(np.where(z2<p2,-pi,-th2))
# 	mns2.append(np.where(z2<p2,0.0,th2))
# 	mns2.append(np.where(z2<p2,0.0,-th2))
# 	ph = np.arccos( (z12**2 - z1**2 - z2**2)/(-2.0*z1*z2) )
# 	ph[np.where( (z1 == 0) | (z2==0) )] = 0.0 
# 	th = np.arccos( (z2**2 - z1**2 - z12**2)/(-2.0*z1*z12) )
# 	ph[np.where( (z1 == 0) | (z12==0) )] = 0.0
# 	om = np.arccos( (p2**2 - p1**2 - z12**2)/(-2.0*p1*z12) )
# 	x1 = z1 - p1*np.cos( th + om )
# 	y1 = p1*np.sin( th + om )
# 	x2 = z1 - p1*np.cos( th - om )
# 	y2 = p1*np.sin( th - om )
# 	rA, thA = polarnew(x1,y1)
# 	rB, thB = polarnew(x2,y2)



def limits(p1,p2,z1,z2,z12):
	starscons = [-pi,pi]
	plpls1cons = []
	plmns1cons = []
	plpls2cons = []
	plmns2cons = []
	##
	if z1 < p1:
		plpls1cons.append(-pi)
		plpls1cons.append(pi)
	elif z1 == p1:
		plpls1cons.append(-pi/2)
		plpls1cons.append(pi/2)
	else:
		th1 = np.arcsin(p1/z1)
		plpls1cons.append(-th1)
		plpls1cons.append(th1)
		plmns1cons.append(-th1)
		plmns1cons.append(th1)
	##
	if z2 < p2:
		plpls2cons.append(-pi)
		plpls2cons.append(pi)
	elif z2 == p2:
		plpls2cons.append(-pi/2)
		plpls2cons.append(pi/2)
	else:
		th2 = np.arcsin(p2/z2)
		plpls2cons.append(-th2)
		plpls2cons.append(th2)
		plmns2cons.append(-th2)
		plmns2cons.append(th2)
	##
	if ( z1 == 0 or z2==0 or z12 + z2 == z1 or z1 + z12 == z2):
		ph = 0
	else:
		ph = np.arccos((z12**2 - z1**2 - z2**2)/(-2.0*z1*z2))
	if z12 + z2 == z1:
		th = 0
	elif z1 + z12 == z2:
		th = pi
	else:
		th = np.arccos(( z2**2 - z1**2 - z12**2 )/( -2.0*z1*z12 ))
	om = np.arccos(( p2**2 - p1**2 - z12**2)/( -2.0*p1*z12 ))
	x1 = z1 - p1*np.cos( th + om )
	y1 = p1*np.sin( th + om )
	x2 = z1 - p1*np.cos( th - om )
	y2 = p1*np.sin( th - om )
	rA, thA = polar(x1,y1)
	rB, thB = polar(x2,y2)
	if z1 <= p1:
		plpls1cons.append(thA)
		plpls1cons.append(thB)
	else:
		if   rA < np.sqrt( z1**2 - p1**2 ):
			plmns1cons.append(thA)
		if   rA > np.sqrt( z1**2 - p1**2 ):
			plpls1cons.append(thA)
		if   rB < np.sqrt( z1**2 - p1**2 ):
			plmns1cons.append(thB)
		if   rB > np.sqrt( z1**2 - p1**2 ):
			plpls1cons.append(thB)
	if z2 <= p2:
		plpls2cons.append(pilim(thA-ph))
		plpls2cons.append(pilim(thB-ph))
	else:
		if   rA < np.sqrt( z2**2 - p2**2 ):
			plmns2cons.append(pilim(thA-ph))
		if   rA > np.sqrt( z2**2 - p2**2 ):
			plpls2cons.append(pilim(thA-ph))
		if   rB < np.sqrt( z2**2 - p2**2 ):
			plmns2cons.append(pilim(thB-ph))
		if   rB > np.sqrt( z2**2 - p2**2 ):
			plpls2cons.append(pilim(thB-ph))
	##
	if z1 == 1.0 - p1:
		plpls1cons.append(0.0)
	if abs(1.0 - p1) < z1 < 1.0 + p1:
		th1 = np.arccos( ( 1.0 - p1**2 + z1**2 )/( 2.0*z1 ) )
		starscons.append(-th1)
		starscons.append(th1)
		if z1 <= p1:
			plpls1cons.append(-th1)
			plpls1cons.append(th1)
		else:
			if   np.sqrt( z1**2 - p1**2 ) > 1.0:
				plmns1cons.append(-th1)
				plmns1cons.append(th1)
			elif np.sqrt( z1**2 - p1**2 ) < 1.0:
				plpls1cons.append(-th1)
				plpls1cons.append(th1)
	##
	if z2 == 1.0 - p2:
		plpls2cons.append(0.0)
	if abs(1.0 - p2) < z2 < 1.0 + p2:
		th2 = np.arccos( ( 1.0 - p2**2 + z2**2 )/( 2.0*z2 ) )
		starscons.append(pilim(ph-th2))
		starscons.append(pilim(ph+th2))
		if z2 <= p2:
			plpls2cons.append(-th2)
			plpls2cons.append(th2)
		else:
			if   np.sqrt( z2**2 - p2**2 ) > 1.0:
				plmns2cons.append(-th2)
				plmns2cons.append(th2)
			elif np.sqrt( z2**2 - p2**2 ) < 1.0:
				plpls2cons.append(-th2)
				plpls2cons.append(th2)
	##
	starscons.sort()
	plpls1cons.sort()
	plmns1cons.sort()
	plpls2cons.sort()
	plmns2cons.sort()
	final_starscons = []
	final_plpls1cons = []
	final_plmns1cons = []
	final_plpls2cons = []
	final_plmns2cons = []
	###
	for i in range(len(starscons)-1):
		testu = (2.0/3)*( starscons[i+1] + 0.5*starscons[i] )
		testr = 1.0
		if testr**2 + z1**2 -2.0*testr*z1*np.cos(abs(testu-0)) < p1**2:
			if testr**2 + z2**2 -2.0*testr*z2*np.cos(abs(testu-ph)) < p2**2:
				final_starscons.append(starscons[i])
				final_starscons.append(starscons[i+1])
	for i in range(len(plpls1cons)-1):
		testu = (2.0/3)*( plpls1cons[i+1] + 0.5*plpls1cons[i] )
		testr = z1*np.cos(testu) + np.sqrt( p1**2 - (z1*np.sin(testu))**2 )
		if  testr < 1.0:
			if testr**2 + z2**2 -2.0*testr*z2*np.cos(abs(testu-ph)) < p2**2:
				final_plpls1cons.append(plpls1cons[i])
				final_plpls1cons.append(plpls1cons[i+1])
	if len(plmns1cons) > 0:
		for i in range(len(plmns1cons)-1):
			testu = (2.0/3)*( plmns1cons[i+1] + 0.5*plmns1cons[i] )
			testr = z1*np.cos(testu) - np.sqrt( p1**2 - (z1*np.sin(testu))**2 )
			if testr < 1.0:
				if testr**2 + z2**2 -2.0*testr*z2*np.cos(abs(testu-ph)) < p2**2:
					final_plmns1cons.append(plmns1cons[i])
					final_plmns1cons.append(plmns1cons[i+1])
	for i in range(len(plpls2cons)-1):
		testu = (2.0/3)*( plpls2cons[i+1] + 0.5*plpls2cons[i] )
		testr = z2*np.cos(testu) + np.sqrt( p2**2 - (z2*np.sin(testu))**2 )
		if  testr < 1.0:
			if testr**2 + z1**2 -2.0*testr*z1*np.cos(abs(testu+ph-0)) < p1**2:
				final_plpls2cons.append(plpls2cons[i])
				final_plpls2cons.append(plpls2cons[i+1])
	if len(plmns2cons) > 0:
		for i in range(len(plmns2cons)-1):
			testu = (2.0/3)*( plmns2cons[i+1] + 0.5*plmns2cons[i] )
			testr = z2*np.cos(testu) - np.sqrt( p2**2 - (z2*np.sin(testu))**2 )
			if testr < 1.0:
				if testr**2 + z1**2 -2.0*testr*z1*np.cos(abs(testu+ph-0)) < p1**2:
					final_plmns2cons.append(plmns2cons[i])
					final_plmns2cons.append(plmns2cons[i+1])
	##
	return (	final_starscons, 
				final_plpls1cons,
				final_plmns1cons,
				final_plpls2cons,
				final_plmns2cons 
				)

def snapshot(p1,p2,z1,z2,z12,ii):
	z1 = z1[ii]
	z2 = z2[ii]
	z12 = z12[ii]
	if ( z1 == 0 or z2==0 or z12 + z2 == z1 or z1 + z12 == z2):
		ph = 0
	else:
		ph = np.arccos((z12**2 - z1**2 - z2**2)/(-2.0*z1*z2))
	fig = plt.figure()
	ax1 = fig.add_subplot(1,1,1)
	circle0=plt.Circle((0,0), 1, color=[0.9,0.9,0], fill=False)
	ax1.add_patch(circle0)
	circle1=plt.Circle((z1,0), p1, color='r', fill=False)
	ax1.add_patch(circle1)
	circle2=plt.Circle((z2*np.cos(ph),z2*np.sin(ph)), p2, color='k', fill=False)
	ax1.add_patch(circle2)
	plt.xlim(-2,2)
	plt.ylim(-2,2)
	plt.axhline(0)
	plt.axvline(0)
	plt.plot([0,2*z2*np.cos(ph)],[0,2*z2*np.sin(ph)],'k-')
	plt.show()


def double_model( ldcoeffs, r_ratio1, r_ratio2, xyz1, xyz2, tt):
	a1, a2, a3, a4 = ldcoeffs
	p1, p2 = r_ratio1, r_ratio2
	## projected distance
	fx1, fy1, fz1 = xyz1
	z1 = np.sqrt(fy1 ** 2 + fz1 ** 2)
	fx2, fy2, fz2 = xyz2
	z2 = np.sqrt(fy2 ** 2 + fz2 ** 2)
	z12=np.sqrt((fy1-fy2)**2+(fz1-fz2)**2)
	## cases
	case1 = np.where( z12 <= p1 - p2 )
	case2 = np.where( z12 <= p2 - p1 )
	case3 = np.where((fx1 > 0) & (fx2 > 0) & (z1 < 1 + p1) & (z2 < 1 + p2) & (z12 > abs(p1-p2)) & (z12 < p1+p2) )
	if len(case3[0]) > 1:
		## flux overlap
		theta_star  = []
		theta_plus1 = []
		theta_mins1 = []
		theta_plus2 = []
		theta_mins2 = []
		for i in case3[0]:
			lims = limits(p1,p2,z1[i],z2[i],z12[i])
			theta_star.append( lims[0])
			theta_plus1.append(lims[1])
			theta_mins1.append(lims[2])
			theta_plus2.append(lims[3])
			theta_mins2.append(lims[4])
		theta_star	=	np.array( map(None,*theta_star), dtype=float )
		theta_plus1	=	np.array( map(None,*theta_plus1), dtype=float )
		theta_mins1	=	np.array( map(None,*theta_mins1), dtype=float )
		theta_plus2	=	np.array( map(None,*theta_plus2), dtype=float )
		theta_mins2	=	np.array( map(None,*theta_mins2), dtype=float )
		starflux	=	np.zeros(len(z1))
		for i in range(0,len(theta_star),2):
			subcase = np.where(np.isnan(theta_star[i])==False)
			case3sub = case3[0][subcase]
			starflux[case3sub] = starflux[case3sub] + INTCENT(a1, a2, a3, a4, 1, theta_star[i][subcase], theta_star[i+1][subcase] )
		plusflux1	=	np.zeros(len(z1))
		for i in range(0,len(theta_plus1),2):
			subcase = np.where(np.isnan(theta_plus1[i])==False)
			case3sub = case3[0][subcase]
			plusflux1[case3sub] = plusflux1[case3sub] + INTPLUS(a1, a2, a3, a4, p1, z1[case3sub], theta_plus1[i][subcase], theta_plus1[i+1][subcase] )
		plusflux2	=	np.zeros(len(z1))
		for i in range(0,len(theta_plus2),2):
			subcase = np.where(np.isnan(theta_plus2[i])==False)
			case3sub = case3[0][subcase]
			plusflux2[case3sub] = plusflux2[case3sub] + INTPLUS(a1, a2, a3, a4, p2, z2[case3sub], theta_plus2[i][subcase], theta_plus2[i+1][subcase] )
		minsflux1	=	np.zeros(len(z1))
		for i in range(0,len(theta_mins1),2):
			subcase = np.where(np.isnan(theta_mins1[i])==False)
			case3sub = case3[0][subcase]
			minsflux1[case3sub] = minsflux1[case3sub] + INTMINS(a1, a2, a3, a4, p1, z1[case3sub], theta_mins1[i][subcase], theta_mins1[i+1][subcase] )
		minsflux2	=	np.zeros(len(z1))
		for i in range(0,len(theta_mins2),2):
			subcase = np.where(np.isnan(theta_mins2[i])==False)
			case3sub = case3[0][subcase]
			minsflux2[case3sub] = minsflux2[case3sub] + INTMINS(a1, a2, a3, a4, p2, z2[case3sub], theta_mins2[i][subcase], theta_mins2[i+1][subcase] )
	else:
		plusflux1 = 0
		plusflux2 = 0
		starflux  = 0
		minsflux1 = 0
		minsflux2 = 0
	F0					=	INTCENT(a1, a2, a3, a4, 1, 0.0, 2.0 * pi)
	single2				=	single_model( ldcoeffs, r_ratio2, xyz2, tt )
	single1				=	single_model( ldcoeffs, r_ratio1, xyz1, tt )
	# flux_final
	finalflux			=	2.0 - single1 - single2
	finalflux[case1]	=	1.0 - single1[case1]
	finalflux[case2]	=	1.0 - single2[case2]
	F0					=	INTCENT(a1, a2, a3, a4, 1, 0.0, 2.0 * pi)
	finalflux			=	finalflux - (plusflux1 + plusflux2 + starflux - minsflux1 - minsflux2)/F0
	# return
	return 1 - finalflux





def transit(ldcoeffs, RpRs, P, a, e, i, w, W, T, tt):
	xyz = position(P, a, e, i * pi / 180, w * pi / 180, W * pi / 180, T, tt )
	return single_model(ldcoeffs, RpRs, xyz, tt)

def eclipse(FpFs, RpRs, P, a, e, i, w, W, T, tt):
	xyz = position(P, -a/RpRs, e, i * pi / 180, w * pi / 180, W * pi / 180, T, tt )
	return (1.0 + FpFs*single_model((0,0,0,0), 1/RpRs, xyz, tt))/(1.0 + FpFs)


def double_transit(ldcoeffs, RpRs1, RpRs2, P1, P2, a1, a2, e1, e2, i1, i2, w1, w2, W1, W2, T1, T2, tt):
	xyz1 = position(P1, a1, e1, i1 * pi / 180, w1 * pi / 180, W1 * pi / 180, T1, tt )
	xyz2 = position(P2, a2, e2, i2 * pi / 180, w2 * pi / 180, W2 * pi / 180, T2, tt )
	return double_model(ldcoeffs, RpRs1, RpRs2, xyz1, xyz2, tt)



def binary_caseA( lratio, ld_prm, r_prm, ld_sec, r_sec, P_sec, a_sec, e_sec, i_sec, w_sec, W_sec, T_sec, r_pln, P_pln, a_pln, e_pln, i_pln, w_pln, W_pln, T_pln, tt ):
	xyz_sec = position(P_sec, a_sec, e_sec, i_sec * pi / 180, w_sec * pi / 180, W_sec * pi / 180, T_sec, tt )
	xyz_pln = position(P_pln, a_pln, e_pln, i_pln * pi / 180, w_pln * pi / 180, W_pln * pi / 180, T_pln, tt )
	#primary
	xyz1 = [ xyz_sec[0]/r_prm, xyz_sec[1]/r_prm, xyz_sec[2]/r_prm ]
	xyz2 = [ xyz_pln[0]/r_prm, xyz_pln[1]/r_prm, xyz_pln[2]/r_prm ]
 	primary = double_model(ld_prm, r_sec/r_prm, r_pln/r_prm, xyz1, xyz2, tt)
	# secondary
	xyz1 = [ -xyz_sec[0]/r_sec, -xyz_sec[1]/r_sec, -xyz_sec[2]/r_sec ]
	xyz2 = [ (xyz_pln[0]-xyz_sec[0])/r_sec, (xyz_pln[1]-xyz_sec[1])/r_sec, (xyz_pln[2]-xyz_sec[2])/r_sec ]
	secondary = double_model(ld_sec, r_prm/r_sec, r_pln/r_sec, xyz1, xyz2, tt)
	return (primary*lratio+secondary)/(1.0+lratio)



def binary_caseB( lratio, ld_prm, r_prm, ld_sec, r_sec, P_sec, a_sec, e_sec, i_sec, w_sec, W_sec, T_sec, r_pln, P_pln, a_pln, e_pln, i_pln, w_pln, W_pln, T_pln, tt ):
	xyz_sec = position(P_sec, a_sec, e_sec, i_sec * pi / 180, w_sec * pi / 180, W_sec * pi / 180, T_sec, tt )
	xyz_pln = position(P_pln, a_pln, e_pln, i_pln * pi / 180, w_pln * pi / 180, W_pln * pi / 180, T_pln, tt )
	#primary
	xyz1 = [ xyz_sec[0]/r_prm, xyz_sec[1]/r_prm, xyz_sec[2]/r_prm ]
	xyz2 = [ (xyz_sec[0]+xyz_pln[0])/r_prm, (xyz_sec[1]+xyz_pln[1])/r_prm, (xyz_sec[2]+xyz_pln[2])/r_prm ]
 	primary = double_model(ld_prm, r_sec/r_prm, r_pln/r_prm, xyz1, xyz2, tt)
	# secondary
	xyz1 = [ -xyz_sec[0]/r_sec, -xyz_sec[1]/r_sec, -xyz_sec[2]/r_sec ]
	xyz2 = [ xyz_pln[0]/r_sec, xyz_pln[1]/r_sec, xyz_pln[2]/r_sec ]
	secondary = double_model(ld_sec, r_prm/r_sec, r_pln/r_sec, xyz1, xyz2, tt)
	return (primary*lratio+secondary)/(1.0+lratio)



def binary( lratio, ld_prm, r_prm, ld_sec, r_sec, P_sec, a_sec, e_sec, i_sec, w_sec, W_sec, T_sec, tt ):
	xyz_sec = position(P_sec, a_sec, e_sec, i_sec * pi / 180, w_sec * pi / 180, W_sec * pi / 180, T_sec, tt )
	#primary
	xyz1 = [ xyz_sec[0]/r_prm, xyz_sec[1]/r_prm, xyz_sec[2]/r_prm ]
 	primary = single_model(ld_prm, r_sec/r_prm, xyz1, tt)
	# secondary
	xyz1 = [ -xyz_sec[0]/r_sec, -xyz_sec[1]/r_sec, -xyz_sec[2]/r_sec ]
	secondary = single_model(ld_sec, r_prm/r_sec, xyz1, tt)
	return (primary*lratio+secondary)/(1.0+lratio)


