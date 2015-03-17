import numpy as np
import scipy.integrate as spi
import matplotlib.pyplot as plt
pi=np.pi






def INTR(a1,a2,a3,a4,r):
	a0 = 1.0 - a1 - a2 - a3 - a4
	RR = (1.0-r**2)**(1.0/4)
	AA4 = - (2.0/4)*a0*(RR**4.0)
	AA5 = - (2.0/5)*a1*(RR**5.0)
	AA6 = - (2.0/6)*a2*(RR**6.0)
	AA7 = - (2.0/7)*a3*(RR**7.0)
	AA8 = - (2.0/8)*a4*(RR**8.0)
	return	 AA4 + AA5 + AA6 + AA7 + AA8

def num(r,a1,a2,a3,a4,p,z):
	return ( 1.0 - a1*(1.0 - (1.0 - r**2)**(1.0/4)) - a2*(1.0 - (1.0 - r**2)**(1.0/2)) - a3*(1.0 - (1.0 - r**2)**(3.0/4)) - a4*(1.0 - (1.0 - r**2)) )*r*np.arccos( ( -p**2 + z**2 +r**2 )/( 2.0*z*r ) )

def intr0(a1,a2,a3,a4,p,z,r1,r2):
	return spi.quad(num, r1, r2, args=(a1,a2,a3,a4,p,z))[0]

intr = np.vectorize(intr0)




def INTCENT(a1,a2,a3,a4,p,z,ww1,ww2):
	w1 = np.min([ww1,ww2],0)
	w2 = np.max([ww1,ww2],0)
	return ( INTR(a1,a2,a3,a4,p) - INTR(a1,a2,a3,a4,0.0) )*( w2 - w1 )



def INTPLUS(a1,a2,a3,a4,p,z,ww1,rr1,ww2,rr2):
	if len(z) == 0:
		return []
	w1 = np.min([ww1,ww2],0)
	r1 = np.min([rr1,rr2],0)
	w2 = np.max([ww1,ww2],0)
	r2 = np.max([rr1,rr2],0)
	PARTA = INTR(a1,a2,a3,a4,0.0)*( w1 - w2 )
	PARTB = INTR(a1,a2,a3,a4,r1 )*(  w2 )
	PARTC = INTR(a1,a2,a3,a4,r2 )*( -w1 )
	PARTD = intr(a1,a2,a3,a4,p,z,r1,r2)
	return  PARTA + PARTB + PARTC + PARTD



def INTMINS(a1,a2,a3,a4,p,z,ww1,rr1,ww2,rr2):
	if len(z) == 0:
		return []
	w1 = np.min([ww1,ww2],0)
	r1 = np.min([rr1,rr2],0)
	w2 = np.max([ww1,ww2],0)
	r2 = np.max([rr1,rr2],0)
	PARTA = INTR(a1,a2,a3,a4,0.0)*( w1 - w2 )
	PARTB = INTR(a1,a2,a3,a4,r1 )*( -w1 )
	PARTC = INTR(a1,a2,a3,a4,r2 )*(  w2 )
	PARTD = intr(a1,a2,a3,a4,p,z,r1,r2)
	return  PARTA + PARTB + PARTC - PARTD



def position(P,A,E,I,W,WW,T0,tt):
	#
	if W<pi/2:
		AA=1.0*pi/2-W
	else:
		AA=5.0*pi/2-W
	BB=2*np.arctan(np.sqrt((1-E)/(1+E))*np.tan(AA/2))
	if BB < 0:
		BB = 2*pi + BB
	T0=T0-(P/2.0/pi)*(BB-E*np.sin(BB))
	#
	M=(tt-T0-np.int_((tt-T0)/P)*P)*2.0*pi/P
	u0=M
	stop=False
	while stop==False:
		u1=u0-(u0-E*np.sin(u0)-M)/(1-E*np.cos(u0))
		stop = (np.abs(u1-u0)<10**(-7)).all()
		if stop==True:
			break
		else:
			u0=u1
	vv=2*np.arctan(np.sqrt((1+E)/(1-E))*np.tan(u1/2))
	#
	rr=A*(1-(E**2))/(np.ones_like(vv)+E*np.cos(vv))
	AA=np.cos(vv+W)
	BB=np.sin(vv+W)
	X=rr*BB*np.sin(I)
	Y=rr*(-AA*np.cos(WW)+BB*np.sin(WW)*np.cos(I))
	Z=rr*(-AA*np.sin(WW)-BB*np.cos(WW)*np.cos(I))
	return [X,Y,Z]



### THE MODEL FUNCTION

def model(a1,a2,a3,a4,p,P,A,E,I,W,WW,T0,tt):
	## projected distance
	pos	=	position(P,A,E,I,W,WW,T0,tt)
	fx	=	pos[0]
	fy	=	pos[1]
	fz	=	pos[2]
	z	=	np.sqrt(fy**2+fz**2)
	## cases
	case0 = np.where(	(fx <= 0)	| 	(z >= 1+p)								)
	case1 = np.where(	(fx >  0)	& 	(z == 0  )								)
	case2 = np.where(	(fx >  0)	& 	(z <  p  )								)
	case3 = np.where(	(fx >  0)	& 	(z == p  )								)
	case4 = np.where(	(fx >  0)	& 	(z >  p  )			&	(z < 1-p)		)
	case5 = np.where(	(fx >  0)	& 	(z == 1-p)								)
	case6 = np.where(	(fx >  0)	& 	(z >  1-p) 			&	(z**2-p**2<1) 	)
	case7 = np.where(	(fx >  0)	& 	(z**2-p**2 == 1)					 	)
	case8 = np.where(	(fx >  0)	& 	(z**2-p**2 >  1)	&	(z < 1+p) 		)
	## cross points
	zero = np.zeros(len(z))
	ones = np.ones_like(z)
	parr = ones*float(p)
	piar = ones*pi
	th = np.arcsin(p/z)
	ro = np.sqrt( z**2 - p**2 )
	ph = np.arccos( ( 1.0 - p**2 + z**2 )/( 2.0*z ) )
	## flux
	plusflux = np.zeros(len(z))
	plusflux[case1]	= INTCENT(a1,a2,a3,a4,p,z[case1]	,zero[case1]	,2*piar[case1]	)
	plusflux[case2]	= INTPLUS(a1,a2,a3,a4,p,z[case2]	,zero[case2]	,p+z[case2]		,piar[case2]	,p-z[case2]		)
	plusflux[case3]	= INTPLUS(a1,a2,a3,a4,p,z[case3]	,zero[case3]	,2*parr[case3]	,piar[case3]/2	,zero[case3]	)
	plusflux[case4]	= INTPLUS(a1,a2,a3,a4,p,z[case4]	,zero[case4]	,p+z[case4]		,th[case4]		,ro[case4]		)
	plusflux[case5]	= INTPLUS(a1,a2,a3,a4,p,z[case5]	,zero[case5]	,ones[case5]	,th[case5]		,ro[case5]		)
	plusflux[case6]	= INTPLUS(a1,a2,a3,a4,p,z[case6]	,ph[case6]		,ones[case6]	,th[case6]		,ro[case6]		)
	minsflux = np.zeros(len(z))
	minsflux[case4]	= INTMINS(a1,a2,a3,a4,p,z[case4]	,zero[case4]	,z[case4]-p		,th[case4]		,ro[case4]		)
	minsflux[case5]	= INTMINS(a1,a2,a3,a4,p,z[case5]	,zero[case5]	,z[case5]-p		,th[case5]		,ro[case5]		)
	minsflux[case6]	= INTMINS(a1,a2,a3,a4,p,z[case6]	,zero[case6]	,z[case6]-p		,th[case6]		,ro[case6]		)
	minsflux[case7]	= INTMINS(a1,a2,a3,a4,p,z[case7]	,zero[case7]	,z[case7]-p		,th[case7]		,ro[case7]		)
	minsflux[case8]	= INTMINS(a1,a2,a3,a4,p,z[case8]	,zero[case8]	,z[case8]-p		,ph[case8]		,ones[case8]	)
	starflux = np.zeros(len(z))
	starflux[case6]	= INTCENT(a1,a2,a3,a4,1,z[case6]	,zero[case6]	,ph[case6]		)
	starflux[case7]	= INTCENT(a1,a2,a3,a4,1,z[case7]	,zero[case7]	,ph[case7]		)
	starflux[case8]	= INTCENT(a1,a2,a3,a4,1,z[case8]	,zero[case8]	,ph[case8]		)
	F0 				= INTCENT(a1,a2,a3,a4,1,0			,0.0			,2.0*pi			)
	plt.plot(tt,plusflux)
	plt.plot(tt,minsflux)
	plt.plot(tt,starflux)
	plt.show()
	return np.array( 1 - 2.0*( plusflux + starflux - minsflux )/F0 )
### THE MODEL FUNCTION





#planet1
prex=0.15
pex=2.2
aex=9.0
eex=0.0
iex=87.5*pi/180
wex=0*pi/180
wwex=0*pi/180
t0ex=132.74052
#coeff
a1ex=0.6023
a2ex=-0.5110
a3ex=0.4655
a4ex=-0.1752
import time
#### test1
xx=np.arange(t0ex-pex/2,t0ex+pex/2,0.001)
t0=time.time()
yy=model(a1ex,a2ex,a3ex,a4ex,prex,pex,aex,eex,iex,wex,wwex,t0ex,xx)
print time.time()-t0
plt.plot(xx,yy,'k-',ms=2)
plt.xlabel(r'$time\,(t)\,[days]$')
plt.ylabel(r'$relative\,flux\,(f(t))$')
plt.ylim((plt.ylim()[0],1.001))
plt.xlim((xx[0],xx[-1]))
plt.show()
#### test1





