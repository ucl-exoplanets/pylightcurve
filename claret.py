import glob
import numpy as np
import glob
import os
__location__ = os.path.realpath(os.path.join(os.getcwd(), os.path.dirname(__file__)))
##########################################################################################

def ldcoeff(Mett,Teff,Logg):
	filtertext			=	'Choose one of the available filters: \n u v b y U B V R I J H K : '
	filterlist			=	[	['u','v','b','y','U','B','V','R','I','J','H','K'],
								[4  ,5  ,6  ,7  ,8  ,9  ,10 ,11 ,12 ,13 ,14 ,15 ]
							]
	Filter				=	filterlist[1][filterlist[0].index(raw_input(filtertext))]
	tables, mett		=	np.loadtxt(glob.glob(__location__+'/*claretinfo*')[0],usecols=(0,4),unpack=True)
	Table				=	str(int(tables[np.argmin(abs(Mett-mett))]))
	File				=	glob.glob(__location__+'/*/TABLE'+Table)[0]
	logglist, tefflist	=	np.loadtxt(File,usecols=(1,2),unpack=True,skiprows=5)
	Teff				=	tefflist[np.argmin(abs(Teff-tefflist))]
	Logg				=	logglist[np.argmin(abs(Logg-logglist))]
	for i in open(File).readlines()[5:]:
		coef = i.split()[0]
		numb = float(i.split()[Filter])
		logg = float(i.split()[1])
		teff = float(i.split()[2])
		if ( logg == Logg and teff == Teff ):
			print '{0}: {1}'.format(coef, numb)
	