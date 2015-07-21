import matplotlib.pyplot as plt
import matplotlib.animation as pltanim
import matplotlib.cm as cm
import numpy as np

import models
import tasks


def animation(ldcoeffs, RpRs, P, a, e, i, W, WW, T0, tt, save, file_name):
	a1, a2, a3, a4 = ldcoeffs
	light_curve = models.transit(ldcoeffs, RpRs, P, a, e, i, W, WW, T0, tt)
	x=np.arange(-2.5,2.5,0.01)
	y=np.arange(-1.5,1.5,0.01)
	xx,yy=np.meshgrid(x,y)
	star = ( 1.0 - a1*(1.0 - (1.0 - (xx**2+yy**2))**(1.0/4)) - a2*(1.0 - (1.0 - (xx**2+yy**2))**(1.0/2)) - a3*(1.0 - (1.0 - (xx**2+yy**2))**(3.0/4)) - a4*(1.0 - (1.0 - (xx**2+yy**2))) )
	star[np.where(np.isnan(star))]=0
	###########
	pos1=tasks.position(P,a,e,i*np.pi/180,W*np.pi/180,WW*np.pi/180,T0,tt)
	fx1=pos1[0]
	fy1=pos1[1]
	fz1=pos1[2]	
	fig = plt.figure()
	ax1 = fig.add_subplot(2,1,1)
	ax2 = fig.add_subplot(2,1,2)
	step=1
	##############
	window=len(tt)-1
	limitx1=-(tt[window]-tt[0])/2
	limitx2=(tt[window]-tt[0])/2
	ax2.set_xlim((limitx1,limitx2)) 
	limity1=min(light_curve)-0.001
	limity2=1.005
	ax2.set_ylim((limity1,limity2))
	ax2.tick_params(axis='x', which='both', bottom='off',top='off', labelbottom='off')
	text10 = ax2.text(0,1.001,'TIME = '+str(tt[0]),verticalalignment='baseline', horizontalalignment='center')
	line1, = ax2.plot([], [], 'r-', lw=4)
	line2, = ax2.plot([], [], 'ko', ms=3)
	###############
	ax1.tick_params(axis='both', which='both', bottom='off',left='off',top='off',right='off', labelbottom='off',labelleft='off')
	limitx1=-2.05
	limitx2=2.05
	ax1.set_xlim((limitx1,limitx2))
	limity1=-1.05
	limity2=1.05
	ax1.set_ylim((limity1,limity2))
	circle0=plt.Circle((0,0), 1, color=[0.9,0.9,0], fc=[0.9,0.9,0])
	ax1.add_patch(circle0)
	circle1=plt.Circle((0,0), RpRs, color='k', fc='k')
	ax1.add_patch(circle1)
	###############
	def animate(i):
		x=tt[max(0,step*i-window/2):step*i]-tt[step*i]
		y=light_curve[max(0,step*i-window/2):step*i]
		line1.set_data(x, y)
		line2.set_data(0, light_curve[step*i])
		text10.set_text('TIME = '+str(tt[step*i]))
		circle1.center=(fy1[step*i],fz1[step*i])
		return circle0,circle1,line1,line2,text10,
	##############
	if not save:
		anim = pltanim.FuncAnimation(fig, animate, frames=len(tt), interval=50, blit=False, repeat=True)
		plt.show()
	else:
		anim = pltanim.FuncAnimation(fig, animate, frames=len(tt), interval=50, blit=False, repeat=False)
		anim.save(file_name + '.mp4', fps=20, dpi=300, bitrate=-1, codec="libx264", extra_args=['-pix_fmt', 'yuv420p'])
		plt.show()
	##############

