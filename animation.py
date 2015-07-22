import matplotlib.pyplot as plt
import matplotlib.animation as pltanim
import numpy as np

import models


def animation(ldcoeffs, rprs, p, a, e, i, w, ww, t0, tt, save, file_name):
    light_curve = models.transit(ldcoeffs, rprs, p, a, e, i, w, ww, t0, tt)
    #
    pos1 = models.position(p, a, e, i * np.pi / 180, w * np.pi / 180, ww * np.pi / 180, t0, tt)
    fy1 = pos1[1]
    fz1 = pos1[2]
    fig = plt.figure()
    ax1 = fig.add_subplot(2, 1, 1)
    ax2 = fig.add_subplot(2, 1, 2)
    step = 1
    #
    window = len(tt) - 1
    limitx1 = -(tt[window] - tt[0]) / 2
    limitx2 = (tt[window] - tt[0]) / 2
    ax2.set_xlim((limitx1, limitx2))
    limity1 = min(light_curve) - 0.001
    limity2 = 1.005
    ax2.set_ylim((limity1, limity2))
    ax2.tick_params(axis='x', which='both', bottom='off', top='off', labelbottom='off')
    text10 = ax2.text(0, 1.001, 'TIME = ' + str(tt[0]), verticalalignment='baseline', horizontalalignment='center')
    line1, = ax2.plot([], [], 'r-', lw=4)
    line2, = ax2.plot([], [], 'ko', ms=3)
    #
    ax1.tick_params(axis='both', which='both',
                    bottom='off', left='off', top='off', right='off', labelbottom='off', labelleft='off')
    limitx1 = -2.05
    limitx2 = 2.05
    ax1.set_xlim((limitx1, limitx2))
    limity1 = -1.05
    limity2 = 1.05
    ax1.set_ylim((limity1, limity2))
    circle0 = plt.Circle((0, 0), 1, color=[0.9, 0.9, 0], fc=[0.9, 0.9, 0])
    ax1.add_patch(circle0)
    circle1 = plt.Circle((0, 0), rprs, color='k', fc='k')
    ax1.add_patch(circle1)
    #
    def animate(an):
        x = tt[max(0, step * an - window / 2):step * an] - tt[step * an]
        y = light_curve[max(0, step * an - window / 2):step * an]
        line1.set_data(x, y)
        line2.set_data(0, light_curve[step * an])
        text10.set_text('TIME = ' + str(tt[step * an]))
        circle1.center = (fy1[step * an], fz1[step * an])
        return circle0, circle1, line1, line2, text10,
    #
    if not save:
        anim = pltanim.FuncAnimation(fig, animate, frames=len(tt), interval=50, blit=False, repeat=True)
        plt.show()
    else:
        anim = pltanim.FuncAnimation(fig, animate, frames=len(tt), interval=50, blit=False, repeat=False)
        anim.save(file_name + '.mp4', fps=20, dpi=300, bitrate=-1, codec="libx264", extra_args=['-pix_fmt', 'yuv420p'])
        plt.show()