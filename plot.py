import argparse
import math

import matplotlib
import numpy as np
from matplotlib.patches import Ellipse
from matplotlib.pyplot import cm

(PLOT_NOPLOT, PLOT_PGF, PLOT_SCREEN, PLOT_PDF)=range(0, 4)

def pgf_setup(plot_status=PLOT_NOPLOT):
    if plot_status==PLOT_PGF:
        matplotlib.use('pgf')
        pgf_with_rc_fonts = {
            "font.family" : "serif",
            "font.serif" : [],
            "font.sans-serif" : ["DejaVu Sans"],
            "font.monospace" : [],
        }
    
        matplotlib.rcParams.update(pgf_with_rc_fonts)
    elif plot_status==PLOT_PDF:
        matplotlib.use('pdf')

def ellipse(I1=1.0, I2=3.0, I3=5.0, T=100.0, nellipses=10):


    # energy ellipsoid projected down

    width=2*math.sqrt(T/I2)
    height=2*math.sqrt(T/I3)

    plt.figure()
    ax = plt.gca()
    ellipse = Ellipse(xy=(0,0), width=width, height=height, 
                        fc='None', lw=2)
    ax.add_patch(ellipse)
    plt.plot([width/2,width/2], [-height/2, height/2], 'k-', lw=2)

    # intersection contours

    rhs_min=0           # LL = I1 * T
    rhs_max=T*(I3 - I1) # LL = I3 * T
    #rhs=LL - I1*T

    print("Energy eccentricity", math.sqrt(1 - I1/I3))
    print("Angular momentum eccentricity", math.sqrt(1 - (I2**2 - I1*I2)/(I3**2 - I1*I3)))
    print("Minimum angular momentum", I1*T)
    print("Maximum angular momentum", I3*T)

    color=iter(cm.rainbow(np.linspace(0, 1, nellipses)))
    for rhs in np.linspace(rhs_min, rhs_max, num=10):

        #print(rhs, LL, I1, I1*T)
        width=2*math.sqrt(rhs/(I2**2 - I1*I2))
        height=2*math.sqrt(rhs/(I3**2 - I1*I3))

        c=next(color)
        LL=rhs + I1*T
        print(LL, width, height)
        ellipse = Ellipse(xy=(0,0), width=width, height=height, 
                        edgecolor=c, fc='None', lw=2)
        ax.add_patch(ellipse)

    #ax.set_xlim(0, 200)
    #ax.set_ylim(0, 200)
    ax.relim()
    ax.autoscale_view(tight=False)

    plt.legend(range(0,10))
    plt.show()


if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Make plots")
    parser.add_argument('-s', '--screen', action="store_true")
    parser.add_argument('-p', '--pdf', action="store_true")
    args=parser.parse_args()
    if args.screen:
        plot_status=PLOT_SCREEN
    elif args.pdf:
        plot_status=PLOT_PDF
    else:
        plot_status=PLOT_PGF

    pgf_setup(plot_status)

    import matplotlib.pyplot as plt

    ellipse()
