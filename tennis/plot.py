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

# density 689 
SWAN_VESTA_BOX=[0.016760206663199998, 29.37, 47.6, 17.4, 1.0] # mass, dimensions, thickness

def mi_rec_prism(m, h, d, w): # b, a, c
    return [
        m * (h**2 + w**2) / 12,
        m * (h**2 + d**2) / 12,
        m * (d**2 + w**2) / 12]

def mi_hollow_rec_prism(m, h, d, w, thickness):
    outer=mi_rec_prism(m, h, d, w)
    inner=mi_rec_prism(m, h - thickness, d - thickness, w - thickness)
    return [a - b for a, b in zip(outer, inner)]
        

def ellipsoid_projection(ax, A, B, n):
    '''Projects the ellipsoid defined by A\dot x = B onto 
    the plane normal to direction n0 into the plane n1, n2.'''

    width=2*math.sqrt(B/A[n[1] - 1])
    height=2*math.sqrt(B/A[n[2] - 1])

    print("Plotting ellipse {} x {}".format(width, height))

    ellipse = Ellipse(xy=(0,0), width=width, height=height, 
                        fc='None', lw=2)
    ax.add_patch(ellipse)

    return(width, height)

def ellipse(I=[1.0, 3.0, 5.0], T=100.0):

    plt.figure()

    ax=plt.gca()
    # energy ellipsoid from 1 onto 2, 3
    (width, height) = ellipsoid_projection(ax, I, T, [1, 2, 3])
    plt.axis('scaled')

    [I1, I2, I3]=I
    for LL in np.linspace(I1*T, I2*T, num=5):

        rhs=LL - I1*T
        width=2*math.sqrt(rhs/(I2**2 - I1*I2))
        height=2*math.sqrt(rhs/(I3**2 - I1*I3))
        ellipse = Ellipse(xy=(0,0), width=width, height=height, 
                        edgecolor="blue", fc='None', lw=2)
        ax.add_patch(ellipse)
#
#    #ax.set_xlim(0, 200)
    #ax.set_ylim(0, 200)

    #ax.relim()
    #ax.autoscale_view(tight=False)

    #plt.legend(range(0,10))
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

    matchbox=[0.1279, 0.1788, 0.2122]
    ellipse(matchbox, T=2*math.pi*0.1279) # realistic rotation
    
