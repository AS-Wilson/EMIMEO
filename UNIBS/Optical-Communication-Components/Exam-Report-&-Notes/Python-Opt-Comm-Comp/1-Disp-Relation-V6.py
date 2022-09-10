#####################################################################################################################
##### THIS IS NECESSARY, THIS GIVES US ALL THE FUNCTIONS WE NEED FOR NUMERICAL IMPLEMENTATION #####

import numpy as np
import cmath
import matplotlib.pyplot as plt
from OCCToolbox import dispersion_relation, find_intersections, make_2d_graph
import OCNToolbox

#####################################################################################################################
##### IGNORE THIS IF YOU WISH, THIS IS SIMPLY TO MAKE PRETTY GRAPHS #####

##### matplotlib graph settings #####
# Setting fonts for all the graphs etc.
plt.rcParams["text.usetex"] = True
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.serif"] = "ptserif"
plt.rcParams["font.size"] = 9

# Match the colouring if my Latex docs
COLOUR = '#2E2E2E'
plt.rcParams['text.color'] = COLOUR
plt.rcParams['axes.labelcolor'] = COLOUR
plt.rcParams['xtick.color'] = COLOUR
plt.rcParams['ytick.color'] = COLOUR
plt.rcParams['pgf.preamble'] = r"\usepackage[T1]{fontenc} \usepackage{mathpazo}"

# Set style to bmh and place the ticks for the axes on the outside to be readable
plt.style.use('bmh')
plt.rcParams['xtick.direction'] = 'out'
plt.rcParams['ytick.direction'] = 'out'




##### DISPERSION RELATION FOR V=6 #####
signal_duration = 10
sample_size = 500
X, LHS, RHS = dispersion_relation(6, signal_duration, LP_Mode=0, no_of_samples=sample_size)

LP_0x_locations = find_intersections(X, LHS, RHS)

##### GRAPHING #####
make_2d_graph(X, LHS, r'LHS', X, RHS, r'RHS, $V=6$', Legend=True,
              X_Label=r'$\frac{- X J_{-1}(X)}{J_{0}(X)}$', Y_Label=r'$\frac{Y K_{-1}(Y)}{K_{0}(Y)}$',
              Title=r'Dispersion Transcendence / Relationship for $l=0$', X_Lim=[0, 7], Y_Lim=[-10, 10])

plt.scatter(X[LP_0x_locations[0]], LHS[LP_0x_locations[0]], color='black', s=15, zorder=2.5)
plt.scatter(X[LP_0x_locations[1]], LHS[LP_0x_locations[1]], color='black', s=15, zorder=2.5)

plt.annotate(r'$LP_{01}$', xy=(X[LP_0x_locations[0]], RHS[LP_0x_locations[0]]),
             xycoords='data', xytext=(2.5, 3), textcoords='data', ha="center", va="center",
             arrowprops=dict(arrowstyle="->", color='black'))
plt.annotate(r'$LP_{11}$', xy=(X[LP_0x_locations[1]], RHS[LP_0x_locations[1]]),
             xycoords='data', xytext=(5.5, 6.25), textcoords='data', ha="center", va="center",
             arrowprops=dict(arrowstyle="->", color='black'))

plt.show()