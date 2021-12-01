import numpy as np
import matplotlib.pyplot as plt
import scipy.special as scsp
from mpl_toolkits import mplot3d as plt3d
import os as sys
from sympy.functions.special.bessel import jn

# from matplotlib import rc

########################################
##### matplotlib graph settings #####
# Produce all graphs as PGFs so that they can be imported natively to latex 
# and still be edited if needed
plt.rcParams["pgf.texsystem"] = "pdflatex"
# Setting fonts for all the graphs etc.
plt.rcParams["text.usetex"] = True
plt.rcParams["font.family"] = "serif"
plt.rcParams["font.serif"] = "ptserif"
plt.rcParams["font.size"] = 11

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




########################################


def step_index_graph():
    
    T1 = 
    P_Sat = 
    eeta_p = 
    delta = 
    
    ########## Creating Data ##########
    angle = np.linspace( 0 , 2 * np.pi , 150 ) # Create an array of angle values
 
    radius_core = 5 # Set the radius of the core's circle
    radius_clad = 62.5 # Set the radius of the cladding's circle
     
    x_core = radius_core * np.cos(angle)
    y_core = radius_core * np.sin(angle)
    
    x_clad = radius_clad * np.cos(angle)
    y_clad = radius_clad * np.sin(angle)
    
    y_rad = [-62.5, -5, -5, 5, 5, 62.5]
    x_ind = [1.453, 1.453, 1.458, 1.458, 1.453, 1.453]
    
    ########## Graphing ##########
    
    fig_csec_ind, (ax_c_sec, ax_ind) = plt.subplots(1, 2, sharey='row')
    
    ## Left hand plot
    ax_c_sec.plot(x_core, y_core, color='black')
    ax_c_sec.plot(x_clad, y_clad, color='black')
    # Set limits and aspect ratio
    ax_c_sec.set_xlim(-65, 65)
    ax_c_sec.set_ylim(-65, 65)
    ax_c_sec.set_aspect(1)
    # Titles and labels
    ax_c_sec.set_title(r'Fiber Cross-Section')
    ax_c_sec.set_ylabel(r'x Dimension ($\mu$m)')
    ax_c_sec.set_xlabel(r'y Dimension ($\mu$m)')
    
    ## Right hand plot
    ax_ind.plot(x_ind, y_rad)
    ax_ind.set_xlim(1.45, 1.46) # Set limits
    # Set the aspect ratio divide the range of x values in the smaller graph by 
    # the range of y or y values in the larger graph.
    # For every x value there are this many y values, 
    # i.e. the y axis increments this number of times for every x value 
    # change in the x axis
    ax_ind.set_aspect(0.00007692307692)
    # Hide top and 
    
    # Position the y-axis label of the RH graph to the right to be more readable
    ax_ind.yaxis.tick_right()
    ax_ind.yaxis.set_label_position("right")
    # Titles and labels
    ax_ind.set_title(r'Refractive Index vs Radius')
    ax_ind.set_ylabel(r'Radius, r ($\mu$m)')
    ax_ind.set_xlabel(r'Refractive Index, $\lambda$')
    
    # plt.rc('text', usetex=True)
    # plt.rc('font', family='paratype')
    # plt.suptitle(r'Fibre Cross-Section and Corresponding Refractive Index', 
                 # size=16)
    
    plt.tight_layout(pad=0.001) # Place everything slightly closer together
    # plt.show() # Print the graph
    
    


def main():
    step_index_graph()


if __name__ == "__main__":
    main()
    
    