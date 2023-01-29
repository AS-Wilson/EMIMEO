import numpy as np
from mpl_toolkits import mplot3d
import os as sys

import matplotlib.pyplot as plt
from cmath import pi, sqrt
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

    plt.savefig('../../../Graphs/cs-&-refrac-ind.pgf', bbox_inches = 'tight', pad_inches = 0) # Save graph in this location
    

def normalised_spatial_freq_calcs():
    n_clad = 1.453
    n_core = [1.456, 1.456, 1.458, 1.458, 1.46]
    radius = [2.5e-6, 3e-6, 2.5e-6, 3e-6, 2.5e-6]
    
    NA = []
    frac_ref_ind_diff = []
    normal_spac_freq = []
    
    lamb_1 = 800e-9
    
    for i in range(len(n_core)):
        NA.append( np.round_((np.sqrt((n_core[i]**2) - (n_clad**2))), decimals=5) )
        frac_ref_ind_diff.append( np.round_( (((n_core[i]**2)-(n_clad**2)) / (2*(n_core[i]**2))), decimals=5) )
        normal_spac_freq.append( np.round_( (((2 * np.pi * radius[i]) / (lamb_1)) * NA[i]), decimals=5) )
        print("Fibre,", i, "has a NA of", NA[i])
        print("Fibre,", i, "has a Delta of", frac_ref_ind_diff[i])
        print("Fibre,", i, "has a V of", normal_spac_freq[i])
        
        NA_req = 0.12
        lamb_2 = 750e-9
        v_req = 2.405
        
        a_min = ((np.round_(((v_req * lamb_2) / (2 * np.pi * NA_req)), decimals=8)) * (10e5))
        a_max = ((np.round_(((v_req * lamb_1) / (2 * np.pi * NA_req)), decimals=8)) * (10e5))
        
        d_min = np.round_((a_min * 2), decimals=2)
        d_max = np.round_((a_max * 2), decimals=2)
        print(a_min)
        print(a_max)
        print(d_min)
        print(d_max)
        

    
    
def complex_3d():
    
    fig_1 = plt.figure()
    ax = plt.axes(projection ='3d')
    
    ########## 3D Line ##########
    ########## REQUIRES: from mpl_toolkits import mplot3d
    fig_2 = plt.figure()
    
    # syntax for 3-D projection
    ax = plt.axes(projection ='3d')
     
    # defining all 3 axes
    z = np.linspace(0, 1, 100)
    x = z * np.sin(25 * z)
    y = z * np.cos(25 * z)
     
    # plotting
    ax.plot3D(x, y, z, 'green')
    ax.set_title('3D line plot geeks for geeks')
    plt.show()
    
    ########## 3D Scatter Graph ##########
    fig_3 = plt.figure()
 
    # syntax for 3-D projection
    ax = plt.axes(projection ='3d')
     
    # defining axes
    z = np.linspace(0, 1, 100)
    x = z * np.sin(25 * z)
    y = z * np.cos(25 * z)
    c = x + y
    ax.scatter(x, y, z, c = c)
     
    # syntax for plotting
    ax.set_title('3d Scatter plot geeks for geeks')
    plt.show()
    
    ########## Surface Graphs and Wireframes ##########
    # defining surface and axes
    x = np.outer(np.linspace(-2, 2, 10), np.ones(10))
    y = x.copy().T
    z = np.cos(x ** 2 + y ** 3)
     
    fig_4 = plt.figure()
     
    # syntax for 3-D plotting
    ax = plt.axes(projection ='3d')
     
    # syntax for plotting
    ax.plot_surface(x, y, z, cmap ='viridis', edgecolor ='green')
    ax.set_title('Surface plot geeks for geeks')
    plt.show()
    
    # function for z axea
    def f(x, y):
        return np.sin(np.sqrt(x ** 2 + y ** 2))
     
    # x and y axis
    x = np.linspace(-1, 5, 10)
    y = np.linspace(-1, 5, 10)
      
    X, Y = np.meshgrid(x, y)
    Z = f(X, Y)
     
    fig_5 = plt.figure()
    ax = plt.axes(projection ='3d')
    ax.plot_wireframe(X, Y, Z, color ='green')
    ax.set_title('wireframe geeks for geeks');


def main():
    step_index_graph()
    normalised_spatial_freq_calcs()


if __name__ == "__main__":
    main()
    
    