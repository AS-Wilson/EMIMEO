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
        ##### 3D Graph of LP21 Mode #####
    bess_3d = plt.figure()
    ax_bess = bess_3d.gca(projection='3d')
        
    x = np.array([r*np.cos(angle) for r in radius])
    y = np.array([r*np.sin(angle) for r in radius])
    field_energy = np.array([field_value(2, 1, r, angle, 0.5) for r in radius])
    
    plt3d.Axes3D.plot_surface(ax_bess, x, y, field_energy, rstride=1, cstride=1, cmap=plt.cm.jet)
    # plt.legend()
    plt.xlabel(r'x axis, $\mu m$')
    plt.ylabel(r'y axis, $\mu m$')
    plt.title(r'3D Energy distribution in a fibre core in $LP_{21}$ mode')
    plt.savefig('../../../Graphs/3d-charge-distribution-LP21-mode.pgf', bbox_inches = 'tight', pad_inches = 0)
    
    
    ##### 3D Graph of LP02 Mode #####
    bess_3d = plt.figure()
    ax_bess = bess_3d.gca(projection='3d')
        
    x = np.array([r*np.cos(angle) for r in radius])
    y = np.array([r*np.sin(angle) for r in radius])
    field_energy = np.array([field_value(0, 2, r, angle, 0.5) for r in radius])
    
    plt3d.Axes3D.plot_surface(ax_bess, x, y, field_energy, rstride=1, cstride=1, cmap=plt.cm.jet)
    # plt.legend()
    plt.xlabel(r'x axis, $\mu m$')
    plt.ylabel(r'y axis, $\mu m$')
    plt.title(r'3D Energy distribution in a fibre core in $LP_{02}$ mode')
    plt.savefig('../../../Graphs/3d-charge-distribution-LP02-mode.pgf', bbox_inches = 'tight', pad_inches = 0)
    
    
    
    
    
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


def bess_funcs_and_fibre_core_graphs():
    ##### plotting multiple zeros of the first three bessel functions #####
    x = np.linspace(0, 12, 1000)
    plt.figure()
    for n in range(3):
        y = scsp.jn(n, x)
        plt.plot(x, y, label=r'$J_%d(x)$'%n)
    plt.legend()
    plt.xlabel(r'x')
    plt.ylabel(r'$J_m(x)$')
    plt.title(r'The First Three Bessel Functions $J_n(x)$')
    plt.savefig('../../../Graphs/bessel-functions.pgf', bbox_inches = 'tight', pad_inches = 0)
    
    ##### Plotting the J0 until the first zero in 2D #####
    plt.figure()
    j0_zero = scsp.jn_zeros(0, 1)
    x = np.linspace((-j0_zero), j0_zero, 1000)
    y = scsp.jn(0, x)
    
    plt.plot(x, y, label=r'$J_0(x)$')
    plt.legend()
    plt.xlabel(r'x axis, $\mu m$')
    plt.ylabel(r'Eletric Field Charge')
    plt.title(r'2D Energy distribution in a fibre core in $LP_{01}$ mode')
    plt.savefig('../../../Graphs/2d-charge-distribution-LP01-mode.pgf', bbox_inches = 'tight', pad_inches = 0)
    
    ##### plotting up until the first zero of a number of the besel functions #####
    #plt.figure()
    #for n in range(3):
    #    nth_zero = scsp.jn_zeros(n, 1)
    #    x = np.linspace((-nth_zero), nth_zero, 1000)
    #    y = scsp.jn(n, x)
    #    
    #    plt.plot(x, y, label=r'$J_%d(x)$'%n)
    #
    #plt.legend()
    #plt.xlabel(r'x axis, $\mu m$')
    #plt.ylabel(r'Eletric Field Charge')
    #plt.title(r'Energy distribution in a fibre core in $LP_{01}$ mode')
    #plt.savefig('../../../Graphs/charge-distribution.pgf', bbox_inches = 'tight', pad_inches = 0)
    
    ##### 3D Graph of LP01 Mode #####
    bess_3d = plt.figure()
    ax_bess = bess_3d.gca(projection='3d')
    angle = np.r_[0:2*np.pi:100j]
    radius = np.r_[0:5:100j]
    
    x = np.array([r*np.cos(angle) for r in radius])
    y = np.array([r*np.sin(angle) for r in radius])
    
    i = 0
    
    def field_energy_calc(x, zero):
        if x <= zero:
            return (np.cos(1)*np.cos(0*angle)*scsp.jn(0, x*zero))
        
        else:
            return 0
    
    
    field_energy = np.array([field_energy_calc(r, j0_zero) for r in radius])
    
    plt3d.Axes3D.plot_surface(ax_bess, x, y, field_energy, rstride=1, cstride=1, cmap=plt.cm.jet)
    # plt.legend()
    plt.xlabel(r'x axis, $\mu m$')
    plt.ylabel(r'y axis, $\mu m$')
    plt.title(r'3D Energy distribution in a fibre core in $LP_{01}$ mode')
    plt.savefig('../../../Graphs/3d-charge-distribution-LP01-mode.pgf', bbox_inches = 'tight', pad_inches = 0)    
    
    ##### Function for getting field intensities of any LP mode #####
    #def field_value(n, k, rad, ang, t):
    #    # n = order of the bessel function
    #    # k = argument of the bessel function
    #    # t = something*******
    #    # rad = radius to calculate for
    #    # ang = angles to calculate for
    #    nth_zero = scsp.jn_zeros(n, k)
    #    return np.cos(t)*np.cos(n*ang)*scsp.jn(n, rad*nth_zero)


    ##### 3D Graph of LP11 Mode #####
    #bess_3d = plt.figure()
    #ax_bess = bess_3d.gca(projection='3d')
    
    #field_energy = np.array([field_value(1, 1, r, angle, 0.5) for r in radius])
    
    #plt3d.Axes3D.plot_surface(ax_bess, x, y, field_energy, rstride=1, cstride=1, cmap=plt.cm.jet)
    # plt.legend()
    #plt.xlabel(r'x axis, $\mu m$')
    #plt.ylabel(r'y axis, $\mu m$')
    #plt.title(r'3D Energy distribution in a fibre core in $LP_{11}$ mode')
    #plt.savefig('../../../Graphs/3d-charge-distribution-LP11-mode.pgf', bbox_inches = 'tight', pad_inches = 0)
    
    



def main():
    # step_index_graph()
    # normalised_spatial_freq_calcs()
    bess_funcs_and_fibre_core_graphs()
    plt.show()


if __name__ == "__main__":
    main()
    
    