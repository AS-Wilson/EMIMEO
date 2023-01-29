#####################################################################################################################
##### THIS IS NECESSARY, THIS GIVES US ALL THE FUNCTIONS WE NEED FOR NUMERICAL IMPLEMENTATION #####

import numpy as np
import matplotlib.pyplot as plt
import cmath

import sympy.core.numbers
from scipy.special import jv, kn, jve, kve
from sympy import Point, Line, Segment
from mpl_toolkits.mplot3d import Axes3D

#####################################################################################################################
##### matplotlib graph settings #####
# Produce all graphs as PGFs so that they can be imported natively to latex
# and still be edited if needed
plt.rcParams["pgf.texsystem"] = "pdflatex"
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



def find_intersections(X, Y1, Y2):
    """Function to find the intersection point (in terms of array index) of two data sets"""

    for idx in range(len(X)-1):

        if Y1[idx] is np.nan or Y1[idx+1] is np.nan or Y2[idx] is np.nan or Y2[idx+1] is np.nan:
            pass

        else:

            # Create two sympy Lines from points given
            # Line_A = Segment(Point(X[idx], Y1[idx]), Point(X[idx+1], Y1[idx+1]))
            # Line_B = Segment(Point(X[idx], Y2[idx]), Point(X[idx+1], Y2[idx+1]))

            # Using intersection() method find if there is an intersection in the two lines, this will return a
            # list of intersections, an empty list if there are no intersections
            Inter_Point = Segment(Point(X[idx], Y1[idx]), Point(X[idx+1], Y1[idx+1])).intersection(Segment(Point(X[idx], Y2[idx]), Point(X[idx+1], Y2[idx+1])))

            if len(Inter_Point) == 0:
                pass

            elif len(Inter_Point) == 1:
                intersect_value = Inter_Point[0]._eval_evalf()

                if abs(X[idx] - intersect_value[0]) < abs(X[idx+1] - intersect_value[0]):
                    Closest_Index = idx

                elif abs(X[idx+1] - intersect_value[0]) < abs(X[idx] - intersect_value[0]):
                    Closest_Index = idx + 1

                try:
                    intersections = np.hstack((intersections, np.array(Closest_Index)))

                except UnboundLocalError:
                    intersections = np.array(Closest_Index)


    if 'intersections' in globals():
        return intersections

    else:
        return None




def dispersion_relation(V, signal_duration, LP_Mode=0, no_of_samples=100):

    X = np.linspace(0, signal_duration, no_of_samples)

    for i in range(no_of_samples):
        try:
            LHS.append(-X[i] * jv(LP_Mode-1, X[i]) / jv(LP_Mode, X[i]))

            RHS.append((cmath.sqrt(V**2 - X[i]**2)).real * kn(LP_Mode-1, (cmath.sqrt(V**2 - X[i]**2)).real) / kn(LP_Mode, (cmath.sqrt(V**2 - X[i]**2)).real))

        except NameError:
            if not 'LHS' in locals():
                LHS = [-X[i] * jv(LP_Mode-1, X[i]) / jv(LP_Mode, X[i])]
                RHS = [(cmath.sqrt(V**2 - X[i]**2)).real * kn(LP_Mode-1, (cmath.sqrt(V**2 - X[i]**2)).real) / kn(LP_Mode, (cmath.sqrt(V**2 - X[i]**2)).real)]

            else:
                raise 'Unknown Error'

        if abs(LHS[i]) > 2*RHS[i]:
            LHS[i] = np.nan

        if cmath.isnan(RHS[i]):
            RHS[i] = 0


    return X, LHS, RHS




def fibre_dispersion(X, V, r_core, r_clad, LP_Location, no_of_samples):
    r_1 = np.linspace(0, r_core, no_of_samples)
    r_2 = np.linspace(r_core, r_clad, no_of_samples)

    k_t = (X[LP_Location]).real / r_core
    gamma = np.sqrt(V**2 - X[LP_Location]**2) / r_core

    amp_core = jv(0, (k_t*r_1))

    A_2 = amp_core[-1] / kn(0, (gamma*r_2[0]))

    amp_clad = A_2 * kn(0, (gamma*r_2))

    r_1 = np.array(r_1)
    r_2 = np.array(r_2)
    radius = np.concatenate((r_1, r_2.T))

    amp_core = np.array(amp_core)
    amp_clad = np.array(amp_clad)
    amplitude_fibre = np.concatenate((amp_core, amp_clad.T))

    return radius, amplitude_fibre




def make_2d_graph(*args, Legend=True, X_Label=r'X Axis', Y_Label=r'Y Axis', Title=r'Title', X_Lim=None, Y_Lim=None):
    plt.figure()

    for i in range(0, len(args), 3):
        plt.plot(args[i], args[i+1], label=args[(i+2)])

    if Legend is True:
        plt.legend()

    plt.xlim(X_Lim)
    plt.ylim(Y_Lim)
    plt.xlabel(X_Label)
    plt.ylabel(Y_Label)
    plt.title(Title)
    # plt.savefig('../../../Graphs/2d-charge-distribution-LP01-mode.pgf', bbox_inches='tight', pad_inches=0)




def get_effective_index(n_core, n_clad, r_core, X_Bar, V=None, Lambda=None, omega=None):
    """Something TODO Later"""

    mu_0 = 4e-7 * np.pi
    sig_0 = 8.854e-12
    c = 299792458

    if V is not None:
        omega = np.sqrt(V**2 / (mu_0 * sig_0 * ((n_core**2) - (n_clad**2)) * r_core**2))

    elif Lambda is not None:
        omega = (2 * np.pi * c) / Lambda


    beta_core = omega * n_core * np.sqrt(mu_0 * sig_0)
    beta_clad = omega * n_clad * np.sqrt(mu_0 * sig_0)

    beta_0 = beta_clad + ((beta_core - beta_clad) / 2)

    beta_l = (((n_core ** 2 * omega ** 2 * mu_0 * sig_0) - (X_Bar / r_core) ** 2) / (2 * beta_0)) + (beta_0 / 2)

    n_eff = beta_l / (omega * np.sqrt(mu_0 * sig_0))

    return n_eff




def main():
    ##### DISPERSION RELATION FOR V=6 #####
    signal_duration_1 = 10
    sample_size_1 = 500
    X, LHS, RHS = dispersion_relation(6, signal_duration_1, LP_Mode=0, no_of_samples=sample_size_1)

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




    ##### DISPERSION RELATIONSHIP AND X_BAR LOCATIONS FOR V=1, V=2 #####
    signal_duration_2 = 5
    sample_size_2 = 500

    V_1 = 1
    V_2 = 2

    X_1, LHS_1, RHS_1 = dispersion_relation(V_1, signal_duration_2, no_of_samples=sample_size_2)
    X_2, LHS_2, RHS_2 = dispersion_relation(V_2, signal_duration_2, no_of_samples=sample_size_2)

    V_1_LP_0x_locations = find_intersections(X_1, LHS_1, RHS_1)
    V_2_LP_0x_locations = find_intersections(X_2, LHS_2, RHS_2)

    ##### GRAPHING #####
    make_2d_graph(X_1, LHS_1, r'LHS', X_1, RHS_1, r'RHS, $V=1$', X_2, RHS_2, r'RHS, $V=2$',
                  Legend=True, X_Label=r'$\frac{- X J_{-1}(X)}{J_{0}(X)}$',
                  Y_Label=r'$\frac{Y K_{-1}(Y)}{K_{0}(Y)}$',
                  Title=r'Dispersion Transcendence / Relationship for $l=0$', X_Lim=[0, 2.5], Y_Lim=[-1, 3])

    plt.scatter(X_1[V_1_LP_0x_locations[0]], RHS_1[V_1_LP_0x_locations[0]], color='black', s=15, zorder=2.5)
    plt.scatter(X_2[V_2_LP_0x_locations[0]], RHS_2[V_2_LP_0x_locations[0]], color='black', s=15, zorder=2.5)

    plt.annotate(r'$LP_{01}$, V=1', xy=(X_1[V_1_LP_0x_locations[0]], RHS_1[V_1_LP_0x_locations[0]]),
                 xycoords='data', xytext=(1.5, 0.5), textcoords='data', ha="center", va="center",
                 arrowprops=dict(arrowstyle="->", color='black'))
    plt.annotate(r'$LP_{01}$, V=2', xy=(X_2[V_2_LP_0x_locations[0]], RHS_2[V_2_LP_0x_locations[0]]),
                 xycoords='data', xytext=(2, 2), textcoords='data', ha="center", va="center",
                 arrowprops=dict(arrowstyle="->", color='black'))

    plt.show()




    ##### TRANSVERSAL DISPERSION IN THE THE FIBRE CORE AND CLADDING #####
    r_core = 4.5e-6
    r_clad = 62.5e-6

    radius_V_1, amplitude_fibre_V_1 = fibre_dispersion(X_1, V_1, r_core, r_clad, V_1_LP_0x_locations[0], len(LHS_1))
    radius_V_2, amplitude_fibre_V_2 = fibre_dispersion(X_2, V_2, r_core, r_clad, V_2_LP_0x_locations[0], len(LHS_2))

    ##### GRAPHING #####
    make_2d_graph(radius_V_1, amplitude_fibre_V_1, r'$V=1$', radius_V_2, amplitude_fibre_V_2,
                  r'$V=2$', Legend=True, X_Label=r'Radius, r',
                  Y_Label=r'Field Intensity, $M_0(r)$',
                  Title=r'Field Dispersion In a Fibre, for $l=0$', X_Lim=[0, 20e-6], Y_Lim=[0, 1])

    plt.vlines(4.5e-6, 2, -2, color='gray', linestyles="dashed")

    plt.show()




    ##### BETA AND REFRACTIVE INDEX CALCULATIONS #####
    n_core = 1.447
    n_clad = 1.443

    effective_index_V1 = round(get_effective_index(n_core, n_clad, r_core, X_1[V_1_LP_0x_locations[0]], V=V_1), ndigits=3)
    print('Effective Index for V=1 and Core Radius of', r_core, 'calculated as: ', effective_index_V1)

    effective_index_V2 = round(get_effective_index(n_core, n_clad, r_core, X_1[V_2_LP_0x_locations[0]], V=V_2), ndigits=4)
    print('Effective Index for V=2 and Core Radius of', r_core, 'calculated as: ', effective_index_V2)


    Wavelengths = np.linspace(116e-9, 5e-6, 500)
    signal_duration = 3.5
    sample_size_3 = 1250

    # for i in range(len(Wavelengths)):
    #     V = ((2 * np.pi * r_core) / Wavelengths[i])**2 * (n_core**2 - n_clad**2)
    #     X, LHS, RHS = dispersion_relation(V, signal_duration, LP_Mode=0, no_of_samples=sample_size_3)
    #     X_Bar = find_intersections(X, LHS, RHS)
    #     effective_index = np.vstack((effective_index, get_effective_index(n_core, n_clad, r_core, X_Bar[0], V=V, Lambda=Wavelengths[i])))
    #
    #
    # make_2d_graph(Wavelengths, effective_index, r'$N_e$', Legend=True,
    #               X_Label=r'Wavelength', Y_Label=r'Effective Index',
    #               Title=r'Effective Index vs Wavelength')


if __name__ == "__main__":
    main()
    plt.show()


