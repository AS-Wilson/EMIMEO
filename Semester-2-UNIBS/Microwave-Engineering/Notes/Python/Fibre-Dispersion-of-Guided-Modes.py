#####################################################################################################################
##### THIS IS NECESSARY, THIS GIVES US ALL THE FUNCTIONS WE NEED FOR NUMERICAL IMPLEMENTATION #####

import numpy as np
import matplotlib.pyplot as plt
import cmath
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



def find_intersections(X, Y1, Y2, accuracy=1e-3):
    """Function to find the intersection point (in terms of array index) of two data sets"""

    def calculate_any_intersection(A1x, A1y, A2x, A2y, B1x, B1y, B2x, B2y):

        # Create four sympy Points, then 2 symoy Lines from those points
        PA1, PA2, PB1, PB2 = Point(A1x, A1y), Point(A2x, A2y), Point(B1x, B1y), Point(B2x, B2y)
        Line_A = Segment(PA1, PA2)
        Line_B = Segment(PB1, PB2)

        # Using intersection() method find if there is an intersection in the two lines, this will return a
        # list of intersections, an empty list if there are no intersections
        Intersection_Point = Line_A.intersection(Line_B)
        return Intersection_Point

    def find_nearest_XY_index_and_store(x1, x2, inter_x_point, i):

        if abs(x1 - inter_x_point) < abs(x2 - inter_x_point):
            Closest_Index = i

        elif abs(x2 - inter_x_point) < abs(x1 - inter_x_point):
            Closest_Index = i+1

        return Closest_Index




    for i in range(len(X)-1):
        Inter_Point = calculate_any_intersection(X[i], Y1[i], X[i+1], Y1[i+1], X[i], Y2[i], X[i+1], Y2[i+1])
        if len(Inter_Point) == 0:
            pass

        elif len(Inter_Point) == 1:
            Intersec_Value  = Inter_Point[0]._eval_evalf()
            find_nearest_XY_index_and_store(X[i], X[i+1], Intersec_Value.x(), i)

            try:
                if Index_to_Store != intersections[-1]:
                    np.concatenate((intersections, np.array(Index_to_Store).T))

                else:
                    pass

            except IndexError:
                intersections = np.array(Index_to_Store)

        else:
            raise ValueError("It should not be possible to have more than one intersection considering we "
                             "should only have two straight lines with four points")





def dispersion_relation(V, signal_duration, l=0, no_of_samples=10000):

    X = np.linspace(0, signal_duration, no_of_samples)

    for i in range(no_of_samples):
        try:
            LHS.append(-X[i] * jv(l-1, X[i]) / jv(l, X[i]))

            RHS.append((cmath.sqrt(V**2 - X[i]**2)).real * kn(l-1, (cmath.sqrt(V**2 - X[i]**2)).real) / kn(l, (cmath.sqrt(V**2 - X[i]**2)).real))

        except NameError:
            if not 'LHS' in locals():
                LHS = [-X[i] * jv(l-1, X[i]) / jv(l, X[i])]
                RHS = [(cmath.sqrt(V**2 - X[i]**2)).real * kn(l-1, (cmath.sqrt(V**2 - X[i]**2)).real) / kn(l, (cmath.sqrt(V**2 - X[i]**2)).real)]

            else:
                raise 'Unknown Error'

        if abs(LHS[i]) > 30:
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

    if not V is None:
        omega = np.sqrt(V**2 / (mu_0 * sig_0 * ((n_core**2) - (n_clad**2)) * r_core**2))

    elif not Lambda is None:
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
    X, LHS, RHS = dispersion_relation(6, signal_duration_1, l=0)

    LP_0x_locations = find_intersections(LHS, RHS, accuracy=0.006)

    ##### GRAPHING #####
    make_2d_graph(X, LHS, r'LHS', X, RHS, r'RHS, $V=6$', Legend=True,
                  X_Label=r'$\frac{- X J_{-1}(X)}{J_{0}(X)}$', Y_Label=r'$\frac{Y K_{-1}(Y)}{K_{0}(Y)}$',
                  Title=r'Dispersion Transcendence / Relationship for $l=0$', X_Lim=[0, 7], Y_Lim=[-10, 10])

    plt.scatter(X[LP_0x_locations[0]], LHS[LP_0x_locations[0]], color='black', s=15, zorder=2.5)
    plt.scatter(X[LP_0x_locations[1]], RHS[LP_0x_locations[1]], color='black', s=15, zorder=2.5)

    plt.annotate(r'$LP_{01}$', xy=(X[LP_0x_locations[0]], RHS[LP_0x_locations[0]]),
                 xycoords='data', xytext=(2.5, 3), textcoords='data', ha="center", va="center",
                 arrowprops=dict(arrowstyle="->", color='black'))
    plt.annotate(r'$LP_{11}$', xy=(X[LP_0x_locations[1]], RHS[LP_0x_locations[1]]),
                 xycoords='data', xytext=(5.5, 6.25), textcoords='data', ha="center", va="center",
                 arrowprops=dict(arrowstyle="->", color='black'))



    ##### DISPERSION RELATIONSHIP AND X_BAR LOCATIONS FOR V=1, V=2 #####
    signal_duration_2 = 5
    V_1 = 1
    V_2 = 2

    X_1, LHS_1, RHS_1 = dispersion_relation(V_1, signal_duration_2)
    X_2, LHS_2, RHS_2 = dispersion_relation(V_2, signal_duration_2)

    V_1_LP_0x_locations = find_intersections(LHS_1, RHS_1, accuracy=0.002)
    V_2_LP_0x_locations = find_intersections(LHS_2, RHS_2, accuracy=0.001)

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



    ##### BETA AND REFRACTIVE INDEX CALCULATIONS #####
    n_core = 1.447
    n_clad = 1.443

    effective_index_V1 = round(get_effective_index(n_core, n_clad, r_core, X_1[V_1_LP_0x_locations[0]], V=V_1), ndigits=3)
    print('Effective Index for V=1 and Core Radius of', r_core, 'calculated as: ', effective_index_V1)

    effective_index_V2 = round(get_effective_index(n_core, n_clad, r_core, X_1[V_2_LP_0x_locations[0]], V=V_2), ndigits=4)
    print('Effective Index for V=2 and Core Radius of', r_core, 'calculated as: ', effective_index_V2)


    Wavelengths = np.linspace(5e-9, 5e-6, 10000)
    signal_duration = 20

    try:
        for i in range(1, 10000):
            V = ((2 * np.pi * r_core) / Wavelengths[i])**2 * (n_core**2 - n_clad**2)
            X, LHS, RHS = dispersion_relation(V, signal_duration, l=0, no_of_samples=10000)
            X_Bar = find_intersections(LHS, RHS, accuracy=0.001)
            effective_index =  np.concatenate((effective_index, get_effective_index(n_core, n_clad, r_core, X_Bar[0], V=V)))

    except UnboundLocalError:
        V = ((2 * np.pi * r_core) / Wavelengths[0]) ** 2 * (n_core ** 2 - n_clad ** 2)
        X, LHS, RHS = dispersion_relation(V, signal_duration, l=0, no_of_samples=10000)
        X_Bar = find_intersections(LHS, RHS, accuracy=0.00001)
        effective_index = np.matrix(get_effective_index(n_core, n_clad, r_core, X_Bar[0], V=V))

    make_2d_graph(Wavelengths, effective_index, r'N_e', Legend=True,
                  X_Label=r'Wavelength', Y_Label=r'Effective Index',
                  Title=r'Effective Index vs Wavelength')


if __name__ == "__main__":
    main()
    plt.show()


