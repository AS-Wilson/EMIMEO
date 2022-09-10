#####################################################################################################################
##### THIS IS NECESSARY, THIS GIVES US ALL THE FUNCTIONS WE NEED FOR NUMERICAL IMPLEMENTATION #####

import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import fft, fftshift, ifft, ifftshift
from mpl_toolkits.mplot3d import Axes3D

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


def propagation_terms(ang_freq_axis, z_step_size, gamma, beta_2):
    """Something Here: TODO Later"""
    ## Factors used in the Split-Step Fourier Method
    # Since gamma is 0 in this code, this will be an array with no value, but the code is still useful
    non_linear_term = gamma * z_step_size

    # Produce all the dispersive terms first as they are frequency dependent
    for i in range(len(ang_freq_axis)):
        try:
            dispersive_terms = np.column_stack((dispersive_terms, np.array(np.exp(complex(0, ((ang_freq_axis[i] ** 2 / 2) * beta_2 * z_step_size))))))

        except NameError:
            dispersive_terms = np.array(np.exp(complex(0, ((ang_freq_axis[i] ** 2 / 2) * beta_2 * z_step_size))))

    return non_linear_term, dispersive_terms


def get_spectrum(signal, sample_period):
    ## Input Spectrum Calculations
    # First the raw intensity FFT transform, that will be arranged in order from -ve frequency to +ve frequency
    # The FFT algorithm in python is missing a component (sample_period/sqrt(2*pi)) see:
    #   https://cvarin.github.io/CSci-Survival-Guide/fft.html
    raw_spectrum = fftshift(fft(signal)) * sample_period / np.sqrt(2 * np.pi)

    # Then one calculates the absolute intensity (|E(omega)|^2) for the input spectrum
    abs_spectrum = abs(raw_spectrum) ** 2

    return raw_spectrum, abs_spectrum


def make_3d_graph(x, y, z, title=r'Graph', x_lab=r'X Axis', y_lab=r'Y Axis',
                  z_lab=r'Z Axis', auto_scale=None, aspect_ratio=None, fig_loc=None):
    """Something here: TODO Later"""
    ## 3D FREQUENCY DOMAIN GRAPH
    # Create a mesh matrix to plot signal values against
    X, Z = np.meshgrid(x, z)

    # Create figure and axis objects
    fig = plt.figure()
    ax = plt.axes(projection='3d')

    # Plot the signal as it propagates
    ax.plot_surface(X, Z, y, rstride=1, cstride=1, cmap='viridis', edgecolor='none')

    # Colouring the background of the graph
    fig.patch.set_facecolor('white')
    ax.set_facecolor('white')
    ax.w_xaxis.set_pane_color((0.95, 0.95, 0.95, 0.95))
    ax.w_yaxis.set_pane_color((0.95, 0.95, 0.95, 0.95))
    ax.w_zaxis.set_pane_color((0.95, 0.95, 0.95, 0.95))

    # Title and axis labels
    ax.set_title(title)
    ax.set_xlabel(x_lab)
    ax.set_ylabel(y_lab)
    ax.set_zlabel(z_lab)

    if auto_scale is None:
        pass

    else:
        # Aspect ratio settings
        ax.auto_scale_xyz(*auto_scale)


    if aspect_ratio is None:
        pass

    else:
        ax.set_box_aspect(aspect_ratio)

    plt.tight_layout(pad=0.8)  # Place everything slightly closer together


    if fig_loc is None:
        pass

    else:
        plt.savefig(fig_loc, format='pgf', bbox_inches='tight', pad_inches=0)




def split_step_fourier_method(input_signal, sample_period, no_of_samples, ang_freq_interval,
                              non_linear_term, dispersive_terms):
    """Something Here, TODO Later"""

    signal_propagation_spectrum, abs_signal_propagation_spectrum = get_spectrum(input_signal, sample_period)

    signal_propagation_time = abs(input_signal.astype(complex))

    # The actual Split-Step Fourier Method Calculations
    for i in range(0, no_of_samples-1): #(len(signal_propagation_time) != no_of_samples):

        try:
            signal_propagation_spectrum = np.vstack((signal_propagation_spectrum,
                                                        signal_propagation_spectrum[-1, :] * dispersive_terms))

            signal_propagation_time = np.vstack((signal_propagation_time,
                                                    ((ifft(ifftshift(signal_propagation_spectrum[-1, :])) *
                                                      ang_freq_interval / np.sqrt(2.0 * np.pi) * no_of_samples) *
                                                     np.exp((non_linear_term *
                                                             (abs(signal_propagation_spectrum[-1, :])**2)).astype(complex)))))

        except IndexError:
            signal_propagation_spectrum = np.vstack((signal_propagation_spectrum,
                                                     (signal_propagation_spectrum * dispersive_terms)))

            signal_propagation_time = np.vstack((signal_propagation_time,
                                                 ((ifft(ifftshift(signal_propagation_spectrum[-1, :])) *
                                                   ang_freq_interval / np.sqrt(2.0 * np.pi) * no_of_samples) *
                                                  np.exp((non_linear_term *
                                                          (abs(signal_propagation_spectrum[-1, :])**2)).astype(complex)))))


    return signal_propagation_time, signal_propagation_spectrum


def gaussian(time, t_0=None, FWHM=None, intensity=1, t_shift=0):
    """Public Method for producing a Gaussian Signal after providing various parameters."""

    try:
        t_0 = FWHM / (2 * np.sqrt(2 * np.log(2)))
        signal = intensity * np.exp(- ((time - t_shift) ** 2 / 2 * t_0 ** 2))
        return signal

    except TypeError:
        signal = intensity * np.exp(- ((time - t_shift) ** 2 / 2 * t_0 ** 2))
        return signal


#####################################################################################################################
##### THE CODE TO IMPLEMENT THE SPLIT-STEP FOURIER METHOD #####

## The Fibre Material Properties
beta_2 = 0.05  #
gamma = 0  #

## Temporal Co-Ordinate Creation
# The higher the number of samples for a given signal duration, the better the curved gaussian shape will
#   be after the FFT:
no_of_samples = 501
signal_duration = 80

# Create an array of time values with X number of equally-spaced values, a larger array size for a given
#   signal duration makes for a better frequency domain transform
time = np.linspace(-signal_duration / 2, signal_duration / 2, no_of_samples)

# This is the period between each discrete sample of the signal in the time domain
sample_period = signal_duration / (no_of_samples - 1)

## Spatial (z-axis) Co-Ordinate Creation
signal_propagation_distance = 10

# Create an array of z-axis values with X number of equally-spaced values, a larger array size for a given
#   signal duration makes for a better frequency domain transform
z_axis = np.linspace(-signal_propagation_distance / 2, signal_propagation_distance / 2, no_of_samples)

# This is the distance between each discrete sample of the signal on the z-axis
z_step_size = signal_propagation_distance / (no_of_samples - 1)

## Frequency Domain Co-Ordinate Creation
# Angular frequency this is complicated but the angular frequency interval (distance between each
#   measured discrete frequency of the FFT) is given by this calculation
ang_freq_interval = (2.0 * np.pi) / signal_duration
# From the above we can then create a frequency axis with defined frequency points based on the number of
#   samples, this will be explained in the notes
ang_freq_axis = np.arange(-no_of_samples / 2, no_of_samples / 2) * ang_freq_interval

## Input Envelope Definition
intensity = 1  # Signal intensity
init_waist = 1  # Used to define the Gaussian signal, this is the half-waist size at 1/e intensity
t_shift = 5

# The Gaussian (envelope) signal
input_signal = gaussian(time, t_0=init_waist, intensity=intensity, t_shift=t_shift) + \
               gaussian(time, t_0=init_waist, intensity=intensity, t_shift=-t_shift)

## Make all the parameters for the split-step Fourier Method from the fibre properties
non_linear_term, dispersive_terms = propagation_terms(ang_freq_axis, z_step_size, gamma, beta_2)

## Get spectrum of input signal for graphing
raw_input_spectrum, abs_input_spectrum = get_spectrum(input_signal, sample_period)


signal_propagation_time, signal_propagation_spectrum = split_step_fourier_method(input_signal,
                                                                                 sample_period,
                                                                                 no_of_samples,
                                                                                 ang_freq_interval,
                                                                                 non_linear_term, dispersive_terms)



#####################################################################################################################
########## GRAPHING ##########

## ENVELOPES IN TIME DOMAIN
plt.figure()

# Plot the input envelope vs time
plt.plot(time, input_signal, color='blue', label='Initial Signal')

# Plot the output envelope vs time
plt.plot(time, abs(signal_propagation_time[-1, :]), color='red', label='Final Signal')

# Titles and labels
plt.title(r'Input and Output Envelope - Time Domain')
plt.ylabel(r'$F(t)$')
plt.xlabel(r'Time, t')
plt.legend(loc="upper left")

# Set limits
plt.xlim(-20, 20)

# Aspect Ratio
axes = plt.gca()
axes.set_aspect(15)

plt.tight_layout(pad=1.2) # Place everything slightly closer together

# Save graph in this location
plt.savefig('../Graphs/split-step-gvd-2d-time-rewrite.pgf', format='pgf', bbox_inches='tight', pad_inches=0)


## ENVELOPES IN FREQUENCY DOMAIN
plt.figure()

# Plot the input envelope intensity vs angular frequency
plt.plot(ang_freq_axis, abs_input_spectrum, color='green', label='Initial Spectrum')

# Plot the output envelope vs time
plt.plot(ang_freq_axis, abs(signal_propagation_spectrum[-1, :]) ** 2, color='orange', label='Final Spectrum')

# Titles and labels
plt.title(r'Input and Output Envelope - Frequency Domain')
plt.ylabel(r'$\tilde{F}(\omega)$')
plt.xlabel(r'Angular Frequency, $\omega$')
plt.legend(loc="upper left")

# Set limits
plt.xlim(-5, 5)

# Aspect Ratio
axes = plt.gca()
axes.set_aspect(1)

plt.tight_layout(pad=1.2) # Place everything slightly closer together

# Save graph in this location
plt.savefig('../Graphs/split-step-gvd-2d-freq-rewrite.pgf', format='pgf', bbox_inches='tight', pad_inches=0)

index_0 = round((-20-time[0]) / sample_period)
index_1 = round((20-time[0]) / sample_period)

## 3D TIME DOMAIN GRAPH
make_3d_graph(time[index_0:index_1], abs(signal_propagation_time[:, index_0:index_1]), z_axis,
              title=r'3D Plot of Signal Propagation in Time and Space', x_lab=r'Time, t',
              y_lab=r'Propagation distance, z', z_lab=r'$F(t)$',
              auto_scale=([-22, 22], [-5, 5], [0, 1.1]),
              aspect_ratio=(1.5, 2.0, 0.6), fig_loc='../Graphs/split-step-gvd-3d-time-rewrite.pgf')

index_0 = round((-4-ang_freq_axis[0]) / ang_freq_interval)
index_1 = round((4-ang_freq_axis[0]) / ang_freq_interval)

## 3D FREQUENCY DOMAIN GRAPH
make_3d_graph(ang_freq_axis[index_0:index_1], abs(signal_propagation_spectrum[:, index_0:index_1]) ** 2, z_axis,
              title=r'3D Plot of Signal Propagation in Freq. and Space', x_lab=r'Angular Frequency, $\omega$',
              y_lab=r'Propagation distance, z', z_lab=r'$\tilde{F}(\omega)$',
              auto_scale=([-5, 5], [-5, 5], [0, 4.1]),
              aspect_ratio=(1.5, 2.0, 0.6), fig_loc='../Graphs/split-step-gvd-3d-freq-rewrite.pgf')


# Display a window of each graph now that we're finished
plt.show()
