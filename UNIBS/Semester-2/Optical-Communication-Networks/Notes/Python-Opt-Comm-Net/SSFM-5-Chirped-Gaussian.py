#####################################################################################################################
##### THIS IS NECESSARY, THIS GIVES US ALL THE FUNCTIONS WE NEED FOR NUMERICAL IMPLEMENTATION #####

import numpy as np
import matplotlib.pyplot as plt
from OCNToolbox import Split_Step_Fibre_Propagation, gaussian, make_3d_graph
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




def chirped_gaussian(time, chirp, t_0, intensity=1):
    """Public Method for producing a Chirped Gaussian Signal after providing various parameters."""

    signal = intensity * np.exp(- (1 + complex(0, 1*chirp) / 2 * (time ** 2 / t_0 ** 2)))
    # np.exp(-(1 + j * C) / 2 * (t.^ 2 / (tau0 ^ 2)));
    return signal


#####################################################################################################################
##### THE CODE TO IMPLEMENT THE SPLIT-STEP FOURIER METHOD #####

## The Fibre Material Properties
beta_2 = 0.05  #
gamma = 0  #

## Temporal Co-Ordinate Creation
# The higher the number of samples for a given signal duration, the better the curved gaussian shape will
#   be after the FFT:
no_of_samples = 1001
signal_duration = 80

# Create an array of time values with X number of equally-spaced values, a larger array size for a given
#   signal duration makes for a better frequency domain transform
time = np.linspace(-signal_duration / 2, signal_duration / 2, no_of_samples)

# This is the period between each discrete sample of the signal in the time domain
sample_period = signal_duration / (no_of_samples - 1)

## Spatial (z-axis) Co-Ordinate Creation
signal_propagation_distance = 10


## Input Envelope Definition
# Signal intensity
intensity = 1
# Initial Waist is used to define the Gaussian signal,
# this is the half-waist size at 1/e intensity
init_waist = 1
t_shift = 5
chirp = -2

# The Gaussian (envelope) signal
input_signal = gaussian(time, t_0=init_waist, intensity=intensity, t_shift=t_shift)
    # chirped_gaussian(time, chirp, init_waist, intensity=intensity)


##
Chirped_Gaussian = Split_Step_Fibre_Propagation(time, signal_propagation_distance,
                                                input_signal, no_of_z_samples=501,
                                                beta_2=1, gamma=None)


#####################################################################################################################
########## GRAPHING ##########

## ENVELOPES IN TIME DOMAIN
plt.figure()

# Plot the input envelope vs time
plt.plot(Chirped_Gaussian.time, input_signal, color='blue', label='Initial Signal')

# Plot the output envelope vs time
plt.plot(Chirped_Gaussian.time, abs(Chirped_Gaussian.signal_propagation_time[-1, :]), color='red', label='Final Signal')

# Titles and labels
plt.title(r'Input and Output Envelope - Time Domain')
plt.ylabel(r'$F(t)$')
plt.xlabel(r'Time, t')
plt.legend(loc="upper left")

# # Set limits
# plt.xlim(-20, 20)
#
# # Aspect Ratio
# axes = plt.gca()
# axes.set_aspect(15)
#
# plt.tight_layout(pad=1.2) # Place everything slightly closer together


## ENVELOPES IN FREQUENCY DOMAIN
plt.figure()

# Plot the input envelope intensity vs angular frequency
plt.plot(Chirped_Gaussian.ang_freq_axis, Chirped_Gaussian.abs_input_spec, color='green', label='Initial Spectrum')

# Plot the output envelope vs time
plt.plot(Chirped_Gaussian.ang_freq_axis, abs(Chirped_Gaussian.signal_propagation_spec[-1, :]) ** 2, color='orange', label='Final Spectrum')

# Titles and labels
plt.title(r'Input and Output Envelope - Frequency Domain')
plt.ylabel(r'$\tilde{F}(\omega)$')
plt.xlabel(r'Angular Frequency, $\omega$')
plt.legend(loc="upper left")

# # Set limits
# plt.xlim(-5, 5)
#
# # Aspect Ratio
# axes = plt.gca()
# axes.set_aspect(20)
#
# plt.tight_layout(pad=1.2) # Place everything slightly closer together


## 3D TIME DOMAIN GRAPH
make_3d_graph(Chirped_Gaussian.time, abs(Chirped_Gaussian.signal_propagation_time), Chirped_Gaussian.z_axis,
              title=r'3D Plot of Signal Propagation in Time and Space', x_lab=r'Time, t',
              y_lab=r'Propagation distance, z', z_lab=r'$F(t)$')

## 3D FREQUENCY DOMAIN GRAPH
make_3d_graph(Chirped_Gaussian.ang_freq_axis, abs(Chirped_Gaussian.signal_propagation_spec) ** 2, Chirped_Gaussian.z_axis,
              title=r'3D Plot of Signal Propagation in Freq. and Space', x_lab=r'Angular Frequency, $\omega$',
              y_lab=r'Propagation distance, z', z_lab=r'$\tilde{F}(\omega)$')

# Display a window of each graph now that we're finished
plt.show()
