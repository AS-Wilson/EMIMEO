#####################################################################################################################
##### THIS IS NECESSARY, THIS GIVES US ALL THE FUNCTIONS WE NEED FOR NUMERICAL IMPLEMENTATION #####

import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import fft, fftshift, ifft, ifftshift
from mpl_toolkits.mplot3d import Axes3D


#####################################################################################################################
##### IGNORE THIS IF YOU WISH, THIS IS SIMPLY TO MAKE PRETTY GRAPHS #####

##### matplotlib graph settings #####
# Produce all graphs as PGFs so that they can be imported natively to latex and still be edited if needed
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




#####################################################################################################################
##### THE CODE TO IMPLEMENT THE SPLIT-STEP FOURIER METHOD #####

## The Fibre Material Properties
beta_2 = 0.05 #
gamma = 0 #


## Temporal Co-Ordinate Creation
# The higher the number of samples for a given signal duration, the better the curved gaussian shape will be after the
#   FFT:
no_of_samples = 1001
signal_duration = 80

# Create an array of time values with X number of equally-spaced values, a larger array size for a given signal
#   duration makes for a better frequency domain transform
time = np.linspace(-signal_duration/2, signal_duration/2, no_of_samples)

# This is the period between each discrete sample of the signal in the time domain
sample_period = signal_duration / (no_of_samples - 1)



## Spatial (z-axis) Co-Ordinate Creation
signal_propagation_distance = 10

# Create an array of z-axis values with X number of equally-spaced values, a larger array size for a given signal
#   duration makes for a better frequency domain transform
z_axis = np.linspace(-signal_propagation_distance/2, signal_propagation_distance/2, no_of_samples)

# This is the distance between each discrete sample of the signal on the z-axis
z_step_distance = signal_propagation_distance / (no_of_samples - 1)



## Frequency Domain Co-Ordinate Creation
# Angular frequency
# This is complicated but the angular frequency interval (distance between each measured discrete frequency of
#   the FFT) is given by this calculation
ang_freq_interval = (2.0 * np.pi) / signal_duration
# From the above we can then create a frequency axis with defined frequency points based on the number of
#   samples, this will be explained in the notes
ang_freq_axis = np.arange(-no_of_samples/2, no_of_samples/2) * ang_freq_interval 



## Input Envelope Definition
intensity = 1 # Signal intensity
init_waist = 1 # Used to define the Gaussian signal, this is the half-waist size at 1/e intensity

# The Gaussian (envelope) signal
envelope = (intensity * np.exp(-(((time - 5) ** 2) / (2 * (init_waist ** 2))))) + \
           (intensity * np.exp(-(((time + 5) ** 2) / (2 * (init_waist ** 2)))))
envelope = np.array(envelope)


## Input Spectrum Calculations
# First the raw intensity FFT transform, that will be arranged in order from -ve frequency to +ve frequency

# The FFT algorithm in python is missing a component (sample_period/sqrt(2*pi)) see:
#   https://cvarin.github.io/CSci-Survival-Guide/fft.html
envelope_raw_spectrum = fftshift(fft(envelope)) * sample_period/np.sqrt(2*np.pi)
envelope_raw_spectrum = np.array(envelope_raw_spectrum)

# Then one calculates the absolute intensity (|E(omega)|^2) for the input spectrum
envelope_abs_spectrum = abs(envelope_raw_spectrum)**2



## Split-Step Fourier Method
# Factors for calculating the evolution of the signal
dispersive_terms = [np.exp(complex(0, (((ang_freq_axis[0] ** 2) / 2) * beta_2 * z_step_distance)))]

for i in range(1, no_of_samples):
    dispersive_terms.append(np.exp(complex(0, (((ang_freq_axis[i] ** 2) / 2) * beta_2 * z_step_distance))))

# Since gamma is 0 in this code, this will be an array with no value, but the code is still useful
non_linear_term = gamma * z_step_distance



envelope_complex = [complex(envelope[0], 0)]

for i in range(1, no_of_samples):
    envelope_complex.append(complex(envelope[i], 0))

pulse_output_time = np.r_[envelope_complex]
pulse_output_spectrum = np.r_[envelope_raw_spectrum]


# Create vectors to manipulate when carrying out our steps
current_step_spectrum = envelope_raw_spectrum

# The actual Split-Step Fourier Method Calculations
for x in range(1, no_of_samples):

    for j in range(no_of_samples):
        # Dispersive (Linear) Step
        current_step_spectrum[j] = current_step_spectrum[j] * dispersive_terms[j]

    # Prep for Non-Linear Step
    current_step_time = ifft(ifftshift(current_step_spectrum)) * ang_freq_interval / np.sqrt(2.0*np.pi) * no_of_samples

    for j in range(no_of_samples):
        # Non-Linear Step
        current_step_time[j] = current_step_time[j] * np.exp(complex(0, (non_linear_term *
                                                                         (abs(current_step_time[j])**2))))

    pulse_output_time = np.row_stack((pulse_output_time, current_step_time))
    pulse_output_spectrum = np.row_stack((pulse_output_spectrum, current_step_spectrum))




#####################################################################################################################
########## GRAPHING ##########

## ENVELOPES IN TIME DOMAIN
plt.figure()

# Plot the input envelope vs time
plt.plot(time, envelope, color='blue', label='Initial Signal')

# Plot the output envelope vs time
plt.plot(time, abs(pulse_output_time[-1, :]), color='red', label='Final Signal')

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
plt.savefig('../Graphs/split-step-gvd-2d-time.pgf', bbox_inches = 'tight', pad_inches = 0)



## ENVELOPES IN FREQUENCY DOMAIN
plt.figure()

# Plot the input envelope intensity vs angular frequency
plt.plot(ang_freq_axis, envelope_abs_spectrum, color='green', label='Initial Spectrum')

# Plot the output envelope vs time
plt.plot(ang_freq_axis, abs(pulse_output_spectrum[-1, :])**2, color='orange', label='Final Spectrum')

# Titles and labels
plt.title(r'Input and Output Envelope - Frequency Domain')
plt.ylabel(r'$\tilde{F}(\omega)$')
plt.xlabel(r'Angular Frequency, $\omega$')
plt.legend(loc="upper left")


# Set limits
plt.xlim(-5, 5)

# Aspect Ratio
axes = plt.gca()
axes.set_aspect(35)

plt.tight_layout(pad=1.2) # Place everything slightly closer together

# Save graph in this location
plt.savefig('../Graphs/split-step-gvd-2d-freq.pgf', bbox_inches = 'tight', pad_inches = 0)




## 3D TIME DOMAIN GRAPH
# Create a mesh matrix to plot signal values against
T, Z = np.meshgrid(time, z_axis)

# Create figure and axis objects
fig = plt.figure()
ax = plt.axes(projection='3d')

# Plot the signal as it propagates
ax.plot_surface(T, Z, abs(pulse_output_time), rstride=1, cstride=1, cmap='viridis', edgecolor='none')

# Colouring the background of the graph
fig.patch.set_facecolor('white')
ax.set_facecolor('white')
ax.w_xaxis.set_pane_color((0.95, 0.95, 0.95, 0.95))
ax.w_yaxis.set_pane_color((0.95, 0.95, 0.95, 0.95))
ax.w_zaxis.set_pane_color((0.95, 0.95, 0.95, 0.95))

# Title and axis labels
ax.set_title(r'3D Plot of Signal Propagation in Time and Space')
ax.set_xlabel(r'Time, t')
ax.set_ylabel(r'Propagation distance, z')
ax.set_zlabel(r'$F(t)$')


# Aspect ratio settings
ax.auto_scale_xyz([-22, 22], [-5, 5], [0, 1.1])
ax.set_box_aspect((1.5, 2.0, 0.6))

plt.tight_layout(pad=0.8) # Place everything slightly closer together

# Save graph in this location
plt.savefig('../Graphs/split-step-gvd-3d-time.pgf', bbox_inches='tight', pad_inches=0)



## 3D FREQUENCY DOMAIN GRAPH
# Create a mesh matrix to plot signal values against
S, Zs = np.meshgrid(ang_freq_axis, z_axis)

# Create figure and axis objects
fig = plt.figure()
ax = plt.axes(projection='3d')

# Plot the signal as it propagates
ax.plot_surface(S, Zs, abs(pulse_output_spectrum)**2, rstride=1, cstride=1, cmap='viridis', edgecolor='none')

# Colouring the background of the graph
fig.patch.set_facecolor('white')
ax.set_facecolor('white')
ax.w_xaxis.set_pane_color((0.95, 0.95, 0.95, 0.95))
ax.w_yaxis.set_pane_color((0.95, 0.95, 0.95, 0.95))
ax.w_zaxis.set_pane_color((0.95, 0.95, 0.95, 0.95))

# Title and axis labels
ax.set_title(r'3D Plot of Signal Propagation in Freq. and Space')
ax.set_xlabel(r'Angular Frequency, $\omega$')
ax.set_ylabel(r'Propagation distance, z')
ax.set_zlabel(r'$\tilde{F}(\omega)$')

# Aspect ratio settings
ax.auto_scale_xyz([-5, 5], [-5, 5], [0, 4.1])
ax.set_box_aspect((1.5, 2.0, 0.6))

plt.tight_layout(pad=0.8) # Place everything slightly closer together

# Save graph in this location
plt.savefig('../Graphs/split-step-gvd-3d-freq.pgf', bbox_inches='tight', pad_inches=0)

# Display a window of each graph now that we're finished
plt.show()

