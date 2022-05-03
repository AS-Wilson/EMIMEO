###################################################################################################
##### THIS IS NECESSARY, THIS GIVES US ALL THE FUNCTIONS WE NEED FOR NUMERICAL IMPLEMENTATION #####
from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt
import os as sys
import cmath
from scipy.fftpack import fft, fftshift, fftfreq, ifft, ifftshift


###################################################################################################
##### IGNORE THIS IF YOU WISH, THIS IS SIMPLY TO MAKE PRETTY GRAPHS #####
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




###################################################################################################
##### THE CODE TO IMPLEMENT THE SPLIT-STEP FOURIER METHOD #####


## The Fibre Material Properties
beta_2 = 1 #
gamma = 0 # 


## Temporal Co-Ordinate Creation
no_of_samples = 500 # The higher the number of samples for a given signal duration, the better the curved gaussian shape will be after the FFT
signal_duration = 40

# Create an array of time values with X number of equally-spaced values, a larger array size for a given signal duration makes for a better frequency domain transform
time = np.linspace(-signal_duration/2, signal_duration/2, no_of_samples)

sample_period = signal_duration / (no_of_samples - 1) # This is the period between each discrete sample of the signal in the time domain



## Spatial (z-axis) Co-Ordinate Creation
signal_propagation_distance = 10

# Create an array of z-axis values with X number of equally-spaced values, a larger array size for a given signal duration makes for a better frequency domain transform
z_axis = np.linspace(-signal_propagation_distance/2, signal_propagation_distance/2, no_of_samples)

z_step_distance = signal_propagation_distance / (no_of_samples - 1) # This is the distance between each discrete sample of the signal on the z-axis



## Frequency Domain Co-Ordinate Creation
# Angular frequency
# This is complicated but the angular frequency interval (distance between each measured discrete frequency of the FFT) is given by this calculation
ang_freq_interval = (2.0 * np.pi) / signal_duration
# From the above we can then create a frequency axis with defined frequency points based on the number of samples, this will be explained in the notes
ang_freq_axis = np.arange(-no_of_samples/2, no_of_samples/2) * ang_freq_interval 



## Input Envelope Definition
intensity = 1 # Signal intensity
init_waist = 1 # Used to define the Gaussian signal, this is the half-waist size at 1/e intensity

# The Gaussian (envelope) signal
envelope = intensity * np.exp( -( ((time)**2) / (2*(init_waist**2)) ))
envelope = np.array(envelope)


## Input Spectrum Calculations
# First the raw intensity FFT transform, that will be arranged in order from -ve frequency to +ve frequency
envelope_raw_spectrum = fftshift(fft(envelope)) * sample_period/np.sqrt(2*np.pi) # The FFT algorithm in python is missing a component (sample_period/sqrt(2*pi)) see: https://cvarin.github.io/CSci-Survival-Guide/fft.html
envelope_raw_spectrum = np.array(envelope_raw_spectrum)

# Then one calculates the absolute intensity (|E(omega)|^2) for the input spectrum
envelope_abs_spectrum = abs(envelope_raw_spectrum)**2



## Split-Step Fourier Method
# Initialise time and frequency matrices for output values
pulse_output_time = envelope
# pulse_output_time.append(envelope)

pulse_output_spectrum = envelope_raw_spectrum


# Factors for calculating the evolution of the signal
dispersive_factor = np.array([0+(((ang_freq_axis[0]**2)/2) * beta_2 * z_step_distance)j])
for i in range(1, no_of_samples):
    dispersive_factor_iterable = np.array([0+(((ang_freq_axis[i]**2)/2) * beta_2 * z_step_distance)j])
    np.vstack([dispersive_factor, dispersive_factor_iterable])   # Since beta_2 is 0 in this code, this will be an array with no value, but the code is still useful

non_linear_factor = complex(0, (gamma * z_step_distance)) # Since gamma is 0 in this code, this will be an array with no value, but the code is still useful


# # Factors for calculating the evolution of the signal
# if beta_2 != 0:
#     dispersive_factor = []
#     for i in range(len(ang_freq_axis)):
#         np.append(dispersive_factor, (complex(0, (((ang_freq_axis[i]**2)/2) * beta_2 * z_step_distance))))   # Since beta_2 is 0 in this code, this will be a array with no value, but the code is still useful
#
# else:
#     dispersive_factor = None
#
# if gamma != 0:
#     non_linear_factor = complex(0, (gamma * z_step_distance)) # Since gamma is 0 in this code, this will be a array with no value, but the code is still useful
#
# else:
#     non_linear_factor = None


# Create vectors to manipulate when carrying out our steps
current_step_time = envelope
current_step_spectrum = envelope_raw_spectrum


# The actual Split-Step Fourier Method Calculations
for i in range(1, no_of_samples):
    for j in range(1, no_of_samples):
        # Dispersive (Linear) Step
        current_step_spectrum[j] = current_step_spectrum[j] * cmath.exp(dispersive_factor[j])

    # Prep for Non-Linear Step
    current_step_time = ifft(ifftshift(current_step_spectrum)) * ang_freq_interval / np.sqrt(2.0*np.pi) * no_of_samples

    for j in range(1, no_of_samples):
        # Non-Linear Step
        current_step_time[j] = current_step_time[j] * cmath.exp(non_linear_factor * abs(current_step_time[j])**2 )

    pulse_output_time.append(current_step_time)
    pulse_output_spectrum.append(current_step_spectrum)
    




########## GRAPHING ##########

## INPUT ENVELOPE IN TIME DOMAIN
plt.figure()

# Plot the input envelope vs time
plt.plot(time, envelope, color='green')

# Titles and labels
plt.title(r'Input Envelope Intensity')
plt.ylabel(r'$F(t)$')
plt.xlabel(r'Time, t')

# Set limits
plt.xlim(-20, 20)

# Aspect Ratio
axes = plt.gca()
axes.set_aspect(16)

plt.tight_layout(pad=1.2) # Place everything slightly closer together

# plt.savefig('../Graphs/split-step-.pgf', bbox_inches = 'tight', pad_inches = 0) # Save graph in this location



## INPUT ENVELOPE SPECTRUM
plt.figure()

# Plot the input envelope intensity vs angular frequency
plt.plot(ang_freq_axis, envelope_abs_spectrum, color='green')

# Titles and labels
plt.title(r'Input Envelope Spectrum')
plt.ylabel(r'$\tilde{F}(\omega)$')
plt.xlabel(r'Angular Frequency, $\omega$')

# Set limits
plt.xlim(-50, 50)

# Aspect Ratio
axes = plt.gca()
axes.set_aspect(16)

plt.tight_layout(pad=1.2) # Place everything slightly closer together

# plt.savefig('../Graphs/split-step-.pgf', bbox_inches = 'tight', pad_inches = 0) # Save graph in this location




plt.figure()

ax = plt.axes(projection='3d')

ax.plot_surface(time, z_axis, pulse_output_time, cmap ='viridis', edgecolor ='green')

ax.set_title('SOMETHING')





plt.show()


#### PULSED GAUSSIAN NORMALISED SPECTRUM
##plt.figure()
##
### Plot the normalised angular frequency (omega/omega_0), we will do the calculation for the x-axis inside the argument of the plot
##plt.plot(ang_freq_axis/ang_freq, envelope_abs_spectrum, color='green')
##
### Titles and labels
##plt.title(r'Pulsed Gaussian Normalised Spectrum')
##plt.ylabel(r'$\tilde{E}(\omega)$')
##plt.xlabel(r'Normalised Angular Frequency, $\frac{\omega}{\omega_0}$')
##
### Plot Delta Omega on graph
##plt.annotate(r'\textbf{${\Delta\omega}\ll{\omega_0}$}', xy=(-1.07, 0.367), xycoords='data', xytext=(-0.75, 0.367), textcoords='data', ha="center", va="center", arrowprops=dict(arrowstyle="->", color='black'))
##plt.annotate(r'', xy=(-0.93, 0.367), xycoords='data', xytext=(-1, 0.367), textcoords='data', ha="center", va="center", arrowprops=dict(arrowstyle="->", color='black'))
##
### Set limits
##plt.xlim(-1.3, 1.3)
##
### Aspect Ratio
###axes = plt.gca()
###axes.set_aspect(0.85)
##
##plt.tight_layout(pad=1.2) # Place everything slightly closer together
##
### plt.savefig('../Graphs/split-step-.pgf', bbox_inches = 'tight', pad_inches = 0) # Save graph in this location


