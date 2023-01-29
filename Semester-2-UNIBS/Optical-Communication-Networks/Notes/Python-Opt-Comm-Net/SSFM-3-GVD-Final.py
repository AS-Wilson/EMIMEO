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


#####################################################################################################################
##### THE CODE TO IMPLEMENT THE SPLIT-STEP FOURIER METHOD #####

## Temporal Co-Ordinate Creation
# The higher the number of samples for a given signal duration, the better the curved gaussian shape will
#   be after the FFT:
no_of_samples = 1501
signal_duration = 100
# Create an array of time values with X number of equally-spaced values, a larger array size for a given
#   signal duration makes for a better frequency domain transform
time = np.linspace(-signal_duration / 2, signal_duration / 2, no_of_samples)
# This is the period between each discrete sample of the signal in the time domain
sample_period = signal_duration / (no_of_samples - 1)

## Spatial (z-axis) Co-Ordinate Creation
signal_propagation_distance = 10

## Input Envelope Definition
intensity = 1  # Signal intensity
init_waist = 1  # Used to define the Gaussian signal, this is the half-waist size at 1/e intensity
t_shift = 5

# The Gaussian (envelope) signal
input_signal = gaussian(time, t_0=init_waist, intensity=intensity, t_shift=t_shift) + \
               gaussian(time, t_0=init_waist, intensity=intensity, t_shift=-t_shift)


GVD = Split_Step_Fibre_Propagation(time, signal_propagation_distance, input_signal, no_of_z_samples=101, beta_2=0.05)


#####################################################################################################################
########## GRAPHING ##########

## ENVELOPES IN TIME DOMAIN
plt.figure()

# Plot the input envelope vs time
plt.plot(GVD.time, input_signal, color='blue', label='Initial Signal')

# Plot the output envelope vs time
plt.plot(GVD.time, abs(GVD.signal_propagation_time[-1, :]), color='red', label='Final Signal')

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
plt.savefig('../Graphs/split-step-2d-GVD-final-time.pgf', bbox_inches='tight', pad_inches=0)




## ENVELOPES IN FREQUENCY DOMAIN
plt.figure()

# Plot the input envelope intensity vs angular frequency
plt.plot(GVD.ang_freq_axis, GVD.abs_input_spec, color='green',
         label='Initial Spectrum')

# Plot the output envelope vs time
plt.plot(GVD.ang_freq_axis, abs(GVD.signal_propagation_spec[-1, :]) ** 2, color='orange', label='Final Spectrum')

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
plt.savefig('../Graphs/split-step-2d-GVD-final-freq.pgf', bbox_inches='tight', pad_inches=0)




start_index = round((-15 - GVD.time[0]) / sample_period)
end_index = round((15 - GVD.time[0]) / sample_period)

## 3D TIME DOMAIN GRAPH
make_3d_graph(GVD.time[start_index:end_index],
              abs(GVD.signal_propagation_time[:, start_index:end_index]),
              GVD.z_axis, title=r'3D Plot of Signal Propagation in Time and Space', x_lab=r'Time, t',
              y_lab=r'Propagation distance, z', z_lab=r'$F(t)$',
              auto_scale=([-15, 15], [0, 10], [0, 1.1]),
              aspect_ratio=(1.5, 2.0, 0.6), fig_loc='../Graphs/split-step-3d-GVD-final-time.pgf')




start_index = round((-5 - GVD.ang_freq_axis[0]) / GVD.ang_freq_interval)
end_index = round((5 - GVD.ang_freq_axis[0]) / GVD.ang_freq_interval)

## 3D FREQUENCY DOMAIN GRAPH
make_3d_graph(GVD.ang_freq_axis[start_index:end_index],
              abs(GVD.signal_propagation_spec[:, start_index:end_index]) ** 2, GVD.z_axis,
              title=r'3D Plot of Signal Propagation in Freq. and Space', x_lab=r'Angular Frequency, $\omega$',
              y_lab=r'Propagation distance, z', z_lab=r'$\tilde{F}(\omega)$',
              auto_scale=([-5, 5], [0, 10], [0, 4.1]),
              aspect_ratio=(1.5, 2.0, 0.6), fig_loc='../Graphs/split-step-3d-GVD-final-freq.pgf')

# Display a window of each graph now that we're finished
plt.show()
