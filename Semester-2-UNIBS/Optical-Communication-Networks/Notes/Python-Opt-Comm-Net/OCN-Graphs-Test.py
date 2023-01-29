#####################################################################################################################
##### THIS IS NECESSARY, THIS GIVES US ALL THE FUNCTIONS WE NEED FOR NUMERICAL IMPLEMENTATION #####

import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import fft, fftshift, ifft, ifftshift
from mpl_toolkits.mplot3d import Axes3D

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


######################################

def gaussian_demo():
    ########## PLOT DATA CREATION ##########
    no_of_samples = 5000000  # This number of samples will give a nice curved gaussian shape in the FFT

    intensity = 2  # Signal intensity
    init_waist = 1  # Used to define the Gaussian signal, this is the half-waist size at 1/e intensity
    freq = 4  # Carrier frequency
    ang_freq = 2 * np.pi * freq  # Carrier angular frequency

    signal_duration = 50

    time = np.linspace(-signal_duration / 2, signal_duration / 2,
                       no_of_samples)  # Create an array of time values with them being 50000 equally-spaced values, the large array size makes for a better frequency domain tranform

    sample_period = signal_duration / no_of_samples  # This is the period between each discrete sample of the signal

    carrier = intensity * np.sin(ang_freq * (time))  # The carrier (cosine) portion of the Pulsed Gaussian
    envelope = intensity * np.exp(-(((time - 3.25) ** 2) / (2 * (
                init_waist ** 2))))  # The Gaussian (envelope) signal itself, this is shifted to the right by the -3.25 term as a real signal might be (it also makes the graph pretty)

    gaussian_pulse = intensity * np.exp(-(((time - 3.25) ** 2) / (2 * (init_waist ** 2)))) * np.sin(ang_freq * (time))

    ## SPECTRUM CALCULATIONS
    # Angular frequency
    # This is complicated but the angular frequency interval (distance between each measured discrete frequency of the FFT) is given by this calculation
    ang_freq_interval = (2.0 * np.pi) / signal_duration
    # From the above we can then create a frequency axis with discrete frequency points based on the number of samples, this will be explained in the notes
    ang_freq_axis = np.arange(-no_of_samples / 2, no_of_samples / 2) * ang_freq_interval

    # Spectrum
    gaussian_pulse_spectrum = fftshift(fft(gaussian_pulse)) * sample_period / np.sqrt(
        2 * np.pi)  # The FFT algorithm in python is missing a component (sample_period/sqrt(2*pi)) see: https://cvarin.github.io/CSci-Survival-Guide/fft.html
    gaussian_pulse_spec = abs(gaussian_pulse_spectrum) ** 2  # This calculates the intensity (|E(omega)|^2) for the spectrum graph

    ########## GRAPHING ##########

    ## DECONSTRUCTED GAUSSIAN PLOT
    fig_gauss_decon, (sub_carrier, sub_envelope) = plt.subplots(1, 2, sharey='row')  # Create a figure with two sub-plots side by side (sub_carrier, sub_envelope)

    ## Left hand plot
    sub_carrier.plot(time, carrier, color='black')
    # Set limits and aspect ratio
    sub_carrier.set_xlim(-1.5, 6.5)
    sub_carrier.set_aspect(1)
    # Titles and labels
    sub_carrier.set_title(r'Carrier Signal')
    sub_carrier.set_xlabel(r'Time (Seconds)')
    sub_carrier.set_ylabel(r'Intensity (Unitless)')

    ## Right hand plot
    sub_envelope.plot(time, envelope)
    sub_envelope.set_xlim(-1.5, 6.5)  # Set limits
    sub_envelope.set_aspect(1)  # Set aspect ratio
    # Titles and labels
    sub_envelope.set_title(r'|Envelope Signal|')
    sub_envelope.set_xlabel(r'Time (Seconds)')

    plt.tight_layout(pad=1.2)  # Place everything slightly closer together

    # plt.savefig('../Graphs/gaussian-decon.pgf', bbox_inches='tight', pad_inches=0)  # Save graph in this location

    ## FULL PULSED GAUSSIAN PLOT
    plt.figure()

    plt.plot(time, gaussian_pulse, color='green')
    plt.plot(time, envelope, color='orange', linestyle='dotted')
    plt.plot(time, -envelope, color='orange', linestyle='dotted')

    # Plot two lines to point out the envelope signal, these will be curved lines
    plt.annotate(r'Envelope signal, A(r,t)', xy=(1.35, 0.5), xycoords='data', xytext=(-0.5, 0.75), textcoords='data', ha="center", va="center",
                 arrowprops=dict(arrowstyle="->", color='black', connectionstyle="arc3,rad=0.15"))
    plt.annotate(r'Envelope signal, A(r,t)', xy=(1.35, -0.5), xycoords='data', xytext=(-0.5, 0.75), textcoords='data', ha="center", va="center",
                 arrowprops=dict(arrowstyle="->", color='black', connectionstyle="arc3,rad=0.5"))

    # Plot one line to point out actually propagating signal
    plt.annotate(r'Signal Propagating in Fibre, E(r,t)', xy=(2.2, -1), xycoords='data', xytext=(0.25, -1.75), textcoords='data', ha="center", va="center",
                 arrowprops=dict(arrowstyle="->", color='black'))

    # Set limits
    plt.xlim(-1.5, 6.5)

    # Aspect Ratio
    axes = plt.gca()
    axes.set_aspect(0.65)

    # Titles and labels
    plt.title(r'Pulsed Gaussian Signal')
    plt.xlabel(r'Time (Seconds)')
    plt.ylabel(r'Intensity (Unitless)')

    plt.tight_layout(pad=1.2)  # Place everything slightly closer together

    # plt.savefig('../Graphs/gaussian-demo.pgf', bbox_inches='tight', pad_inches=0)  # Save graph in this location

    ## PULSED GAUSSIAN SPECTRUM
    plt.figure()

    # Plot the normalised angular frequency (omega/omega_0), we will do the calculation for the x-axis inside the argument of the plot
    plt.plot(ang_freq_axis, np.abs(gaussian_pulse_spectrum), color='green')

    # Titles and labels
    plt.title(r'Pulsed Gaussian Spectrum')
    plt.ylabel(r'$\tilde{E}(\omega)$')
    plt.xlabel(r'Angular Frequency, $\omega$')

    # Set limits
    plt.xlim(-30, 30)

    # Aspect Ratio
    axes = plt.gca()
    axes.set_aspect(16)

    plt.tight_layout(pad=1.2)  # Place everything slightly closer together

    # plt.savefig('../Graphs/gaussian-demo-spectrum.pgf', bbox_inches='tight', pad_inches=0)  # Save graph in this location

    ## PULSED GAUSSIAN NORMALISED SPECTRUM
    plt.figure()

    # Plot the normalised angular frequency (omega/omega_0), we will do the calculation for the x-axis inside the argument of the plot
    plt.plot(ang_freq_axis / ang_freq, np.abs(gaussian_pulse_spectrum), color='green')

    # Titles and labels
    plt.title(r'Pulsed Gaussian Normalised Spectrum')
    plt.ylabel(r'$\tilde{E}(\omega)$')
    plt.xlabel(r'Normalised Angular Frequency, $\frac{\omega}{\omega_0}$')

    # Plot Delta Omega on graph
    plt.annotate(r'\textbf{${\Delta\omega}\ll{\omega_0}$}', xy=(-1.07, 0.367), xycoords='data', xytext=(-0.75, 0.367), textcoords='data', ha="center", va="center",
                 arrowprops=dict(arrowstyle="->", color='black'))
    plt.annotate(r'', xy=(-0.93, 0.367), xycoords='data', xytext=(-1, 0.367), textcoords='data', ha="center", va="center",
                 arrowprops=dict(arrowstyle="->", color='black'))

    # Set limits
    plt.xlim(-1.3, 1.3)

    # Aspect Ratio
    # axes = plt.gca()
    # axes.set_aspect(0.85)

    plt.tight_layout(pad=1.2)  # Place everything slightly closer together

    # plt.savefig('../Graphs/gaussian-demo-normalised-spectrum.pgf', bbox_inches='tight', pad_inches=0)  # Save graph in this location



def discrete_vs_continous_2d():
    ########## PLOT DATA CREATION ##########
    no_of_samples = 5000  # This number of samples will give a nice curved gaussian shape in the FFT

    intensity = 2  # Signal intensity
    init_waist = 1  # Used to define the Gaussian signal, this is the half-waist size at 1/e intensity

    signal_duration = 16

    time_samp_int = np.zeros((no_of_samples,), dtype=int)

    time = np.linspace(-signal_duration / 2, signal_duration / 2,
                       no_of_samples)  # Create an array of time values with them being 50000 equally-spaced values, the large array size makes for a better frequency domain tranform

    envelope = intensity * np.exp(-(((time - 3.25) ** 2) / (2 * (
            init_waist ** 2))))  # The Gaussian (envelope) signal itself, this is shifted to the right by the -3.25 term as a real signal might be (it also makes the graph pretty)

    ########## GRAPHING ##########
    ## TIME VECTOR DEMONSTRATION
    plt.figure()

    plt.plot(time, time_samp_int, color='orange', linestyle='dotted')

    # Set limits
    plt.xlim(-0.5, 6.5)

    # Aspect Ratio
    # axes = plt.gca()
    # axes.set_aspect(0.65)

    # Titles and labels
    plt.title(r'Sampling Vector')
    plt.xlabel(r'Time (Seconds)')

    plt.tight_layout(pad=1.2)  # Place everything slightly closer together

    # plt.savefig('../Graphs/envelope-discrete-vs-continuous.pgf', bbox_inches='tight', pad_inches=0)  # Save graph in this location

    ## ENVELOPE CONTINUOUS VS DISCRETE
    plt.figure()

    plt.plot(time, envelope, color='green')
    plt.plot(time, -envelope, color='green')
    plt.plot(time, envelope, color='orange', linestyle='dotted')
    plt.plot(time, -envelope, color='orange', linestyle='dotted')

    # Set limits
    plt.xlim(-0.5, 6.5)

    # Aspect Ratio
    axes = plt.gca()
    axes.set_aspect(0.65)

    # Titles and labels
    plt.title(r'Envelope of the Signal')
    plt.xlabel(r'Time (Seconds)')
    plt.ylabel(r'Intensity (Unitless)')

    plt.tight_layout(pad=1.2)  # Place everything slightly closer together

    # plt.savefig('../Graphs/envelope-discrete-vs-continuous.pgf', bbox_inches='tight', pad_inches=0)  # Save graph in this location




def discrete_vs_continous_3d():
    #####################################################################################################################
    ##### THE CODE TO IMPLEMENT THE SPLIT-STEP FOURIER METHOD #####

    ## The Fibre Material Properties
    beta_2 = 0.25  #
    gamma = 0  #

    ## Temporal Co-Ordinate Creation
    # The higher the number of samples for a given signal duration, the better the curved gaussian shape will be after the
    #   FFT:
    no_of_samples = 21
    signal_duration = 40

    # Create an array of time values with X number of equally-spaced values, a larger array size for a given signal
    #   duration makes for a better frequency domain transform
    time = np.linspace(-signal_duration / 2, signal_duration / 2, no_of_samples)

    # This is the period between each discrete sample of the signal in the time domain
    sample_period = signal_duration / (no_of_samples - 1)

    ## Spatial (z-axis) Co-Ordinate Creation
    signal_propagation_distance = 10

    # Create an array of z-axis values with X number of equally-spaced values, a larger array size for a given signal
    #   duration makes for a better frequency domain transform
    z_axis = np.linspace(-signal_propagation_distance / 2, signal_propagation_distance / 2, no_of_samples)

    # This is the distance between each discrete sample of the signal on the z-axis
    z_step_distance = signal_propagation_distance / (no_of_samples - 1)

    ## Frequency Domain Co-Ordinate Creation
    # Angular frequency
    # This is complicated but the angular frequency interval (distance between each measured discrete frequency of
    #   the FFT) is given by this calculation
    ang_freq_interval = (2.0 * np.pi) / signal_duration
    # From the above we can then create a frequency axis with defined frequency points based on the number of
    #   samples, this will be explained in the notes
    ang_freq_axis = np.arange(-no_of_samples / 2, no_of_samples / 2) * ang_freq_interval

    ## Input Envelope Definition
    intensity = 1  # Signal intensity
    init_waist = 1  # Used to define the Gaussian signal, this is the half-waist size at 1/e intensity

    # The Gaussian (envelope) signal
    envelope = intensity * np.exp(-((time ** 2) / (2 * (init_waist ** 2))))
    envelope = np.array(envelope)

    # Define null value vector for sampling mesh demonstration
    time_samp_int = np.zeros((no_of_samples,), dtype=int)


    ## Input Spectrum Calculations
    # First the raw intensity FFT transform, that will be arranged in order from -ve frequency to +ve frequency

    # The FFT algorithm in python is missing a component (sample_period/sqrt(2*pi)) see:
    #   https://cvarin.github.io/CSci-Survival-Guide/fft.html
    envelope_raw_spectrum = fftshift(fft(envelope)) * sample_period / np.sqrt(2 * np.pi)
    envelope_raw_spectrum = np.array(envelope_raw_spectrum)

    # Then one calculates the absolute intensity (|E(omega)|^2) for the input spectrum
    envelope_abs_spectrum = abs(envelope_raw_spectrum) ** 2

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
        current_step_time = ifft(ifftshift(current_step_spectrum)) * ang_freq_interval / np.sqrt(2.0 * np.pi) * no_of_samples

        for j in range(no_of_samples):
            # Non-Linear Step
            current_step_time[j] = current_step_time[j] * np.exp(complex(0, (non_linear_term *
                                                                             (abs(current_step_time[j]) ** 2))))

        pulse_output_time = np.row_stack((pulse_output_time, current_step_time))
        pulse_output_spectrum = np.row_stack((pulse_output_spectrum, current_step_spectrum))

    #####################################################################################################################
    ########## GRAPHING ##########

    ## SAMPLING MESH DEMONSTRATION

    T, Z = np.meshgrid(time, z_axis)

    # Create figure and axis objects
    fig = plt.figure()
    ax = plt.axes(projection='3d')

    for i in range(no_of_samples):
        for j in range(no_of_samples):
            ax.scatter(T[j, :], Z[i, :], time_samp_int)

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
    ax.set_zlabel(r'Sampling Points')

    # Aspect ratio settings
    ax.auto_scale_xyz([-22, 22], [-5, 5], [0, 0.1])
    ax.set_box_aspect((1.5, 2.0, 0.6))

    plt.tight_layout(pad=0.8)  # Place everything slightly closer together

    # Save graph in this location
    # plt.savefig('../Graphs/split-step-intro-3d-time.pgf', bbox_inches = 'tight', pad_inches = 0)

    #
    #
    # ## 3D FREQUENCY DOMAIN GRAPH
    # # Create a mesh matrix to plot signal values against
    # S, Zs = np.meshgrid(ang_freq_axis, z_axis)
    #
    # # Create figure and axis objects
    # fig = plt.figure()
    # ax = plt.axes(projection='3d')
    #
    # # Plot the signal as it propagates
    # ax.plot_surface(S, Zs, abs(pulse_output_spectrum)**2, rstride=1, cstride=1, cmap='viridis', edgecolor='none')
    #
    # # Colouring the background of the graph
    # fig.patch.set_facecolor('white')
    # ax.set_facecolor('white')
    # ax.w_xaxis.set_pane_color((0.95, 0.95, 0.95, 0.95))
    # ax.w_yaxis.set_pane_color((0.95, 0.95, 0.95, 0.95))
    # ax.w_zaxis.set_pane_color((0.95, 0.95, 0.95, 0.95))
    #
    # # Title and axis labels
    # ax.set_title(r'3D Plot of Signal Propagation in Freq. and Space')
    # ax.set_xlabel(r'Angular Frequency, $\omega$')
    # ax.set_ylabel(r'Propagation distance, z')
    # ax.set_zlabel(r'$\tilde{F}(\omega)$')
    #
    # # Aspect ratio settings
    # ax.auto_scale_xyz([-22, 22], [-5, 5], [0, 1.1])
    # ax.set_box_aspect((1.5, 2.0, 0.6))
    #
    # plt.tight_layout(pad=0.8) # Place everything slightly closer together
    #
    # # Save graph in this location
    # # plt.savefig('../Graphs/split-step-intro-3d-freq.pgf', bbox_inches = 'tight', pad_inches = 0)







def discrete_vs_continous():
    discrete_vs_continous_2d()
    discrete_vs_continous_3d()





def main():
    # gaussian_demo()
    discrete_vs_continous()
    plt.show()


if __name__ == "__main__":
    main()