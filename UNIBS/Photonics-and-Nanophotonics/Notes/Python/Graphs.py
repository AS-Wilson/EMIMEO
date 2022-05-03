import numpy as np
import matplotlib.pyplot as plt
import os as sys
from scipy.fftpack import fft, fftshift, fftfreq

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
    no_of_samples = 5000000 # This number of samples will give a nice curved gaussian shape in the FFT
    
    intensity = 2 # Signal intensity
    init_waist = 1 # Used to define the Gaussian signal, this is the half-waist size at 1/e intensity
    freq = 4 # Carrier frequency
    ang_freq = 2 * np.pi * freq # Carrier angular frequency

    signal_duration = 50
    
    time = np.linspace(-signal_duration/2, signal_duration/2, no_of_samples) # Create an array of time values with them being 50000 equally-spaced values, the large array size makes for a better frequency domain tranform
    
    sample_period = signal_duration / no_of_samples # This is the period between each discrete sample of the signal
    
    carrier = intensity * np.sin(ang_freq * (time)) # The carrier (cosine) portion of the Pulsed Gaussian
    envelope = intensity * np.exp( -( ((time-3.25)**2) / (2*(init_waist**2)) )) # The Gaussian (envelope) signal itself, this is shifted to the right by the -3.25 term as a real signal might be (it also makes the graph pretty)

    gaussian_pulse = intensity * np.exp( -( ((time-3.25)**2) / (2*(init_waist**2)) )) * np.sin(ang_freq * (time))
    
    
    ## SPECTRUM CALCULATIONS
    # Angular frequency
    # This is complicated but the angular frequency interval (distance between each measured discrete frequency of the FFT) is given by this calculation
    ang_freq_interval = (2.0 * np.pi) / signal_duration 
    # From the above we can then create a frequency axis with discrete frequency points based on the number of samples, this will be explained in the notes
    ang_freq_axis = np.arange(-no_of_samples/2, no_of_samples/2) * ang_freq_interval 

    # Spectrum
    gaussian_pulse_spectrum = fftshift(fft(gaussian_pulse)) * sample_period/np.sqrt(2*np.pi) # The FFT algorithm in python is missing a component (sample_period/sqrt(2*pi)) see: https://cvarin.github.io/CSci-Survival-Guide/fft.html
    gaussian_pulse_spec = abs(gaussian_pulse_spectrum)**2 # This calculates the intensity (|E(omega)|^2) for the spectrum graph

    
    
    
    ########## GRAPHING ##########

    ## DECONSTRUCTED GAUSSIAN PLOT
    fig_gauss_decon, (sub_carrier, sub_envelope) = plt.subplots(1, 2, sharey='row') # Create a figure with two sub-plots side by side (sub_carrier, sub_envelope)

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
    sub_envelope.set_xlim(-1.5, 6.5) # Set limits
    sub_envelope.set_aspect(1) # Set aspect ratio
    # Titles and labels
    sub_envelope.set_title(r'|Envelope Signal|')
    sub_envelope.set_xlabel(r'Time (Seconds)')    
    
    plt.tight_layout(pad=1.2) # Place everything slightly closer together

    plt.savefig('../Graphs/gaussian-decon.pgf', bbox_inches = 'tight', pad_inches = 0) # Save graph in this location

        
    ## FULL PULSED GAUSSIAN PLOT
    plt.figure()
    
    plt.plot(time, gaussian_pulse, color='green')
    plt.plot(time, envelope, color='orange', linestyle='dotted')
    plt.plot(time, -envelope, color='orange', linestyle='dotted')

    # Plot two lines to point out the envelope signal, these will be curved lines
    plt.annotate(r'Envelope signal, A(r,t)', xy=(1.35, 0.5), xycoords='data', xytext=(-0.5, 0.75), textcoords='data', ha="center", va="center", arrowprops=dict(arrowstyle="->", color='black', connectionstyle="arc3,rad=0.15"))
    plt.annotate(r'Envelope signal, A(r,t)', xy=(1.35, -0.5), xycoords='data', xytext=(-0.5, 0.75), textcoords='data', ha="center", va="center", arrowprops=dict(arrowstyle="->", color='black', connectionstyle="arc3,rad=0.5"))

    # Plot one line to point out actually propagating signal
    plt.annotate(r'Signal Propagating in Fibre, E(r,t)', xy=(2.2, -1), xycoords='data', xytext=(0.25, -1.75), textcoords='data', ha="center", va="center", arrowprops=dict(arrowstyle="->", color='black'))
    
    # Set limits
    plt.xlim(-1.5, 6.5)
    
    # Aspect Ratio
    axes = plt.gca()
    axes.set_aspect(0.65)
    
    # Titles and labels
    plt.title(r'Pulsed Gaussian Signal')
    plt.xlabel(r'Time (Seconds)')
    plt.ylabel(r'Intensity (Unitless)')

    plt.tight_layout(pad=1.2) # Place everything slightly closer together

    plt.savefig('../Graphs/gaussian-demo.pgf', bbox_inches = 'tight', pad_inches = 0) # Save graph in this location
    
    
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

    plt.tight_layout(pad=1.2) # Place everything slightly closer together

    plt.savefig('../Graphs/gaussian-demo-spectrum.pgf', bbox_inches = 'tight', pad_inches = 0) # Save graph in this location
    
    
    ## PULSED GAUSSIAN NORMALISED SPECTRUM
    plt.figure()

    # Plot the normalised angular frequency (omega/omega_0), we will do the calculation for the x-axis inside the argument of the plot
    plt.plot(ang_freq_axis/ang_freq, np.abs(gaussian_pulse_spectrum), color='green')
    
    # Titles and labels
    plt.title(r'Pulsed Gaussian Normalised Spectrum')
    plt.ylabel(r'$\tilde{E}(\omega)$')
    plt.xlabel(r'Normalised Angular Frequency, $\frac{\omega}{\omega_0}$')

    # Plot Delta Omega on graph
    plt.annotate(r'\textbf{${\Delta\omega}\ll{\omega_0}$}', xy=(-1.07, 0.367), xycoords='data', xytext=(-0.75, 0.367), textcoords='data', ha="center", va="center", arrowprops=dict(arrowstyle="->", color='black'))
    plt.annotate(r'', xy=(-0.93, 0.367), xycoords='data', xytext=(-1, 0.367), textcoords='data', ha="center", va="center", arrowprops=dict(arrowstyle="->", color='black'))

    # Set limits
    plt.xlim(-1.3, 1.3)
    
    # Aspect Ratio
    #axes = plt.gca()
    #axes.set_aspect(0.85)

    plt.tight_layout(pad=1.2) # Place everything slightly closer together

    plt.savefig('../Graphs/gaussian-demo-normalised-spectrum.pgf', bbox_inches = 'tight', pad_inches = 0) # Save graph in this location
    
    


def discrete_vs_continous():
    
    ########## PLOT DATA CREATION ##########
    no_of_samples = 5000 # This number of samples will give a nice curved gaussian shape in the FFT
    
    intensity = 2 # Signal intensity
    init_waist = 1 # Used to define the Gaussian signal, this is the half-waist size at 1/e intensity
    freq = 4 # Carrier frequency
    ang_freq = 2 * np.pi * freq # Carrier angular frequency

    signal_duration = 16
    
    time = np.linspace(-signal_duration/2, signal_duration/2, no_of_samples) # Create an array of time values with them being 50000 equally-spaced values, the large array size makes for a better frequency domain tranform
    
    envelope = intensity * np.exp( -( ((time-3.25)**2) / (2*(init_waist**2)) )) # The Gaussian (envelope) signal itself, this is shifted to the right by the -3.25 term as a real signal might be (it also makes the graph pretty)
    
    
    ########## GRAPHING ##########
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

    plt.tight_layout(pad=1.2) # Place everything slightly closer together

    plt.savefig('../Graphs/envelope-discrete-vs-continuous.pgf', bbox_inches = 'tight', pad_inches = 0) # Save graph in this location
    
    
    
    
    
def main():
    # gaussian_demo()
    discrete_vs_continous()
    plt.show()


if __name__ == "__main__":
    main()
    
