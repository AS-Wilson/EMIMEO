import numpy as np
import matplotlib.pyplot as plt
import os as sys
from scipy.fftpack import fft, fftshift, fftfreq


def gaussian_demo():
    ########## PLOT DATA CREATION ##########
    no_of_samples = 5000000
    
    intensity = 2 # Signal intensity
    init_waist = 1 # Used to define the Gaussian signal, this is the half-waist size at 1/e intensity
    freq = 2 # Carrier frequency
    ang_freq = 2 * np.pi * freq # Carrier angular frequency

    signal_duration = 50
    
    time = np.linspace(0, signal_duration, no_of_samples) # Create an array of time values between -5 and 5 seconds, with 50000 equally-spaced values
    
    sample_period = time[1] - time[0] # This is the period between each discrete sample of the signal
    
    gaussian_pulse = intensity * np.exp( -( ((time-25)**2) / (2*(init_waist**2)) )) * np.sin(ang_freq * (time))

    
    
    ## SPECTRUM CALCULATIONS
    # Angular frequency
    ang_freq_interval = (2.0 * np.pi) / signal_duration
    ang_freq_axis = np.arange(-no_of_samples/2, no_of_samples/2) * ang_freq_interval

    # Spectrum
    gaussian_pulse_spectrum = fftshift(fft(gaussian_pulse)) * sample_period/np.sqrt(2*np.pi) # The FFT algorithm is missing a component (sample_period/sqrt(2*pi))
    gaussian_pulse_spec = abs(gaussian_pulse_spectrum)**2 # This calculates the intensity (|E(omega)|^2)
    
    
    ## PULSED GAUSSIAN SPECTRUM
    plt.figure()
    
    plt.plot(ang_freq_axis/ang_freq, np.abs(gaussian_pulse_spectrum), color='green')
    
    # Titles and labels
    plt.title('Pulsed Gaussian Signal')
    plt.ylabel('E(omega)')
    plt.xlabel('Normalised Frequency')

    # Set limits
    plt.xlim(-1.3, 1.3)
    
    # Aspect Ratio
    axes = plt.gca()
    axes.set_aspect(0.85)

    plt.tight_layout(pad=1.2) # Place everything slightly closer together
    
    
    
    
def main():
    gaussian_demo()
    plt.show()


if __name__ == "__main__":
    main()
    
