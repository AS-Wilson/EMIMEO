# To use this module just type include at start of your file:
# from OCNToolbox import *

import numpy as np
import matplotlib.pyplot as plt
from scipy.fftpack import fft, fftshift, ifft, ifftshift
from mpl_toolkits.mplot3d import Axes3D


def get_spectrum(signal, sample_period):
    ## Input Spectrum Calculations
    # First the raw intensity FFT transform, that will be arranged in order from -ve
    # frequency to +ve frequency
    # The FFT algorithm in python is missing a component (sample_period/sqrt(2*pi)) see:
    #   https://cvarin.github.io/CSci-Survival-Guide/fft.html
    raw_spectrum = fftshift(fft(signal)) * sample_period

    # Then one calculates the absolute intensity (|E(omega)|^2) for the input spectrum
    abs_spectrum = abs(raw_spectrum)

    return raw_spectrum, abs_spectrum


def gaussian(time, t_0=None, FWHM=None, intensity=1, t_shift=0):
    """Public Method for producing a Gaussian Signal after providing various parameters."""

    try:
        t_0 = FWHM / (2 * np.sqrt(2 * np.log(2)))
        signal = intensity * np.exp(- ((time + t_shift) ** 2 / 2 * t_0 ** 2))
        return signal

    except TypeError:
        signal = intensity * np.exp(- ((time + t_shift) ** 2 / 2 * t_0 ** 2))
        return signal


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


class Split_Step_Fibre_Propagation:
    """The propagation of a given input signal down a fibre in a given amount of
        propagation distance and time.
        Based on the Split-Step Fourier Method

        Parameters
        ==========


        time : The matrix of time values used to calculate the input signal
            and using which each step of the split-step fourier method will
            be calculated.
        signal_distance : The distance of fibre along which to calculate the
            Split-Step Fourier method. This will be used to produce the array of
            z-axis (propagation axis) values and is required to caluclate the step
            size of the Split-Step Fourier method.
        input_signal : This is the matrix of intensity values for the input
            signal under study. Length of the matrix must match the number of samples.
        beta_2 : The dispersion term value of the fibre under study. This is an
            inherent material property of the fibre. If undefined it will be set to 0.
        gamma : The non-linear term of the fibre under study. This is an
            inherent material property of the fibre. If undefined it will be set to 0.
        no_of_z_samples : This is the number of discrete steps used to
            calculate for the space axis for the given input pulse. If this is undefined
            the time matrix length is used instead.


        Attributes
        ==========

        signal_duration : The duration (amount of time), in seconds, that the input
            signal was calulated through. This is required to calculate the sample
            period and angular frequency of the results from the class.
        sample_period : The duration between each discrete time sample, required for
            frequency calculations.
        z_axis : The propagation axis values
        z_step_size : This is the distance between each discrete sample of the signal on
            the z-axis.
        ang_freq_interval : The distance between each measured discrete frequency of the
            FFT, required to calculate each discrete frequency in the ang_freq_axis.
        ang_freq_axis : The frequency axis which is a set of discrete frequency points
            based on the number of samples in the time matrix.
        non_linear_term :
        dispersive_terms :
        raw_input_spec :
        abs_input_spec :
        sig_propagation_time :
        sig_propagation_freq :


        Raises
        ======

        TypeError : When instantiating with anything but a Point or sequence
        ValueError : when instantiating with a sequence with length < 2 or
            when trying to reduce dimensions if keyword `on_morph='error'` is
            set.

        Example
        ========

        >>> import numpy as np
        >>> import matplotlib.pyplot as plt
        >>> import OCNToolbox
        # Time axis between -10 and 10 seconds with 1001 samples
        >>> time = np.linspace(-10 / 2, 10 / 2, 1001)
        # Create a Gaussian which is time shifted by 5 seconds and has width of 1
        >>> input_signal = 1 * np.exp(- ((time - 5) ** 2 / 2 * 1 ** 2))
        # Calculate the propagation of light down a fibre with no properties
        >>> normal = Split_Step_Fibre_Propagation(time, 10, input_signal)

        """

    def __init__(self, time, signal_distance, input_signal, no_of_z_samples=None,
                 beta_2=None, gamma=None):

        self.signal_propagation_time = None
        self.signal_propagation_spec = None
        self.abs_signal_propagation_spec = None
        self.input_signal = input_signal

        ## Create attributes related to Time Co-Ordinates
        self.time = time
        self.signal_duration = time[-1] - time[0]
        self.sample_period = self.signal_duration / (len(time) - 1)

        ## Spatial (z-axis) Co-Ordinate Creation
        # If no_of_z_samples is undefined we shall just use the time matrix length

        if no_of_z_samples is None:
            self.z_axis = np.linspace(0, signal_distance, len(time))
            self.z_step_size = signal_distance / (len(time) - 1)
            self.no_of_z_samples = len(time)

        else:
            self.z_axis = np.linspace(0, signal_distance, no_of_z_samples)
            self.z_step_size = signal_distance / (no_of_z_samples - 1)
            self.no_of_z_samples = no_of_z_samples

        ## Angular Frequency Co-Ordinate Creation
        # Defining the angular frequency axis and associated values
        self.ang_freq_interval = (2.0 * np.pi) / self.signal_duration
        self.ang_freq_axis = np.arange(-len(time) / 2, len(time) / 2) * self.ang_freq_interval

        if beta_2 is None and gamma is None:
            self.dispersive_terms = 0
            self.non_linear_term = 0
            pass

        elif beta_2 is None and gamma is not None:
            self.beta_2 = 0
            self.gamma = gamma

            self.calculate_propagation_terms()

        elif beta_2 is not None and gamma is None:
            self.beta_2 = beta_2
            self.gamma = 0

            self.calculate_propagation_terms()

        else:
            self.beta_2 = beta_2
            self.gamma = gamma

            self.calculate_propagation_terms()

        ## Get spectrum of input signal for graphing
        self.raw_input_spec, self.abs_input_spec = \
            get_spectrum(self.input_signal, self.sample_period)

        self.split_step_fourier_method()

    def calculate_propagation_terms(self):
        """Something Here: TODO Later"""
        ## Factors used in the Split-Step Fourier Method
        self.non_linear_term = complex(0, (self.gamma * self.z_step_size))

        # Produce all the dispersive terms, they are frequency dependent
        for i in range(len(self.ang_freq_axis)):
            try:
                self.dispersive_terms = np.column_stack((self.dispersive_terms, np.array(
                    np.exp(complex(0, ((self.ang_freq_axis[i] ** 2 / 2) * self.beta_2 *
                                       self.z_step_size))))))

            except AttributeError:
                self.dispersive_terms = np.array(np.exp(
                    complex(0, ((self.ang_freq_axis[i] ** 2 / 2) * self.beta_2 *
                                self.z_step_size))))

    def split_step_fourier_method(self):
        """Something Here, TODO Later"""

        signal_propagation_spec, abs_signal_propagation_spec = \
            get_spectrum(self.input_signal, self.sample_period)

        # The actual Split-Step Fourier Method Calculations
        for i in range(0, self.no_of_z_samples - 1):

            try:
                signal_propagation_spec = np.vstack((signal_propagation_spec, ((fftshift(fft(signal_propagation_time[-1, :])) * self.sample_period) * self.dispersive_terms)))

            except UnboundLocalError:
                signal_propagation_spec = np.vstack((signal_propagation_spec, (signal_propagation_spec * self.dispersive_terms)))
                signal_propagation_time = self.input_signal.astype(complex)

            inter_prop_time = ifft(ifftshift(signal_propagation_spec[-1, :])) * (1/self.sample_period)

            signal_propagation_time = np.vstack((signal_propagation_time, (inter_prop_time * np.exp((self.non_linear_term * (abs(inter_prop_time) ** 2))))))

        self.signal_propagation_spec = signal_propagation_spec
        self.signal_propagation_time = signal_propagation_time
