import numpy as np
import scipy.signal as signal
from pyNN.parameters import Sequence


def generate_poisson_spike_times(
    pop_size, start_time, duration, fr, timestep, random_seed
):
    """generate_population_spike_times generates (N = pop_size) Poisson
    distributed spiketrains with firing rate fr.

    Example inputs:
        pop_size = 10
        start_time = 0.0		# ms
        end_time = 6000.0		# ms
        timestep = 1  			# ms
        fr = 1					# Hz
    """

    # Convert to sec for calculating the spikes matrix
    dt = float(timestep) / 1000.0  # sec
    sim_time = float(((start_time + duration) - start_time) / 1000.0)  # sec
    n_bins = int(np.floor(sim_time / dt))

    spike_matrix = np.where(np.random.uniform(0, 1, (pop_size, n_bins)) < fr * dt)

    # Create time vector - ms
    t_vec = np.arange(start_time, start_time + duration, timestep)

    # Make array of spike times
    for neuron_index in np.arange(pop_size):
        neuron_spike_times = t_vec[
            spike_matrix[1][np.where(spike_matrix[0][:] == neuron_index)]
        ]
        if neuron_index == 0:
            spike_times = Sequence(neuron_spike_times)
        else:
            spike_times = np.vstack((spike_times, Sequence(neuron_spike_times)))

    return spike_times


def make_beta_cheby1_filter(fs, n, rp, low, high):
    """Calculate bandpass filter coefficients (1st Order Chebyshev Filter)"""
    nyq = 0.5 * fs
    lowcut = low / nyq
    highcut = high / nyq

    b, a = signal.cheby1(n, rp, [lowcut, highcut], "band")

    return b, a


def generate_dbs_pulses(start_time, stop_time, stim_time_points, dt, amplitude, pulse_width):

    # need same units as the output from generate_dbs_signal - probably output DBS pulse in nA and time in ms
    # dt = 0.01 ms
    times = np.round(np.arange(0, stop_time - start_time, dt), 2) # ms 
    DBS_signal = np.zeros(np.size(times))

    # shifting the times of of stim_time_points (to have them relative the begining of stimulation)
    stim_time_points = stim_time_points - start_time 
    # we remove the times that go over 30s in running time (problem due to the extra sample
    # when the controller is called -> for the first stimulation phase -> 17 last to discard)
    #print("DEBUG: overlap", list(set(stim_time_points.tolist()) - set(times.tolist())))
    #print("DEBUG: times[-50:]",times[-50:])
    # get the time index for which the phase is reached
    # do a shift of 30 µs in the past to know when to start stimulating
    #print("DEBUG: stim_time_point last 50", stim_time_points[-50:])
    # print("DEBUG: stim_time_point first 50 ", stim_time_points[:50])

    stim_time_idx = [np.where(times == stim_time_point - pulse_width/2.0)[0][0] for stim_time_point in stim_time_points]
    for idx in stim_time_idx:
       DBS_signal[idx: idx + np.intc(pulse_width/dt)] = 1
    
    # if cathodic stimulation amplitude < 0 (we are interest in anodic stimulations)
    DBS_signal *= amplitude

    return DBS_signal, times


def calculate_avg_beta_power(lfp_signal, tail_length, beta_b, beta_a):
    """Calculate the average power in the beta-band for the current LFP signal
    window, i.e. beta Average Rectified Value (ARV)

    Inputs:
        lfp_signal          - window of LFP signal (samples)

        tail_length         - tail length which will be discarded due to
                              filtering artifact (samples)

        beta_b, beta_a      - filter coefficients for filtering the beta-band
                              from the signal
    """

    lfp_beta_signal = signal.filtfilt(beta_b, beta_a, lfp_signal)
    lfp_beta_signal_rectified = np.absolute(lfp_beta_signal)
    avg_beta_power = np.mean(lfp_beta_signal_rectified[-2 * tail_length : -tail_length])

    return avg_beta_power