#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Description:
    Cortico-Basal Ganglia Network Model implemented in PyNN using the NEURON simulator.
    This version of the model runs to a steady state and implements either DBS
    ampltiude or frequency modulation controllers, where the beta ARV from the STN LFP
    is calculated at each controller call and used to update the amplitude/frequency of
    the DBS waveform that is applied to the network.

    Full documentation of the model and controllers used is given in:
    https://www.frontiersin.org/articles/10.3389/fnins.2020.00166/

Original author: John Fleming, john.fleming@ucdconnect.ie
"""

import os
from pathlib import Path

# Change working directory so that imports works (save old)
oldpwd = Path(os.getcwd()).resolve()
newpwd = Path(__file__).resolve().parent
os.chdir(newpwd)

# No GUI please
opts = os.environ.get("NEURON_MODULE_OPTIONS", "")
if "nogui" not in opts:
    os.environ["NEURON_MODULE_OPTIONS"] = opts + " -nogui"

from mpi4py import MPI
import neuron
from pyNN.neuron import setup, run_until, end, simulator, reset 
from pyNN.parameters import Sequence
from Controllers import (
    ZeroController,
    StandardPIDController,
    IterativeFeedbackTuningPIController,
)

import neo.io
import quantities as pq
import numpy as np
import scipy.io as sio
import math
import argparse
from utils import make_beta_cheby1_filter, calculate_avg_beta_power, generate_dbs_pulses
from model import load_network, electrode_distance
from config import Config, get_controller_kwargs

# Import global variables for GPe DBS
import Global_Variables as GV

h = neuron.h
comm = MPI.COMM_WORLD


if __name__ == "__main__":

    os.chdir(oldpwd)
    parser = argparse.ArgumentParser(description="CBG Model")
    parser.add_argument("config_file", nargs="?", help="yaml configuration file")
    parser.add_argument(
        "-o", "--output-dir", default="RESULTS", help="output directory name"
    )
    args, unknown = parser.parse_known_args()

    #config_file = Path("conf_PTS_30s.yml").resolve()
    output_dir = Path(args.output_dir).resolve()
    c = Config("conf_baseline_30s.yml")
    os.chdir(newpwd)

    simulation_output_dir = output_dir

    simulation_runtime = c.RunTime
    controller_type = c.Controller
    rng_seed = c.RandomSeed
    timestep = c.TimeStep
    steady_state_duration = c.SteadyStateDuration
    save_stn_voltage = c.save_stn_voltage
    beta_burst_modulation_scale = c.beta_burst_modulation_scale
    burst_modulation_offset = c.modulation_offset

    sim_total_time = (
        steady_state_duration + simulation_runtime + timestep
    )  # Total simulation time
    rec_sampling_interval = 0.5  # Signals are sampled every 0.5 ms

    # Setup simulation
    rank = setup(timestep=timestep, rngseed=rng_seed)

    if rank == 0:
            print("\n------ Configuration ------")
            print(c, "\n")

    if rank == 0:
            print(f"Output directory: {simulation_output_dir}")
            simulation_output_dir.mkdir(parents=True, exist_ok=True)

    # Make beta band filter centred on 25Hz (cutoff frequencies are 21-29 Hz)
    # for biomarker estimation
    fs = 1000.0 / rec_sampling_interval
    beta_b, beta_a = make_beta_cheby1_filter(fs=fs, n=4, rp=0.5, low=21, high=29)

    
    def run_PTS_simulation(amplitude, phase):

        print("Simulation for stimulation amplitude {:.1f} mA and phase {:0.1f} deg.".format(amplitude, phase*360/12))
        print("rank: ", rank)

        # Use CVode to calculate i_membrane_ for fast LFP calculation
        cvode = h.CVode()
        cvode.active(0)

        # Second spatial derivative (the segment current) for the collateral
        cvode.use_fast_imem(1)

        # Set initial values for cell membrane voltages
        v_init = -68

        if rank == 0:
            print("Loading network...")
        (
            Pop_size,
            striatal_spike_times,
            Cortical_Pop,
            Interneuron_Pop,
            STN_Pop,
            GPe_Pop,
            GPi_Pop,
            Striatal_Pop,
            Thalamic_Pop,
            prj_CorticalAxon_Interneuron,
            prj_Interneuron_CorticalSoma,
            prj_CorticalSTN,
            prj_STNGPe,
            prj_GPeGPe,
            prj_GPeSTN,
            prj_StriatalGPe,
            prj_STNGPi,
            prj_GPeGPi,
            prj_GPiThalamic,
            prj_ThalamicCortical,
            prj_CorticalThalamic,
            GPe_stimulation_order,
        ) = load_network(
            steady_state_duration,
            sim_total_time,
            simulation_runtime,
            v_init,
            rng_seed,
            beta_burst_modulation_scale,
            burst_modulation_offset,
        )
        if rank == 0:
            print("Network loaded.")

        # Define state variables to record from each population
        Cortical_Pop.record("soma(0.5).v", sampling_interval=rec_sampling_interval)
        Cortical_Pop.record("collateral(0.5).v", sampling_interval=rec_sampling_interval)
        Interneuron_Pop.record("soma(0.5).v", sampling_interval=rec_sampling_interval)
        STN_Pop.record("soma(0.5).v", sampling_interval=rec_sampling_interval)
        STN_Pop.record("AMPA.i", sampling_interval=rec_sampling_interval)
        STN_Pop.record("GABAa.i", sampling_interval=rec_sampling_interval)
        Striatal_Pop.record("spikes")
        GPe_Pop.record("soma(0.5).v", sampling_interval=rec_sampling_interval)
        GPi_Pop.record("soma(0.5).v", sampling_interval=rec_sampling_interval)
        Thalamic_Pop.record("soma(0.5).v", sampling_interval=rec_sampling_interval)

        # Assign Positions for recording and stimulating electrode point sources
        recording_electrode_1_position = np.array([0, -1500, 250])
        recording_electrode_2_position = np.array([0, 1500, 250])
        stimulating_electrode_position = np.array([0, 0, 250])

        (
            STN_recording_electrode_1_distances,
            STN_recording_electrode_2_distances,
            Cortical_Collateral_stimulating_electrode_distances,
        ) = electrode_distance(
            recording_electrode_1_position,
            recording_electrode_2_position,
            STN_Pop,
            stimulating_electrode_position,
            Cortical_Pop,
        )

        # Conductivity and resistivity values for homogenous, isotropic medium
        sigma = 0.27  # Latikka et al. 2001 - Conductivity of Brain tissue S/m
        # rho needs units of ohm cm for xtra mechanism (S/m -> S/cm)
        rho = 1 / (sigma * 1e-2)

        # Calculate transfer resistances for each collateral segment for xtra
        # units are Mohms
        collateral_rx = (
            0.01
            * (rho / (4 * math.pi))
            * (1 / Cortical_Collateral_stimulating_electrode_distances)
        )

        # Convert ndarray to array of Sequence objects - needed to set cortical
        # collateral transfer resistances
        collateral_rx_seq = np.ndarray(
            shape=(1, Cortical_Pop.local_size), dtype=Sequence
        ).flatten()
        for ii in range(0, Cortical_Pop.local_size):
            collateral_rx_seq[ii] = Sequence(collateral_rx[ii, :].flatten())

        # Assign transfer resistances values to collaterals
        for ii, cell in enumerate(Cortical_Pop):
            cell.collateral_rx = collateral_rx_seq[ii]

        """
        # Create times for when the DBS controller will be called
        # Window length for filtering biomarker
        controller_window_length = 2000.0  # ms
        controller_window_length_no_samples = int(
            controller_window_length / rec_sampling_interval
        )
        # Window Tail length - removed post filtering, prior to
        # biomarker calculation
        controller_window_tail_length = 100.0  # ms
        controller_window_tail_length_no_samples = int(
            controller_window_tail_length / rec_sampling_interval
        )

        # Assuming that each time a stimulation is applied the controller must be called 
        controller_sampling_time = 20.0  # ms
        controller_start = (
            steady_state_duration + controller_window_length + controller_sampling_time
        )
        controller_call_times = np.arange(
            controller_start, sim_total_time, controller_sampling_time
        )

        if len(controller_call_times) == 0:
            controller_call_times = np.array([controller_start])
        """

        mat_dic = sio.loadmat('phase_t.mat') 
        stim_points = mat_dic['phase_t'][phase][0]
        stim_time_points = stim_points*rec_sampling_interval # convert the points into time points in ms 
        after_steady_state = np.array([stim_time_point > steady_state_duration for stim_time_point in stim_time_points])
        call_times = stim_time_points[after_steady_state]
        
        call_times = np.round(call_times, 1) 
        
        if controller_type == "ZERO":
            Controller = ZeroController
        elif controller_type == "PID":
            Controller = StandardPIDController
        elif controller_type == "IFT":
            Controller = IterativeFeedbackTuningPIController
        else:
            raise RuntimeError("Bad choice of Controller")

        controller_kwargs = get_controller_kwargs(c)
        controller = Controller(**controller_kwargs)

        # Generate a square wave which represents the DBS signal
        # Needs to be initialized to zero when unused to prevent
        # open-circuit of cortical collateral extracellular mechanism
        (
            DBS_Signal, 
            DBS_times
        ) = generate_dbs_pulses(steady_state_duration,
            stop_time=sim_total_time,
            stim_time_points=call_times,
            dt=simulator.state.dt,
            amplitude=amplitude, 
            pulse_width=0.06,        
        )
        
        # not sure of the size of the DBS_Signal so it can be concatenated in this way
        DBS_Signal = np.hstack((np.array([0, 0]), DBS_Signal))
        DBS_times = np.hstack((np.array([0, steady_state_duration + 10]), DBS_times))

        DBS_Signal_neuron = h.Vector(DBS_Signal)
        DBS_times_neuron = h.Vector(DBS_times)

        # Play DBS signal to global variable is_xtra
        DBS_Signal_neuron.play(h._ref_is_xtra, DBS_times_neuron, 1)

        # GPe DBS current stimulations - precalculated for % of collaterals
        # entrained for varying DBS amplitude
        interp_DBS_amplitudes = np.array(
            [0.0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2, 2.25, 2.50, 3, 4, 5]
        )
        interp_collaterals_entrained = np.array(
            [0, 0, 0, 1, 4, 8, 19, 30, 43, 59, 82, 100, 100, 100]
        )

        # Make new GPe DBS vector for each GPe neuron - each GPe neuron needs a
        # pointer to its own DBS signal
        GPe_DBS_Signal_neuron = []
        GPe_DBS_times_neuron = []
        updated_GPe_DBS_signal = []
        for i in range(0, Cortical_Pop.local_size):
            
            print("GPe neuron", i)

            (
                GPe_DBS_Signal, 
                GPe_DBS_times
            ) = generate_dbs_pulses(steady_state_duration,
                stop_time=sim_total_time,
                stim_time_points=call_times,
                dt=simulator.state.dt,
                amplitude=100, 
                pulse_width=0.06,        
            )

            GPe_DBS_Signal = np.hstack((np.array([0, 0]), GPe_DBS_Signal))
            GPe_DBS_times = np.hstack(
                (np.array([0, steady_state_duration + 10]), GPe_DBS_times)
            )

            # Set the GPe DBS signals to zero amplitude
            GPe_DBS_Signal[0:] = 0

            # Neuron vector of GPe DBS signals
            GPe_DBS_Signal_neuron.append(h.Vector(GPe_DBS_Signal))
            GPe_DBS_times_neuron.append(h.Vector(GPe_DBS_times))

            # Play the stimulation into each GPe neuron
            GPe_DBS_Signal_neuron[i].play(
                GV.GPe_stimulation_iclamps[i]._ref_amp, GPe_DBS_times_neuron[i], 1
            )

            # Hold a reference to the signal as a numpy array, and append to list
            # of GPe stimulation signals
            updated_GPe_DBS_signal.append(GPe_DBS_Signal_neuron[i].as_numpy())

        # DBS GPe neuron stimulation
        num_GPe_Neurons_entrained = int(
            np.interp(
                amplitude, interp_DBS_amplitudes, interp_collaterals_entrained
            )
        )
        # Stimulate the entrained GPe neurons
        for j in np.arange(0, num_GPe_Neurons_entrained):
            cellid = Cortical_Pop[GPe_stimulation_order[j]]

            if Cortical_Pop.is_local(cellid):
                index = Cortical_Pop.id_to_local_index(cellid)

                print("stimulated GPe neuron", j,'', index)

                (
                    GPe_DBS_Signal, 
                    GPe_DBS_times
                ) = generate_dbs_pulses(steady_state_duration,
                    stop_time=sim_total_time,
                    stim_time_points=call_times,
                    dt=simulator.state.dt,
                    amplitude=100, 
                    pulse_width=0.06,        
                )

                """
                GPe_DBS_Signal = np.hstack((np.array([0, 0]), GPe_DBS_Signal))
                GPe_DBS_times = np.hstack(
                    (np.array([0, steady_state_duration + 10]), GPe_DBS_times)
                )

                # Neuron vector of GPe DBS signals
                GPe_DBS_Signal_neuron.append(h.Vector(GPe_DBS_Signal))
                GPe_DBS_times_neuron.append(h.Vector(GPe_DBS_times))
                """

                # Play the stimulation into each GPe neuron
                GPe_DBS_Signal_neuron[index].play(
                    GV.GPe_stimulation_iclamps[index]._ref_amp, GPe_DBS_times_neuron[index], 1
                )

                # Hold a reference to the signal as a numpy array, and append to list
                # of GPe stimulation signals
                updated_GPe_DBS_signal.append(GPe_DBS_Signal_neuron[index].as_numpy())

        # Initialise STN LFP list
        STN_LFP = []
        STN_LFP_AMPA = []
        STN_LFP_GABAa = []

        # Variables for writing simulation data
        last_write_time = steady_state_duration

        if rank == 0:
            print(
                f"\n---> Running simulation to steady state ({steady_state_duration} ms) ..."
            )
        # Load the steady state
        run_until(steady_state_duration + simulator.state.dt, run_from_steady_state=False)
        if rank == 0:
            print("Steady state finished.")
            print(
                "\n---> Running simulation for %.0f ms after steady state (%.0f ms) with %s control"
                % (simulation_runtime, steady_state_duration, controller_type)
            )

        # Reload striatal spike times after loading the steady state
        Striatal_Pop.set(spike_times=striatal_spike_times[:, 0])
        
        run_until(sim_total_time - simulator.state.dt)

        if rank == 0:
            print("Controller Called at t: %.2f" % simulator.state.t)

        # Calculate the LFP and biomarkers, etc.
        STN_AMPA_i = np.array(
            STN_Pop.get_data("AMPA.i", gather=False).segments[0].analogsignals[0]
        )
        STN_GABAa_i = np.array(
            STN_Pop.get_data("GABAa.i", gather=False).segments[0].analogsignals[0]
        )
        STN_Syn_i = STN_AMPA_i + STN_GABAa_i


        """
        print("DEBUG STN_AMPA_i:", STN_AMPA_i.shape)
        print("DEBUG STN_AMPA_i[0]:", STN_AMPA_i.shape[0])
        #print("DEBUG STN_GABAa_i:", STN_GABAa_i.shape)
        print("DEBUG STN_AMPA_i[0][0]", STN_AMPA_i[0][0])
        print("DEBUG STN_AMPA_i[0]", STN_AMPA_i[:,0])
        #print("DEBUG STN_GABAa_i[0] at t: %.2f:" % simulator.state.t, STN_GABAa_i[0][0])
        print("DEBUG STN_AMPA_i[-1][0]", STN_AMPA_i[-1][0])
        #print("DEBUG STN_GABAa_i[-1] at t: %.2f :" % simulator.state.t, STN_GABAa_i[0][-1])
        print("DEBUG STN_Pop:", type(STN_Pop))
        """

        # Computing the STN_LFP recorded by the two point of the electrode and doing some sort of average perhpas
        # STN LFP Calculation - Syn_i is in units of nA -> LFP units are mV
        STN_LFP_1 = (
            (1 / (4 * math.pi * sigma))
            * np.sum(
                (1 / (STN_recording_electrode_1_distances * 1e-6))
                * STN_Syn_i.transpose(),
                axis=0,
            )
            * 1e-6
        )
        
        STN_LFP_2 = (
            (1 / (4 * math.pi * sigma))
            * np.sum(
                (1 / (STN_recording_electrode_2_distances * 1e-6))
                * STN_Syn_i.transpose(),
                axis=0,
            )
            * 1e-6
        )

        STN_LFP = np.hstack(
            (STN_LFP, comm.allreduce(STN_LFP_1 - STN_LFP_2, op=MPI.SUM))
        )

        # STN LFP AMPA and GABAa Contributions
        STN_LFP_AMPA_1 = (
            (1 / (4 * math.pi * sigma))
            * np.sum(
                (1 / (STN_recording_electrode_1_distances * 1e-6))
                * STN_AMPA_i.transpose(),
                axis=0,
            )
            * 1e-6
        )
        STN_LFP_AMPA_2 = (
            (1 / (4 * math.pi * sigma))
            * np.sum(
                (1 / (STN_recording_electrode_2_distances * 1e-6))
                * STN_AMPA_i.transpose(),
                axis=0,
            )
            * 1e-6
        )
        STN_LFP_AMPA = np.hstack(
            (STN_LFP_AMPA, comm.allreduce(STN_LFP_AMPA_1 - STN_LFP_AMPA_2, op=MPI.SUM))
        )

        STN_LFP_GABAa_1 = (
            (1 / (4 * math.pi * sigma))
            * np.sum(
                (1 / (STN_recording_electrode_1_distances * 1e-6))
                * STN_GABAa_i.transpose(),
                axis=0,
            )
            * 1e-6
        )
        STN_LFP_GABAa_2 = (
            (1 / (4 * math.pi * sigma))
            * np.sum(
                (1 / (STN_recording_electrode_2_distances * 1e-6))
                * STN_GABAa_i.transpose(),
                axis=0,
            )
            * 1e-6
        )
        STN_LFP_GABAa = np.hstack(
            (
                STN_LFP_GABAa,
                comm.allreduce(STN_LFP_GABAa_1 - STN_LFP_GABAa_2, op=MPI.SUM),
            )
        )

        # Min Jae isn't interested in this for the time being
        """
        # Biomarker Calculation:
        lfp_beta_average_value = calculate_avg_beta_power(
            lfp_signal=STN_LFP[-controller_window_length_no_samples:],
            tail_length=controller_window_tail_length_no_samples,
            beta_b=beta_b,
            beta_a=beta_a,
        )
        print("DEBUG len(STN_LFP):", len(STN_LFP))
        print("DEBUG simulator.state.t:", simulator.state.t)


        if rank == 0:
            print("Beta Average: %f" % lfp_beta_average_value)

        """ 

        # Write population data to file
        if save_stn_voltage:
            # comment 3
            #write_index = "{:.0f}_".format(call_index)
            write_index = "{:.0f}_".format(0)
            suffix = "_{:.0f}ms-{:.0f}ms".format(
                last_write_time, simulator.state.t)
            fname = write_index + "STN_Soma_v" + suffix + ".mat"
            dir_suffix = "_{:.1f}mA-{:.0f}deg".format(
                amplitude, phase*360/12)
            dir_name = "STN_POP" + dir_suffix 
            # testing without clear=True to check if we still have extra samples 
            STN_Pop.write_data(
                str(simulation_output_dir / dir_name / fname),
                "soma(0.5).v", clear=True
            )
        else:
            STN_Pop.get_data("soma(0.5).v", clear=True)

        last_write_time = simulator.state.t


        # Write population membrane voltage data to file
        if c.save_ctx_voltage:
            Cortical_Pop.write_data(str(simulation_output_dir / "Cortical_Pop/Cortical_Collateral_v.mat"), 'collateral(0.5).v', clear=False)
            Cortical_Pop.write_data(str(simulation_output_dir / "Cortical_Pop/Cortical_Soma_v.mat"), 'soma(0.5).v', clear=True)
        # Interneuron_Pop.write_data(str(simulation_output_dir / "Interneuron_Pop/Interneuron_Soma_v.mat"), 'soma(0.5).v', clear=True)
        # GPe_Pop.write_data(str(simulation_output_dir / "GPe_Pop/GPe_Soma_v.mat", 'soma(0.5).v'), clear=True)
        # GPi_Pop.write_data(str(simulation_output_dir / "GPi_Pop/GPi_Soma_v.mat", 'soma(0.5).v'), clear=True)
        # Thalamic_Pop.write_data(str(simulation_output_dir / "Thalamic_Pop/Thalamic_Soma_v.mat"), 'soma(0.5).v', clear=True)

        suffix = "_{:.1f}mA-{:.1f}deg".format(amplitude, phase*360/12)

        # Write the STN LFP to .mat file
        STN_LFP_Block = neo.Block(name="STN_LFP")
        STN_LFP_seg = neo.Segment(name="segment_0")
        STN_LFP_Block.segments.append(STN_LFP_seg)
        STN_LFP_signal = neo.AnalogSignal(
            STN_LFP,
            units="mV",
            t_start=0 * pq.ms,
            sampling_rate=pq.Quantity(1.0 / rec_sampling_interval, "1/ms"),
        )
        STN_LFP_seg.analogsignals.append(STN_LFP_signal)

        print("DEBUG final STN_LFP len:", len(STN_LFP))
        print("DEBUG simulator.state.t:", simulator.state.t)

        fname = "STN_LFP" + suffix +  ".mat"
        w = neo.io.NeoMatlabIO(filename=str(simulation_output_dir / fname))
        w.write_block(STN_LFP_Block)

        # Write the DBS Signal to .mat file
        # DBS Amplitude
        DBS_Block = neo.Block(name="DBS_Signal")
        DBS_Signal_seg = neo.Segment(name="segment_0")
        DBS_Block.segments.append(DBS_Signal_seg)
        DBS_signal = neo.AnalogSignal(
            DBS_Signal_neuron,
            units="mA",
            t_start=0 * pq.ms,
            sampling_rate=pq.Quantity(1.0 / simulator.state.dt, "1/ms"),
        )
        DBS_Signal_seg.analogsignals.append(DBS_signal)
        DBS_times = neo.AnalogSignal(
            DBS_times_neuron,
            units="ms",
            t_start=DBS_times_neuron * pq.ms,
            sampling_rate=pq.Quantity(1.0 / simulator.state.dt, "1/ms"),
        )
        DBS_Signal_seg.analogsignals.append(DBS_times)

        fname = "DBS_Signal" + suffix +  ".mat"
        w = neo.io.NeoMatlabIO(filename=str(simulation_output_dir / fname))
        w.write_block(DBS_Block)

        if rank == 0:
            print("Simulation Done!")
        
        reset()

    # Defining the list of amplitudes and phases we want to stimulate for
    amplitudes = np.arange(4, 5, 1)
    # 0 and 5 corresond respectively to stimulating at the peaks and at the trough
    phases = np.arange(0, 1, 1)

    # Iterate over different stimulation amplitudes
    for amplitude in amplitudes:  
        for phase in phases:
            run_PTS_simulation(amplitude, phase)

end()