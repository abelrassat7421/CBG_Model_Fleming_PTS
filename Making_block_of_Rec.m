% Set up parameters
amplitudes = 0:0.5:4;
phases = 0:30:330;
num_amplitudes = length(amplitudes);
num_phases = length(phases);
signal_size = 72001;

% Create empty block
STN_LFP_Recordings = zeros([num_phases, num_amplitudes, signal_size]);

% Loop over recordings for different amplitudes and phases, loading data into block
for jj = 1:num_amplitudes
    for kk = 1:num_phases
        directory_name = 'Full_simulations';
        folder_name = sprintf('Rec_%.1fmA_%ddeg', amplitudes(jj), phases(kk));
        file_name = sprintf('STN_LFP_%.1fmA-%.1fdeg.mat', amplitudes(jj), phases(kk));  
        full_path = fullfile(directory_name, folder_name, file_name);
        data = load(full_path);
        STN_LFP_Recordings(kk, jj, :) = data.block.segments{1, 1}.analogsignals{1, 1}.signal;  
    end
end