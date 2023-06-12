% Set up parameters
amplitudes = 3;
phases = 0:30:330;
num_amplitudes = length(amplitudes);
num_phases = length(phases);
signal_size = 72001;

seeds = [1825,410,4507,4013,3658,2287,1680,8936,1425,9675,6913,521,489,1536,3583,3812,8280,9864,435,9196,3258,8929,6874,3612,7360,9655,4558,107,2616,6925,5575,4553,2548,3528,5515,1675,1520,6225,1585,5882,5636,9892,4334,712,7528,8786,2046,6202,1292,9045];
num_seeds = length(seeds);

% Create empty block
STN_LFP_Recordings_3mA = zeros([num_phases, num_seeds, signal_size]);
i = 600;

% Loop over recordings for different amplitudes and phases, loading data into block
for kk = 1:num_phases
    for jj = 1:num_seeds
        i = i+1;
        directory_name = 'PTS_1_3mA';
        folder_name = sprintf('PTS-%d', i);
        file_name = sprintf('STN_LFP_3.0mA-%.1fdeg-%dseed.mat', phases(kk), seeds(jj));  
        full_path = fullfile(directory_name, folder_name, file_name);
        data = load(full_path);
        STN_LFP_Recordings_3mA(kk, jj, :) = data.block.segments{1, 1}.analogsignals{1, 1}.signal;  
    end
end