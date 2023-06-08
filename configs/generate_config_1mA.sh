#!/bin/bash

# Read the line from the file and convert it to an array
IFS=',' read -r -a seeds <<< "$(cat random_seeds_50.txt)"

# Loop over the array of seeds
for (( j=0; j<=11; j++ )); do
    for i in "${!seeds[@]}"; do
        # Adjust index to start from 51
        index=$((i + j*50 + 1))

        # Create a new config file with the current seed
        cat > "pts_test_${index}.yml" << EOF
---
# Simulation
Controller: Zero
TimeStep: 0.01
RandomSeed: ${seeds[$i]}
SteadyStateDuration: 6000
RunTime: 30000.0

# Controller
SetPoint: 0
Ts: 0

# PTS 
Amplitude: 1.0
Phase: ${j}
EOF
    done
done

