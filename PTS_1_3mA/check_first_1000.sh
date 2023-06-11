#!/bin/bash

# loop from 1 to 600
for i in $(seq 1 600)
do
    # generate the folder name
    foldername="PTS-${i}"

    # check if the folder exists
    if [ ! -d "$foldername" ]; then
        # if the folder does not exist, print out a message
        echo "Folder $foldername does not exist"
    fi
done

