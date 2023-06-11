#!/bin/bash

# loop from 1 to 600
for i in $(seq 1 600)
do
    # generate the filename
    filename="pts_test_${i}.yml"

    # check if file exists
    if [ -f "$filename" ]; then
        # if file exists, delete it
        rm "$filename"
        echo "Deleted $filename"
    else
        echo "File $filename does not exist"
    fi
done

