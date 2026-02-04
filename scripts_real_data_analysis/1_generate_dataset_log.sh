#!/bin/bash

# Check if the directory is provided as a command line argument
if [[ -z $1 ]]; then
    echo "Please provide a directory as a command line argument."
    exit 1
fi

# Specify the directory and output file
directory=$1
base_annt=$2
obo_file=$3
output_file="1_datasets_log.ndjson"

if ! [ -e $output_file ]; then
    touch $output_file
fi


# Check if the first command line argument is 'regenerate'
if [[ $4 == "regenerate" ]]; then
    # If so, empty the output file
    > $output_file
    # And initialize the counter to 1
    counter=1
    echo "Regenerating $output_file."
else
    # Otherwise, read the last line of the output file to get the last counter
    last_counter=$(tail -n 1 "$output_file" | grep -oP '(?<="ID_dataset": )\d+')
    echo "Last counter $last_counter" 

    # If the output file is empty, initialize the counter to 1
    if [[ -z "$last_counter" ]]; then
        counter=1
    else
        # Otherwise, increment the last counter by 1
        counter=$((last_counter + 1))
    fi
fi

# Loop through each file in the directory
for filename in "$directory"/*; do
    # Check if the file name is already in the output file
    if ! grep -q -P "\"ID_dataset\": \d+, \"gene_list\": \"$filename\", " "$output_file"; then
        # If not, check if annt_file has been generated, if not then run the code
        # if [ ! -f "save_${counter}_0to100000_annotations.tsv" ]; then
            annt_file="save_${counter}_0to100000_annotations.tsv"
            # Record the directory and basename file name in the output file
            python3 scripts_real_data_analysis/1_preprocess.py -g $(realpath "$filename") -i $counter -a $base_annt -o $obo_file  && echo -e "{\"ID_dataset\": $counter, \"gene_list\": \"$filename\", \"annotations\": \"$annt_file\"}" >> "1_datasets_log.ndjson" && ((counter++))
            # append the counter and file name to the output file
            # echo -e "$counter\t$(basename "$filename")\t$annt_file" >> $output_file
            # Print these three fields in json format instead, with each line as a separate entry in one json object, spearated by commas, and the whole thing enclosed in square brackets
        else
            echo "Already existed"
        # fi
    fi
done

echo "New file names have been appended to $output_file."
