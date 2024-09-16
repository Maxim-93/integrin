#!/bin/bash

# Set the output directory where the submission subdirectories are located
output_dir="/sc/arion/projects/H_filizm02a/epstem06/glide-dock_SP/no_MIDAS/split_ZINC"

# Iterate through all the split subdirectories in the output directory
for subdir in "${output_dir}/split_"*; do
    # Extract the split number from the subdirectory name
    split_number=$(basename "${subdir}" | grep -o '[0-9]\+')

    # Define the path to the submission script
    submission_script="${subdir}/glide-dock_SP.sh"

    # Change the working directory to the current subdirectory
    cd "${subdir}"

    # Submit the job using bsub with the specified options, including a job name with the split number
    bsub < "${submission_script}" -J "glideSP_${split_number}" -P acc_filizm02a -W 48:00 -n 20 -R "rusage[mem=1000]" -R "span[hosts=1]" -q premium

    echo "Job submitted for ${split_number}."

    # Change back to the original directory
    cd - > /dev/null
done
