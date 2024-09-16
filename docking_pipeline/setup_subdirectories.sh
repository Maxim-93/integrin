#!/bin/bash

# Set necessary directories
ligprep_dir="/sc/arion/work/epstem06/ligprep_ZINC/split_maegz_files"
output_dir="/sc/arion/projects/H_filizm02a/epstem06/glide-dock_SP/no_MIDAS/split_ZINC"
gridfile="/sc/arion/work/epstem06/glide-dock_SP/no_MIDAS/glide-grid_no_MIDAS.zip"
template_in="glide-dock_SP.in"

# Iterate through all split directories
for split_dir in "${ligprep_dir}/split_"*; do
    # Extract the split number (e.g., split_1, split_2)
    split_name=$(basename "${split_dir}")

    # Extract the number from the split directory name (e.g., 1, 2)
    split_number=$(echo "${split_name}" | grep -o '[0-9]\+')

    # Define paths for the ligand file, output subdirectory, and input file
    ligand_file="${split_dir}/${split_name}.maegz"
    subdir="${output_dir}/${split_name}"
    input_file="${subdir}/${template_in}"

    # Create the output subdirectory if it doesn't exist
    mkdir -p "${subdir}"

    # Copy the template .in file to the output subdirectory
    cp "${template_in}" "${input_file}"

    # Modify the input file to point to the correct ligand file
    sed -i "s|^LIGANDFILE.*|LIGANDFILE   ${ligand_file}|" "${input_file}"
    sed -i "s|^GRIDFILE.*|GRIDFILE   ${gridfile}|" "${input_file}"

    # Create a unique temporary directory for this job
    temp_dir="/sc/arion/projects/H_filizm02a/epstem06/tmp_schrodinger/${split_name}"
    mkdir -p "${temp_dir}"

    # Create the submission script
    submission_script="${subdir}/glide-dock_SP.sh"
    echo "#!/bin/bash" > "${submission_script}"
    echo "export SCHRODINGER_TEMP=${temp_dir}" >> "${submission_script}"
    echo "\"${SCHRODINGER}/glide\" ${input_file} -OVERWRITE -adjust -TMPDIR \$SCHRODINGER_TEMP -HOST localhost:20" >> "${submission_script}"

    # Make the submission script executable
    chmod +x "${submission_script}"

    echo "Setup completed for ${split_name}."
done
