import os
import csv

# Parent directory containing the split_* directories
parent_dir = "/sc/arion/projects/H_filizm02a/epstem06/glide-dock_SP/with_MIDAS/split_ZINC"  # Change this to the actual path
output_file = "combined_glide_dock_SP.csv"

# Prepare to write the combined CSV file
with open(output_file, 'w') as outfile:  # Remove newline argument for Python 2 compatibility
    writer = None  # CSV writer
    header_written = False  # Track if header has been written

    # Loop through all directories matching "split_*"
    for subdir in sorted(os.listdir(parent_dir)):
        if subdir.startswith("split_") and os.path.isdir(os.path.join(parent_dir, subdir)):
            split_number = subdir.split('_')[1]  # Get the split number (e.g., '1' from 'split_1')
            csv_file = os.path.join(parent_dir, subdir, "glide-dock_SP.csv")

            if os.path.exists(csv_file):
                with open(csv_file) as infile:  # Open the file without 'newline' argument
                    reader = csv.reader(infile)
                    header = next(reader)  # Read the header

                    # If header not written to the output file, write it and add 'split file' column
                    if not header_written:
                        writer = csv.writer(outfile)
                        header.append('split file')
                        writer.writerow(header)
                        header_written = True

                    # Write the rest of the rows with the split number appended
                    for row in reader:
                        row.append(split_number)
                        writer.writerow(row)

print("Concatenation complete! Combined file saved as '{}'".format(output_file))
