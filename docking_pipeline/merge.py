import os
from tqdm import tqdm
from schrodinger.structure import StructureReader, StructureWriter

# Define the parent directory and the output file path
parent_directory = '/media/E/integrin_docking/no_midas'
output_file = '/media/E/integrin_docking/no_midas/merged.maegz'
batch_size = 10  # Adjust based on your system's memory capacity

def get_maegz_files(parent_directory):
    """Retrieve all .maegz files from subdirectories named split_X."""
    maegz_files = []
    for subdir, _, files in os.walk(parent_directory):
        # Check if the directory name matches the pattern split_X
        if os.path.basename(subdir).startswith('split_'):
            for file_name in files:
                if file_name.endswith('.maegz'):
                    maegz_files.append(os.path.join(subdir, file_name))
    return maegz_files

def process_batch(maegz_files, start_index, end_index, writer, first_file):
    """Process a batch of .maegz files."""
    for i in range(start_index, end_index):
        file_path = maegz_files[i]
        print(f"Processing file: {file_path}")

        with StructureReader(file_path) as reader:
            if first_file:
                # Read the first .maegz file to get the protein
                first_file_protein = next(reader)
                writer.append(first_file_protein)
                first_file = False
            else:
                # Skip the first structure (assumed to be protein) in subsequent files
                next(reader)  # Skip the first entry
                for struct in reader:
                    writer.append(struct)

def concatenate_maegz_files(parent_directory, output_file):
    # Retrieve all .maegz files
    maegz_files = get_maegz_files(parent_directory)
    
    # Open the output file for writing
    with StructureWriter(output_file) as writer:
        first_file = True
        num_files = len(maegz_files)
        num_batches = (num_files + batch_size - 1) // batch_size

        # Initialize progress bar
        with tqdm(total=num_files, desc="Processing .maegz files") as pbar:
            for batch in range(num_batches):
                start_index = batch * batch_size
                end_index = min(start_index + batch_size, num_files)
                
                # Process the current batch
                process_batch(maegz_files, start_index, end_index, writer, first_file)
                
                # Update the progress bar
                pbar.update(end_index - start_index)

    print(f"All .maegz files have been merged into {output_file}")

# Run the concatenation function
concatenate_maegz_files(parent_directory, output_file)

