import os
from schrodinger import structure
from tqdm import tqdm

def split_maegz_file(input_file, max_ligands_per_file=10000, output_dir='split_maegz_files_10000'):
    # Create the output directory if it doesn't exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Initialize counters
    ligand_count = 0
    file_count = 1
     
    print("Open the input .maegz file")
    st_reader = structure.StructureReader(input_file)
    
    print("Estimate total number of structures for progress bar")
    total_structures = sum(1 for _ in structure.StructureReader(input_file))
    
    print("Initialize the first output file")
    output_subdir = os.path.join(output_dir, f'split_{file_count}')
    os.makedirs(output_subdir, exist_ok=True)
    output_file = os.path.join(output_subdir, f'split_{file_count}.maegz')
    st_writer = structure.StructureWriter(output_file)
    
    # Initialize the tqdm progress bar
    with tqdm(total=total_structures, desc="Splitting .maegz file", unit=" ligand") as pbar:
        print("Iterate over the structures in the input file")
        for st in st_reader:
            st_writer.append(st)
            ligand_count += 1
            pbar.update(1)
            
            # If the max ligands per file is reached, close the current file and start a new one
            if ligand_count >= max_ligands_per_file:
                st_writer.close()
                file_count += 1
                ligand_count = 0
                
                output_subdir = os.path.join(output_dir, f'split_{file_count}')
                os.makedirs(output_subdir, exist_ok=True)
                output_file = os.path.join(output_subdir, f'split_{file_count}.maegz')
                st_writer = structure.StructureWriter(output_file)
    
    # Close the last writer
    st_writer.close()
    st_reader.close()

    print(f"Splitting complete. Files saved in {output_dir}")

if __name__ == "__main__":
    input_file = 'ligprep_all_ZINC-out.maegz'
    split_maegz_file(input_file)
