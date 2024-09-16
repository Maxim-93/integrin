import os
import csv
from schrodinger import structure
from schrodinger.structutils import interactions
from tqdm import tqdm  # Import tqdm for the progress bar

def save_ligand_if_interaction_with_residue(protein, maegz_file, output_maegz, score_csv):
    """
    Iterate through ligands in the .maegz file and save the ones that have a salt bridge with
    residue A:218 in the protein to a new .maegz file. Save the docking score for each ligand.

    Parameters:
    protein (Structure object): The protein structure to analyze.
    maegz_file (str): Path to the .maegz file containing ligands.
    output_maegz (str): Path to the output .maegz file where selected ligands will be saved.
    score_csv (str): Path to the CSV file where docking scores will be saved.
    """
    # Get the total number of ligands in the .maegz file excluding the first entry (protein)
    num_ligands = sum(1 for _ in structure.StructureReader(maegz_file)) - 1

    with structure.StructureReader(maegz_file) as reader, structure.StructureWriter(output_maegz) as writer, open(score_csv, 'w', newline='') as csvfile:
        # Setup CSV writer
        csv_writer = csv.writer(csvfile)
        csv_writer.writerow(['Ligand Title', 'r_i_docking_score'])  # CSV header

        # Skip the first entry (protein)
        next(reader)

        # Wrap the loop with tqdm to show progress
        for ligand in tqdm(reader, total=num_ligands, desc="Filtering ligands"):
            # Get the ligand title and docking score
            title = ligand.title
            docking_score = ligand.property.get('r_i_docking_score', 'N/A')

            # Iterate over salt bridges between the protein and the ligand
            for atom1, atom2 in interactions.get_salt_bridges(protein, struc2=ligand):
                # Check if atom1 is part of residue A:218
                if str(atom1.getResidue()) == 'A:218':
                    # Append ligand to the new .maegz file
                    writer.append(ligand)

                    # Save docking score to CSV
                    csv_writer.writerow([title, docking_score])
                    print(f"Saved ligand {title} interacting with A:218 to {output_maegz} with docking score {docking_score}")
                    break  # Save each ligand only once if an interaction is found


def convert_maegz_to_pdb(maegz_file, output_dir):
    """
    Convert ligands from a .maegz file to individual .pdb files, placing no more than 1000 .pdb files
    in each subdirectory.

    Parameters:
    maegz_file (str): Path to the input .maegz file.
    output_dir (str): Directory where .pdb files will be saved.
    """
    # Create the output directory if it does not exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # Get the total number of ligands in the .maegz file excluding the first entry (protein)
    num_ligands = sum(1 for _ in structure.StructureReader(maegz_file)) - 1

    # Load the .maegz file and convert each ligand to .pdb
    with structure.StructureReader(maegz_file) as reader:
        # Skip the first entry (protein)
        next(reader)

        index = 1
        subdirectory_index = 1
        subdirectory = os.path.join(output_dir, f'split_{subdirectory_index}')
        if not os.path.exists(subdirectory):
            os.makedirs(subdirectory)

        # Wrap the loop with tqdm to show progress
        for ligand in tqdm(reader, total=num_ligands, desc="Converting ligands to PDB"):
            # Construct the output filename
            pdb_filename = os.path.join(subdirectory, f'ligand_{index}.pdb')
            index += 1

            # Write the ligand to a .pdb file
            with structure.StructureWriter(pdb_filename) as writer:
                writer.append(ligand)

            # Create a new subdirectory if the limit of 1000 files is reached
            if index % 1000 == 1:
                subdirectory_index += 1
                subdirectory = os.path.join(output_dir, f'split_{subdirectory_index}')
                if not os.path.exists(subdirectory):
                    os.makedirs(subdirectory)


if __name__ == "__main__":
    # Define the input .maegz file and output directories
    maegz_file = '../glide-dock_SP_pv.maegz'
    output_dir = '/home/mepstein/integrin_AI/docking/glide-dock_SP_known_actives/test_interaction_filter'
    output_maegz = '/home/mepstein/integrin_AI/docking/glide-dock_SP_known_actives/test_interaction_filter/ligands_interacting_with_D218_2.maegz'
    score_csv = '/home/mepstein/integrin_AI/docking/glide-dock_SP_known_actives/test_interaction_filter/docking_scores.csv'

    # Read the first entry (protein) separately
    with structure.StructureReader(maegz_file) as reader:
        protein = next(reader)  # Load protein structure

    # Save ligands interacting with residue A:218 and their docking scores
    save_ligand_if_interaction_with_residue(protein, maegz_file, output_maegz, score_csv)

    # Convert ligands to PDB files
    convert_maegz_to_pdb(maegz_file, output_dir)
