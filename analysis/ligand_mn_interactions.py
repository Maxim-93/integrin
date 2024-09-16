import os
import glob
import MDAnalysis as mda
from MDAnalysis.analysis.distances import distance_array

# Define the paths
protein_path = '/home/mepstein/integrin_AI/docking/MMGBSA/filter_based_on_interaction/6MK0_prepped_no_ligand.pdb'

# Use glob to find all split directories (i.e., directories that match 'split_*')
base_dir = '/home/mepstein/integrin_AI/docking/glide-dock_SP_known_actives/test_interaction_filter'
ligand_dirs = glob.glob(os.path.join(base_dir, 'split_*'))

# Load the protein
u_protein = mda.Universe(protein_path)

# Select Mn ion from the protein (resid 708)
mn_ion = u_protein.select_atoms('resname MN and resid 708')

# Define a cutoff distance (e.g., 5 Å for close proximity)
cutoff_distance = 5.0

# Open the output file to write the names of interacting ligands
with open('mn_interacting_ligands.txt', 'w') as output_file:
    # Loop over the found split directories
    for ligand_dir in ligand_dirs:
        # Find the ligand .pdb file in the current subdirectory
        ligand_file = [f for f in os.listdir(ligand_dir) if f.endswith('.pdb')][0]
        ligand_path = os.path.join(ligand_dir, ligand_file)
        
        # Load the ligand
        u_ligand = mda.Universe(ligand_path)
        
        # Select ligand atoms
        ligand = u_ligand.select_atoms('all')
        
        # Compute distance between all ligand atoms and Mn ion
        distances = distance_array(ligand.positions, mn_ion.positions)
        
        # Check if any distance is less than the cutoff
        close_contacts = (distances < cutoff_distance)
        
        # If ligand is within the cutoff distance, write its name to the file
        if close_contacts.any():
            output_file.write(f"{ligand_file}\n")
            print(f"Ligand {ligand_file} is in close proximity to the Mn ion within {cutoff_distance} Å")
        else:
            print(f"Ligand {ligand_file} is not in close proximity to the Mn ion within {cutoff_distance} Å")
