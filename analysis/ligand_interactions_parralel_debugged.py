import os
import MDAnalysis as mda
from tqdm import tqdm
import numpy as np
from concurrent.futures import ThreadPoolExecutor, as_completed
import re

# Paths
protein_pdb = '/home/mepstein/integrin_AI/docking/MMGBSA/filter_based_on_interaction/6MK0_prepped_no_ligand.pdb'
ligand_dir = '/home/mepstein/integrin_AI/docking/glide-dock_XP_filtered_actives/ligands_converted'
interaction_results_dir = '/home/mepstein/integrin_AI/docking/glide-dock_XP_filtered_actives/interaction_results'
filtered_ligands_file = '/home/mepstein/integrin_AI/docking/glide-dock_XP_filtered_actives/filtered_ligands_2.txt'

# Create a directory for interaction results if it doesn't exist
if not os.path.exists(interaction_results_dir):
    os.makedirs(interaction_results_dir)


# Function to perform interaction analysis
def perform_interaction_analysis(ligand_path):
    asp_contact = []
    MIDAS_contact = []
    ligand_pdb = os.path.basename(ligand_path)
    
    # Extract numeric part from ligand PDB filename
    match = re.search(r'(\d+)', ligand_pdb)
    if match:
        ligand_number = int(match.group(1))
        adjusted_number = ligand_number - 1
    else:
        adjusted_number = None  # Handle cases where numeric part is not found

    # Load protein and ligand into MDAnalysis Universe
    protein = mda.Universe(protein_pdb)
    ligand = mda.Universe(ligand_path)

    # Merge protein and ligand into a single universe
    merged = mda.Merge(protein.atoms, ligand.atoms)

    # Select the specific residues and atoms
    MIDAS = merged.select_atoms("resname MN and resid 708")
    Asp218 = merged.select_atoms("resname ASP and resid 218 and (name OD1 or name OD2)")
    Asp218_coords = Asp218.positions
    MIDAS_coords = MIDAS.positions

    # Ligand selection
    electronegative_ligand_atoms = merged.select_atoms(
        "resname UNK and (name O* or name N* or name F* or name Cl* or name Br*)")

    # Asp contact selection
    non_carbon_atoms = merged.select_atoms('resname UNK and not name C*')
    ligand_hydrogens = merged.select_atoms('resname UNK and name H*')

    polar_hydrogens = []
    for non_carbon in non_carbon_atoms:
        for hydrogen in ligand_hydrogens:
            distance = np.linalg.norm(non_carbon.position - hydrogen.position)
            if distance < 1.05:
                polar_hydrogens.append(hydrogen)

    asp_distance=[]
    for asp_oxygen in Asp218_coords:
        for hydrogen in polar_hydrogens:
            if np.linalg.norm(hydrogen.position - asp_oxygen) < 3.5:
                asp_distance.append((np.linalg.norm(hydrogen.position - asp_oxygen)))
                asp_contact.append(1)

    mn_distance=[]
    # MIDAS contact selection
    for atom in electronegative_ligand_atoms:
        if np.linalg.norm(atom.position - MIDAS_coords[0]) < 3.5:
            mn_distance.append(np.linalg.norm(atom.position - MIDAS_coords[0]))
            MIDAS_contact.append(1)

    if sum(asp_contact) > 0 and sum(MIDAS_contact) > 0:
        print('distances are as follows:')
        print('asp distance is:'+str(asp_distance))
        print('mn distance is'+str(mn_distance))
        with open(filtered_ligands_file, 'a') as f:
            if adjusted_number is not None:
                f.write(f"{adjusted_number}\n")
            else:
                print(f"Warning: Numeric part not found in {ligand_pdb}")
        return adjusted_number

    return None

# Collect all ligand PDB files
ligand_files = []
for subdir in sorted(os.listdir(ligand_dir)):
    subdir_path = os.path.join(ligand_dir, subdir)
    if os.path.isdir(subdir_path):
        for ligand_pdb in sorted(os.listdir(subdir_path)):
            if ligand_pdb.endswith('.pdb'):
                ligand_path = os.path.join(subdir_path, ligand_pdb)
                ligand_files.append(ligand_path)

# Parallel processing using ThreadPoolExecutor
with ThreadPoolExecutor(max_workers=os.cpu_count()) as executor:
    futures = [executor.submit(perform_interaction_analysis, ligand_path) for ligand_path in ligand_files]

    for future in tqdm(as_completed(futures), total=len(futures), desc="Processing Ligands"):
        ligand_number = future.result()
        if ligand_number is not None:
            print(f"Ligand number {ligand_number} has been added to {filtered_ligands_file}")

print('All ligands have been analyzed.')
