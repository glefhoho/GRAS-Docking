import os
import subprocess
import glob

# ================= CONFIGURATION =================
PROTEINS_DIR = "test_pdb" # the pdbs of the proteins; you can just run a few for your test     
LIGANDS_DIR = "ligands"   # output from step 2 script 
OUTPUT_FILE = "docking_results.csv"
# =================================================

def get_pocket_center(pdb_file):
    """
    Reads the pocket center from the existing _out folders you created.
    """
    base_name = os.path.basename(pdb_file) # e.g. pdb3sk3_ButyrateKinase.pdb
    root_name = os.path.splitext(base_name)[0] # e.g. pdb3sk3_ButyrateKinase
    
    # Based on your screenshot, the folder is named 'filename_out'
    pocket_dir = os.path.join(os.path.dirname(pdb_file), f"{root_name}_out")
    pocket_file = os.path.join(pocket_dir, "pockets", "pocket1_atm.pdb")
    
    if not os.path.exists(pocket_file):
        print(f"    [WARNING] Could not find pocket file: {pocket_file}")
        return None

    # Calculate Center
    coords = []
    with open(pocket_file, 'r') as f:
        for line in f:
            if line.startswith("ATOM"):
                try:
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])
                    coords.append((x, y, z))
                except: continue
    
    if not coords:
        return None

    avg_x = sum(c[0] for c in coords) / len(coords)
    avg_y = sum(c[1] for c in coords) / len(coords)
    avg_z = sum(c[2] for c in coords) / len(coords)
    
    return (avg_x, avg_y, avg_z)

def main():
    # 1. Setup
    print(f"--- Starting Docking Pilot ---")
    
    # Check if ligands folder exists
    if not os.path.exists(LIGANDS_DIR):
        print(f"[ERROR] '{LIGANDS_DIR}' folder is missing! Please create it and put your chemical .pdbqt files inside.")
        return

    proteins = glob.glob(f"{PROTEINS_DIR}/*.pdb")
    ligands = glob.glob(f"{LIGANDS_DIR}/*.pdbqt")

    print(f"Found {len(proteins)} Proteins and {len(ligands)} Ligands.")
    
    if len(ligands) == 0:
        print("[ERROR] No .pdbqt files found in ligands folder!")
        return

    # Prepare Output
    with open(OUTPUT_FILE, "w") as f:
        f.write("Chemical,Protein,Score(kcal/mol)\n")

    # 2. Loop through proteins
    for prot_pdb in proteins:
        prot_name = os.path.basename(prot_pdb).replace(".pdb", "")
        print(f"\nProcessing: {prot_name}")

        # A. Find Pocket Center (Using your existing output)
        center = get_pocket_center(prot_pdb)
        if not center:
            print("  Skipping (No pocket data found). Did you run fpocket?")
            continue

        # B. Convert Protein to PDBQT (Crucial step missing in screenshot)
        prot_pdbqt = prot_pdb.replace(".pdb", ".pdbqt")
        if not os.path.exists(prot_pdbqt):
            print("  Converting PDB to PDBQT...")
            subprocess.run(f"obabel {prot_pdb} -O {prot_pdbqt} -xr -p 6.5", shell=True)

        # C. Dock every ligand
        for lig in ligands:
            lig_name = os.path.basename(lig).replace(".pdbqt", "")
            
            # Run Smina
            # Note: CPU set to 4. Increase if you have a powerful server.
            cmd = (f"smina -r {prot_pdbqt} -l {lig} "
                   f"--center_x {center[0]:.3f} --center_y {center[1]:.3f} --center_z {center[2]:.3f} "
                   f"--size_x 20 --size_y 20 --size_z 20 "
                   f"--num_modes 1 --cpu 4 --out results/{prot_name}_{lig_name}.pdbqt")
            
            # Ensure results dir exists
            os.makedirs("results", exist_ok=True)
            
            # Execute and Capture
            result = subprocess.run(cmd, shell=True, capture_output=True, text=True)
            
            # Parse Score
            score = 100.0
            for line in result.stdout.splitlines():
                if line.strip().startswith("1"):
                    try:
                        score = float(line.split()[1])
                    except: pass
            
            print(f"  > Docked vs {lig_name}: {score}")
            with open(OUTPUT_FILE, "a") as f:
                f.write(f"{lig_name},{prot_name},{score}\n")

if __name__ == "__main__":
    main()