import pandas as pd
import pubchempy as pcp
import time
import os

# ================= CONFIGURATION =================
# Input: chemical Name and CAS
INPUT_FILE = "chemicals_list.csv" # see example file format and requirement (2 col csv file)
# Output: CSV ready for 3D conversion script
OUTPUT_FILE = "chemicals_prepared.csv"
# =================================================

def get_clean_smiles(query):
    """
    Fetch SMILES from PubChem and perform auto-desalting/filtering.
    """
    try:
        # Search PubChem by Name or CAS
        compounds = pcp.get_compounds(query, 'name')
        if not compounds:
            return None
        
        # Get the first match
        smiles = compounds[0].isomeric_smiles
        
        # --- 1. Auto-Desalting ---
        # If a dot is present, it's a salt or mixture.
        if "." in smiles:
            parts = smiles.split(".")
            # Keep parts longer than 3 chars (filters out ions like Na+, Cl-, Ca2+)
            valid_parts = [p for p in parts if len(p) > 3]
            if not valid_parts:
                return None
            # Select the longest part as the main organic ligand
            smiles = max(valid_parts, key=len)
            
        # --- 2. Size Filtering ---
        # Skip overly large molecules (Polymers/Proteins) that cause 3D scripts to hang.
        if len(smiles) > 350:
            return "TOO_LARGE"
            
        return smiles
    except Exception:
        return None

def main():
    if not os.path.exists(INPUT_FILE):
        print(f"Error: {INPUT_FILE} not found.")
        return

    print(f"Reading input file: {INPUT_FILE}")
    try:
        # Read the two-column CSV
        df = pd.read_csv(INPUT_FILE)
    except Exception as e:
        print(f"Read Error: {e}")
        return

    processed_list = []
    total_rows = len(df)

    print(f"Starting processing for {total_rows} entries...")

    for index, row in df.iterrows():
        # Data Cleaning: remove surrounding whitespace
        name = str(row.get('Name', '')).strip()
        cas = str(row.get('CAS', '')).strip()
        
        if not name or name.lower() == 'nan':
            continue

        print(f"[{index+1}/{total_rows}] Processing: {name[:30]}...", end=" ", flush=True)

        # Search Strategy: Try CAS first (more specific), then Name
        found_smiles = None
        for identifier in [cas, name]:
            if identifier and identifier.lower() != 'nan' and len(identifier) > 2:
                found_smiles = get_clean_smiles(identifier)
                if found_smiles:
                    break
        
        if found_smiles == "TOO_LARGE":
            print("-> SKIPPED (Polymer/Large Molecule)")
        elif found_smiles:
            print("-> SUCCESS")
            processed_list.append({"Name": name, "SMILES": found_smiles})
        else:
            print("-> FAILED")

        time.sleep(0.2)

    # Save to CSV
    output_df = pd.DataFrame(processed_list)
    output_df.to_csv(OUTPUT_FILE, index=False)
    print(f"\nCompleted! {len(output_df)} compounds saved to {OUTPUT_FILE}.")

if __name__ == "__main__":
    main()