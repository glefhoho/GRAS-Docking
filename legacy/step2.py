import pandas as pd
import subprocess
import os
import sys

# ================= CONFIGURATION =================
INPUT_FILE = "chemicals_prepared.csv" # the output from step 1
OUTPUT_FOLDER = "ligands" # chemical structures in this folder
PH_LEVEL = 6.5 # change this value based on your research needs
TIME_LIMIT = "20s"  # Kills the process if it takes > 20 seconds
# =================================================

def main():
    if not os.path.exists(OUTPUT_FOLDER):
        os.makedirs(OUTPUT_FOLDER)
    
    try:
        df = pd.read_csv(INPUT_FILE)
    except Exception as e:
        print("Error reading file: {}".format(e))
        return

    print("Starting conversion with a {} timeout per chemical...".format(TIME_LIMIT))
    print("-" * 40)

    success_count = 0
    fail_count = 0
    timeout_count = 0

    for index, row in df.iterrows():
        # Clean Data
        name = str(row['Name']).strip()
        smiles = str(row['SMILES']).strip()
        safe_name = name.replace(" ", "_").replace("/", "-").replace("(", "").replace(")", "").replace(",", "")
        output_file = os.path.join(OUTPUT_FOLDER, "{}.pdbqt".format(safe_name))

        # Check if already done (Optional: Remove this if you want to overwrite)
        # if os.path.exists(output_file) and os.path.getsize(output_file) > 0:
        #    print("[{}] Skipping (Already exists): {}".format(index+1, safe_name))
        #    continue

        # THE MAGIC COMMAND: "timeout 20s obabel ..."
        # This forces Linux to kill OpenBabel if it hangs.
        cmd = 'timeout {} obabel -:"{}" -O "{}" --gen3d -p {}'.format(TIME_LIMIT, smiles, output_file, PH_LEVEL)

        # Print BEFORE running so you know exactly which one is working
        sys.stdout.write("[{}] Processing: {} ... ".format(index+1, safe_name))
        sys.stdout.flush()

        try:
            # Run command
            retcode = subprocess.call(cmd, shell=True)
            
            # Check Result
            if retcode == 0:
                print("Done.")
                success_count += 1
            elif retcode == 124: # 124 is the Linux exit code for "Timed Out"
                print("SKIPPED (Timed Out).")
                # Remove empty file if created
                if os.path.exists(output_file):
                    os.remove(output_file)
                timeout_count += 1
            else:
                print("FAILED (Error Code {}).".format(retcode))
                fail_count += 1
                
        except Exception as e:
            print("Error: {}".format(e))

    print("-" * 40)
    print("Finished.")
    print("Success: {}".format(success_count))
    print("Timed Out: {} (Too complex for 3D)".format(timeout_count))
    print("Failed: {}".format(fail_count))

if __name__ == "__main__":
    main()