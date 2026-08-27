import pandas as pd

# ================= CONFIGURATION =================
INPUT_FILE = "docking_results.csv"
OUTPUT_FILE = "druggable.csv"

# Thresholds:
Z_THRESHOLD = -1.5     # Must be 1.5 standard deviations better than average
STD_THRESHOLD = 1.5     # Must show significant variation across proteins
SCORE_LEAD_MIN = -8.0   # Must hit AT LEAST ONE target with strong affinity
# =================================================

def main():
    print(f"Loading results from {INPUT_FILE}...")
    try:
        df = pd.read_csv(INPUT_FILE)
    except Exception as e:
        print(f"Error: {e}")
        return

    # 1. Keep chemicals that have at least one strong hit (<= -8)
    # This ensures we are only looking at "real drugs", not junk molecules.
    strong_chemicals = df[df['Score(kcal/mol)'] <= SCORE_LEAD_MIN]['Chemical'].unique()
    df_filtered = df[df['Chemical'].isin(strong_chemicals)].copy()

    # 2. Data Cleaning
    df_filtered = df_filtered[df_filtered['Score(kcal/mol)'] < 0] # Remove positive clashes

    # 3. Horizontal Statistics (Per Chemical)
    stats = df_filtered.groupby('Chemical')['Score(kcal/mol)'].agg(['mean', 'std']).reset_index()
    stats.columns = ['Chemical', 'Chem_Mean', 'Chem_Std']

    # 4. Merge and Calculate Z-Score
    df_final = pd.merge(df_filtered, stats, on='Chemical')
    df_final['Z_Score'] = (df_final['Score(kcal/mol)'] - df_final['Chem_Mean']) / df_final['Chem_Std']

    # 5. FINAL FILTERING LOGIC
    # Rule 1: High Variance (Std > 1.5) -> Kills Dextrins (sticky but flat profile)
    # Rule 2: High Specificity (Z <= -1.5) -> Only shows the targets it truly "likes"
    final_report = df_final[
        (df_final['Chem_Std'] > STD_THRESHOLD) & 
        (df_final['Z_Score'] <= Z_THRESHOLD)
    ].copy()

    # Sort for clarity
    final_report = final_report.sort_values(by=['Chemical', 'Z_Score'])

    # 6. Display Results
    print("\n" + "="*85)
    print(f"{'CHEMICAL':<25} | {'TARGET PDB':<12} | {'SCORE':<8} | {'Z-SCORE':<8} | {'STD'}")
    print("-" * 85)
    
    for chemical in final_report['Chemical'].unique():
        data = final_report[final_report['Chemical'] == chemical]
        std_val = data['Chem_Std'].iloc[0]
        for _, row in data.iterrows():
            print(f"{row['Chemical']:<25} | {row['Protein']:<12} | {row['Score(kcal/mol)']:<8.1f} | {row['Z_Score']:<8.2f} | {std_val:.2f}")
    
    final_report.to_csv(OUTPUT_FILE, index=False)
    print("="*85)
    print(f"\nSaved specific profiles to {OUTPUT_FILE}")

if __name__ == "__main__":
    main()