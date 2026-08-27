#!/usr/bin/env python3
"""
Step 4 - Turn the docking matrix into a ranked, *specificity-aware* shortlist.

Input : docking_results.csv (from step 3)
Output: druggable.csv   - the filtered hit list
        docking_matrix.csv - full chemical x protein score matrix (wide form)

The idea: a raw docking score says little on its own, because it scales with
ligand size. What is informative is whether a chemical scores unusually well
against *one* protein compared with (a) its own average over all proteins
(row z-score) and (b) how that protein scores every other chemical
(column z-score). A real specific binder should stand out on both axes.

Read docs/CAVEATS.md before using any of this as evidence. In particular this
script cannot tell you that a hit is real - only that it is unusual.
"""

from __future__ import annotations

import argparse
import sys

import numpy as np
import pandas as pd

SCORE = "Score(kcal/mol)"


def load(path: str, min_proteins: int) -> tuple[pd.DataFrame, dict]:
    df = pd.read_csv(path)
    for col in ("Chemical", "Protein", SCORE):
        if col not in df.columns:
            sys.exit(f"'{path}' must contain a '{col}' column. "
                     f"Found: {list(df.columns)}")

    audit = {"rows_in": len(df)}

    df[SCORE] = pd.to_numeric(df[SCORE], errors="coerce")
    df = df.dropna(subset=[SCORE])
    audit["rows_with_score"] = len(df)

    # Positive scores are steric clashes, not weak binding - they are failures.
    df = df[df[SCORE] < 0]
    audit["rows_negative"] = len(df)

    # One score per pair (keep the best if step 3 was re-run).
    df = df.sort_values(SCORE).drop_duplicates(["Chemical", "Protein"], keep="first")
    audit["unique_pairs"] = len(df)

    counts = df.groupby("Chemical")["Protein"].nunique()
    keep = counts[counts >= min_proteins].index
    audit["chemicals_in"] = int(counts.size)
    audit["chemicals_with_enough_proteins"] = int(len(keep))
    audit["min_proteins"] = min_proteins

    return df[df["Chemical"].isin(keep)].copy(), audit


def add_zscores(df: pd.DataFrame) -> pd.DataFrame:
    """Row-wise (per chemical) and column-wise (per protein) z-scores."""
    g_chem = df.groupby("Chemical")[SCORE]
    df["Chem_Mean"] = g_chem.transform("mean")
    df["Chem_Std"] = g_chem.transform(lambda s: s.std(ddof=1))
    df["Z_Chem"] = (df[SCORE] - df["Chem_Mean"]) / df["Chem_Std"].replace(0, np.nan)

    g_prot = df.groupby("Protein")[SCORE]
    df["Prot_Mean"] = g_prot.transform("mean")
    df["Prot_Std"] = g_prot.transform(lambda s: s.std(ddof=1))
    df["Z_Prot"] = (df[SCORE] - df["Prot_Mean"]) / df["Prot_Std"].replace(0, np.nan)
    return df


def add_ligand_efficiency(df: pd.DataFrame, index_path: str | None) -> pd.DataFrame:
    """LE = score / heavy atoms. The single best guard against size bias."""
    if not index_path:
        return df
    try:
        idx = pd.read_csv(index_path)
    except Exception as exc:
        print(f"[warn] could not read {index_path}: {exc}")
        return df
    if "SafeName" not in idx.columns or "SMILES" not in idx.columns:
        print(f"[warn] {index_path} has no SafeName/SMILES columns - skipping LE")
        return df

    try:
        from rdkit import Chem, RDLogger
        RDLogger.DisableLog("rdApp.*")
    except ImportError:
        print("[warn] RDKit not installed - ligand efficiency not computed "
              "(pip install rdkit)")
        return df

    heavy = {}
    for _, r in idx.iterrows():
        mol = Chem.MolFromSmiles(str(r["SMILES"]))
        if mol is not None:
            heavy[str(r["SafeName"])] = mol.GetNumHeavyAtoms()

    df["HeavyAtoms"] = df["Chemical"].map(heavy)
    df["LigandEfficiency"] = df[SCORE] / df["HeavyAtoms"]
    return df


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-i", "--input", default="docking_results.csv")
    p.add_argument("-o", "--output", default="druggable.csv")
    p.add_argument("--matrix", default="docking_matrix.csv",
                   help="wide chemical x protein score table (default: %(default)s)")
    p.add_argument("--ligand-index", default="ligands/ligand_index.csv",
                   help="step-2 index, used for ligand efficiency; '' to skip")
    p.add_argument("--z-threshold", type=float, default=-1.5,
                   help="per-chemical z-score a hit must beat (default: %(default)s)")
    p.add_argument("--protein-z-threshold", type=float, default=None,
                   help="optionally also require this per-protein z-score "
                        "(e.g. -1.5). Off by default.")
    p.add_argument("--std-threshold", type=float, default=1.5,
                   help="minimum spread of a chemical's scores across proteins; "
                        "flat profiles are promiscuous binders (default: %(default)s)")
    p.add_argument("--score-lead-min", type=float, default=-8.0,
                   help="a chemical must reach at least this score somewhere "
                        "(default: %(default)s). SIZE-BIASED - see docs/CAVEATS.md.")
    p.add_argument("--min-proteins", type=int, default=5,
                   help="minimum proteins per chemical before a z-score is "
                        "meaningful (default: %(default)s)")
    p.add_argument("--no-lead-filter", action="store_true",
                   help="skip the absolute score cutoff and rank on specificity only")
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)

    df, audit = load(args.input, args.min_proteins)
    if df.empty:
        sys.exit(f"No chemical has scores against >= {args.min_proteins} proteins. "
                 f"Lower --min-proteins or dock more targets.")

    # Wide matrix - useful for clustering / heatmaps / the supplement.
    df.pivot_table(index="Chemical", columns="Protein", values=SCORE) \
      .to_csv(args.matrix)

    if not args.no_lead_filter:
        leads = df[df[SCORE] <= args.score_lead_min]["Chemical"].unique()
        df = df[df["Chemical"].isin(leads)].copy()
        audit["chemicals_passing_lead_filter"] = int(len(leads))

    if df.empty:
        sys.exit(f"No chemical reached {args.score_lead_min} kcal/mol. "
                 f"Try --no-lead-filter or a softer --score-lead-min.")

    df = add_zscores(df)
    df = add_ligand_efficiency(df, args.ligand_index or None)

    mask = (df["Chem_Std"] > args.std_threshold) & (df["Z_Chem"] <= args.z_threshold)
    if args.protein_z_threshold is not None:
        mask &= df["Z_Prot"] <= args.protein_z_threshold

    report = df[mask].sort_values(["Z_Chem", SCORE]).copy()

    cols = ["Chemical", "Protein", SCORE, "Z_Chem", "Z_Prot",
            "Chem_Mean", "Chem_Std"]
    if "LigandEfficiency" in report.columns:
        cols += ["HeavyAtoms", "LigandEfficiency"]
    report = report[cols]
    report.to_csv(args.output, index=False)

    # ---------------- report ----------------
    print("\nFiltering audit")
    print("-" * 62)
    for k, v in audit.items():
        print(f"  {k:<36} {v}")
    print(f"  {'hits after specificity filter':<36} {len(report)}"
          f" ({report['Chemical'].nunique()} chemicals)")

    if report.empty:
        print("\nNo hits survived. Loosen --z-threshold / --std-threshold, "
              "or dock more proteins.")
        return

    header = f"{'CHEMICAL':<26} | {'TARGET':<20} | {'SCORE':>7} | {'Z_chem':>7} | {'Z_prot':>7}"
    if "LigandEfficiency" in report.columns:
        header += f" | {'LE':>6}"
    print("\n" + "=" * len(header))
    print(header)
    print("-" * len(header))
    for _, r in report.iterrows():
        line = (f"{str(r['Chemical'])[:26]:<26} | {str(r['Protein'])[:20]:<20} | "
                f"{r[SCORE]:>7.1f} | {r['Z_Chem']:>7.2f} | "
                f"{r['Z_Prot'] if pd.notna(r['Z_Prot']) else float('nan'):>7.2f}")
        if "LigandEfficiency" in report.columns:
            le = r["LigandEfficiency"]
            line += f" | {le:>6.2f}" if pd.notna(le) else f" | {'n/a':>6}"
        print(line)
    print("=" * len(header))
    print(f"\nHits   -> {args.output}")
    print(f"Matrix -> {args.matrix}")
    print("\nReminder: these are hypotheses to test, not evidence of binding. "
          "See docs/CAVEATS.md.")


if __name__ == "__main__":
    main()
