#!/usr/bin/env python3
"""
Step 2 - Turn SMILES into 3D, protonated, docking-ready ligand files.

Input : chemicals_prepared.csv (from step 1) with columns Name,SMILES
Output: ligands/<safe_name>.pdbqt  +  ligands/ligand_index.csv

`ligand_index.csv` maps the sanitized file name back to the original chemical
name, SMILES and CAS. Steps 3 and 4 only ever see the sanitized name, so this
file is what lets you trace a docking hit back to a real compound.

See docs/CHANGES.md for what differs from legacy/step2.py.
"""

from __future__ import annotations

import argparse
import os
import re
import shutil
import subprocess
import sys

import pandas as pd


def safe_filename(name: str) -> str:
    """ASCII-safe, shell-safe file stem."""
    stem = re.sub(r"[^A-Za-z0-9._-]+", "_", str(name).strip())
    return stem.strip("._-") or "unnamed"


def unique_names(names) -> list[str]:
    """Sanitize names and disambiguate collisions with _2, _3, ..."""
    seen: dict[str, int] = {}
    out = []
    for n in names:
        stem = safe_filename(n)
        seen[stem] = seen.get(stem, 0) + 1
        out.append(stem if seen[stem] == 1 else f"{stem}_{seen[stem]}")
    return out


def pdbqt_is_valid(path: str, min_atoms: int = 3) -> bool:
    """Open Babel can exit 0 and still leave an empty or 1-atom file behind."""
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return False
    n = 0
    with open(path, "r", errors="ignore") as fh:
        for line in fh:
            if line.startswith(("ATOM", "HETATM")):
                n += 1
                if n >= min_atoms:
                    return True
    return False


def convert(smiles: str, out_path: str, ph: float, gen3d: str, timeout: int,
            forcefield: str | None, steps: int) -> tuple[bool, str]:
    """Run obabel for one molecule. Returns (ok, status)."""
    cmd = ["obabel", f"-:{smiles}", "-O", out_path, "--gen3d", gen3d, "-p", str(ph)]
    if forcefield:
        cmd += ["--minimize", "--ff", forcefield, "--steps", str(steps)]
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        if os.path.exists(out_path):
            os.remove(out_path)
        return False, "timeout"

    if proc.returncode != 0:
        if os.path.exists(out_path) and not pdbqt_is_valid(out_path):
            os.remove(out_path)
        err = (proc.stderr or "").strip().splitlines()
        return False, f"obabel_rc{proc.returncode}: {err[-1][:60] if err else ''}"

    if not pdbqt_is_valid(out_path):
        if os.path.exists(out_path):
            os.remove(out_path)
        return False, "empty_output"

    return True, "ok"


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-i", "--input", default="chemicals_prepared.csv",
                   help="output of step 1 (default: %(default)s)")
    p.add_argument("-o", "--outdir", default="ligands",
                   help="folder for .pdbqt files (default: %(default)s)")
    p.add_argument("--ph", type=float, default=6.5,
                   help="pH used for protonation; set it to the pH of your system "
                        "(default: %(default)s)")
    p.add_argument("--gen3d", default="med", choices=["fastest", "fast", "med",
                                                      "better", "best"],
                   help="Open Babel 3D generation effort (default: %(default)s)")
    p.add_argument("--minimize", metavar="FF", default=None,
                   choices=["MMFF94", "MMFF94s", "GAFF", "UFF", "Ghemical"],
                   help="optionally relax the generated conformer with this force field")
    p.add_argument("--steps", type=int, default=500,
                   help="minimization steps when --minimize is used (default: %(default)s)")
    p.add_argument("--timeout", type=int, default=60,
                   help="seconds allowed per molecule (default: %(default)s)")
    p.add_argument("--overwrite", action="store_true",
                   help="rebuild ligands that already exist")
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)

    if shutil.which("obabel") is None:
        sys.exit("Open Babel ('obabel') not found on PATH.\n"
                 "  conda install -c conda-forge openbabel\n"
                 "  or: sudo apt-get install openbabel")

    if not os.path.exists(args.input):
        sys.exit(f"Input file not found: {args.input}")

    df = pd.read_csv(args.input)
    for col in ("Name", "SMILES"):
        if col not in df.columns:
            sys.exit(f"Input must contain a '{col}' column. Found: {list(df.columns)}")
    df = df[df["SMILES"].notna()].reset_index(drop=True)

    os.makedirs(args.outdir, exist_ok=True)
    df["SafeName"] = unique_names(df["Name"])

    index_rows = []
    n_ok = n_skip = n_timeout = n_fail = 0
    total = len(df)

    print(f"Building {total} ligands at pH {args.ph} "
          f"(gen3d={args.gen3d}, timeout={args.timeout}s)")
    print("-" * 62)

    for i, row in df.iterrows():
        stem = row["SafeName"]
        smiles = str(row["SMILES"]).strip()
        out_path = os.path.join(args.outdir, f"{stem}.pdbqt")

        if pdbqt_is_valid(out_path) and not args.overwrite:
            print(f"[{i + 1}/{total}] {stem[:36]:<36} exists")
            status, ok, n_skip = "exists", True, n_skip + 1
        else:
            print(f"[{i + 1}/{total}] {stem[:36]:<36}", end=" ", flush=True)
            ok, status = convert(smiles, out_path, args.ph, args.gen3d,
                                 args.timeout, args.minimize, args.steps)
            print(status)
            if ok:
                n_ok += 1
            elif status == "timeout":
                n_timeout += 1
            else:
                n_fail += 1

        index_rows.append({
            "SafeName": stem,
            "Name": row["Name"],
            "CAS": row.get("CAS", ""),
            "CID": row.get("CID", ""),
            "SMILES": smiles,
            "MolecularWeight": row.get("MolecularWeight", ""),
            "File": out_path if ok else "",
            "Status": status,
        })

    index_path = os.path.join(args.outdir, "ligand_index.csv")
    pd.DataFrame(index_rows).to_csv(index_path, index=False)

    print("-" * 62)
    print(f"Built    : {n_ok}")
    print(f"Reused   : {n_skip}")
    print(f"Timed out: {n_timeout}  (3D embedding too slow - usually flexible polymers)")
    print(f"Failed   : {n_fail}")
    print(f"Index    : {index_path}")


if __name__ == "__main__":
    main()
