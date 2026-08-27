#!/usr/bin/env python3
"""
Step 1 - Resolve chemical names / CAS numbers to docking-ready SMILES.

Input : 2-column CSV with headers `Name,CAS` (see data/example_chemicals.csv)
Output: `chemicals_prepared.csv`  -> Name,CAS,CID,SMILES,MolecularWeight,Source
        `chemicals_failed.csv`    -> everything that could NOT be resolved, with a reason

Design notes (see docs/CHANGES.md for what differs from legacy/step1.py):
  * SMILES attribute name changed across PubChemPy versions (isomeric_smiles ->
    smiles). We try several, so the script works on 1.0.4 and 1.0.5+.
  * Desalting uses RDKit's largest-fragment rule when RDKit is installed and
    falls back to a heavy-atom heuristic otherwise. Never a raw string length.
  * Failures are written to a file instead of being silently dropped.
  * Resumable: existing output rows are re-used unless --overwrite is given.
"""

from __future__ import annotations

import argparse
import os
import re
import sys
import time
import warnings

import pandas as pd

warnings.filterwarnings("ignore", category=DeprecationWarning)

try:
    import pubchempy as pcp
except ImportError:  # pragma: no cover
    sys.exit("pubchempy is not installed.  ->  pip install pubchempy")

try:
    from rdkit import Chem
    from rdkit import RDLogger
    from rdkit.Chem.Descriptors import MolWt

    RDLogger.DisableLog("rdApp.*")
    HAVE_RDKIT = True
except ImportError:  # pragma: no cover
    HAVE_RDKIT = False


# --------------------------------------------------------------------------
# SMILES retrieval
# --------------------------------------------------------------------------

# PubChemPy renamed these in 1.0.5; older releases only have the first two.
SMILES_ATTRS = ("smiles", "isomeric_smiles", "connectivity_smiles", "canonical_smiles")

# Counterions / small ions we never want to keep as "the ligand".
COMMON_COUNTERIONS = {
    "[Na+]", "[K+]", "[Li+]", "[Cl-]", "[Br-]", "[I-]", "[F-]",
    "[Ca+2]", "[Ca++]", "[Mg+2]", "[Mg++]", "[Zn+2]", "[Fe+2]", "[Fe+3]",
    "[NH4+]", "[OH-]", "O", "[H+]", "[Al+3]", "[Cu+2]", "[Mn+2]",
    "OS(=O)(=O)O", "Cl", "[N+](=O)[O-]",
}

_ATOM_TOKEN = re.compile(r"Cl|Br|Si|Se|[BCNOPSFIbcnops]|\[[^\]]+\]")


def heavy_atom_estimate(smiles: str) -> int:
    """Rough heavy-atom count from a SMILES string (RDKit-free fallback)."""
    n = 0
    for tok in _ATOM_TOKEN.findall(smiles):
        if tok.startswith("["):
            inner = tok[1:-1]
            if re.fullmatch(r"\d*H\d*[+-]?\d*", inner):  # bare hydrogen
                continue
        n += 1
    return n


def smiles_from_compound(compound) -> str | None:
    """Read a SMILES string off a PubChemPy Compound across library versions."""
    for attr in SMILES_ATTRS:
        try:
            value = getattr(compound, attr, None)
        except Exception:
            continue
        if value:
            return str(value)
    # Last resort: ask PUG-REST directly for whichever property still exists.
    for prop in ("SMILES", "IsomericSMILES", "CanonicalSMILES"):
        try:
            rows = pcp.get_properties([prop], compound.cid, "cid")
        except Exception:
            continue
        if rows and rows[0].get(prop):
            return str(rows[0][prop])
    return None


def largest_fragment(smiles: str) -> str | None:
    """Strip salts / mixtures and return the main organic fragment."""
    if "." not in smiles:
        return smiles

    parts = [p for p in smiles.split(".") if p]
    parts = [p for p in parts if p not in COMMON_COUNTERIONS] or parts

    if HAVE_RDKIT:
        best, best_n = None, -1
        for p in parts:
            mol = Chem.MolFromSmiles(p, sanitize=False)
            if mol is None:
                continue
            n = mol.GetNumHeavyAtoms()
            if n > best_n:
                best, best_n = p, n
        if best is not None and best_n >= 2:
            return best
        return None

    best = max(parts, key=heavy_atom_estimate)
    return best if heavy_atom_estimate(best) >= 2 else None


def canonicalize(smiles: str) -> tuple[str | None, float | None]:
    """Validate + canonicalize with RDKit when available. Returns (smiles, MW)."""
    if not HAVE_RDKIT:
        return smiles, None
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, None
    return Chem.MolToSmiles(mol), round(MolWt(mol), 2)


def query_pubchem(query: str, retries: int, pause: float):
    """pcp.get_compounds with retry/backoff on transient PubChem errors."""
    last_err = None
    for attempt in range(retries):
        try:
            return pcp.get_compounds(query, "name"), None
        except Exception as exc:  # network hiccups, 503 throttling, ...
            last_err = exc
            time.sleep(pause * (2 ** attempt))
    return [], last_err


def resolve(query: str, args) -> tuple[str | None, int | None, float | None, str]:
    """Return (smiles, cid, mw, status) for a single identifier."""
    compounds, err = query_pubchem(query, args.retries, args.pause)
    if err is not None:
        return None, None, None, f"pubchem_error: {type(err).__name__}"
    if not compounds:
        return None, None, None, "not_found"

    compound = compounds[0]
    raw = smiles_from_compound(compound)
    if not raw:
        return None, compound.cid, None, "no_smiles_in_record"

    main = largest_fragment(raw)
    if not main:
        return None, compound.cid, None, "only_counterions"

    smiles, mw = canonicalize(main)
    if smiles is None:
        return None, compound.cid, None, "rdkit_parse_failed"

    if mw is None:  # RDKit absent -> fall back to PubChem's MW
        try:
            mw = float(compound.molecular_weight) if compound.molecular_weight else None
        except Exception:
            mw = None

    n_heavy = heavy_atom_estimate(smiles)
    if mw is not None and mw > args.max_mw:
        return None, compound.cid, mw, f"too_large_mw>{args.max_mw}"
    if mw is None and n_heavy > args.max_heavy:
        return None, compound.cid, None, f"too_large_heavy>{args.max_heavy}"
    if n_heavy < args.min_heavy:
        return None, compound.cid, mw, f"too_small_heavy<{args.min_heavy}"

    return smiles, compound.cid, mw, "ok"


# --------------------------------------------------------------------------
# Driver
# --------------------------------------------------------------------------

def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-i", "--input", default="chemicals_list.csv",
                   help="2-column CSV with headers Name,CAS (default: %(default)s)")
    p.add_argument("-o", "--output", default="chemicals_prepared.csv",
                   help="resolved compounds (default: %(default)s)")
    p.add_argument("--failed", default="chemicals_failed.csv",
                   help="unresolved compounds + reason (default: %(default)s)")
    p.add_argument("--max-mw", type=float, default=900.0,
                   help="skip molecules heavier than this, in Da (default: %(default)s). "
                        "Docking scores are unreliable well before this limit.")
    p.add_argument("--max-heavy", type=int, default=70,
                   help="heavy-atom ceiling used when no MW is available (default: %(default)s)")
    p.add_argument("--min-heavy", type=int, default=5,
                   help="skip fragments smaller than this (default: %(default)s)")
    p.add_argument("--pause", type=float, default=0.25,
                   help="seconds between PubChem requests; PubChem asks for <=5/s "
                        "(default: %(default)s)")
    p.add_argument("--retries", type=int, default=3,
                   help="retries per query with exponential backoff (default: %(default)s)")
    p.add_argument("--overwrite", action="store_true",
                   help="ignore an existing output file instead of resuming")
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)

    if not os.path.exists(args.input):
        sys.exit(f"Input file not found: {args.input}")

    df = pd.read_csv(args.input)
    df.columns = [str(c).strip() for c in df.columns]
    if "Name" not in df.columns:
        sys.exit(f"Input must have a 'Name' column. Found: {list(df.columns)}")
    if "CAS" not in df.columns:
        df["CAS"] = ""

    done: dict[str, dict] = {}
    if os.path.exists(args.output) and not args.overwrite:
        prev = pd.read_csv(args.output)
        done = {str(r["Name"]).strip(): r.to_dict() for _, r in prev.iterrows()}
        if done:
            print(f"Resuming: {len(done)} compounds already resolved in {args.output}")

    ok_rows, bad_rows = list(done.values()), []
    total = len(df)

    for idx, row in df.iterrows():
        name = str(row.get("Name", "")).strip()
        cas = str(row.get("CAS", "")).strip()
        if not name or name.lower() == "nan":
            continue
        if name in done:
            continue

        print(f"[{idx + 1}/{total}] {name[:40]:<40}", end=" ", flush=True)

        smiles = cid = mw = None
        status = "not_found"
        # CAS first: it is unambiguous; the trivial name often is not.
        for identifier in (cas, name):
            if not identifier or identifier.lower() == "nan" or len(identifier) < 3:
                continue
            smiles, cid, mw, status = resolve(identifier, args)
            time.sleep(args.pause)
            if smiles:
                break

        if smiles:
            print(f"OK   CID {cid}")
            ok_rows.append({"Name": name, "CAS": cas, "CID": cid, "SMILES": smiles,
                            "MolecularWeight": mw, "Source": "PubChem"})
        else:
            print(f"MISS ({status})")
            bad_rows.append({"Name": name, "CAS": cas, "CID": cid,
                             "MolecularWeight": mw, "Reason": status})

    pd.DataFrame(ok_rows).to_csv(args.output, index=False)
    pd.DataFrame(bad_rows).to_csv(args.failed, index=False)

    print("-" * 60)
    print(f"Resolved : {len(ok_rows):>4}  -> {args.output}")
    print(f"Failed   : {len(bad_rows):>4}  -> {args.failed}")
    if not HAVE_RDKIT:
        print("NOTE: RDKit not installed - desalting and validation used the "
              "weaker string heuristic. `pip install rdkit` is recommended.")


if __name__ == "__main__":
    main()
