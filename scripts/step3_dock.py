#!/usr/bin/env python3
"""
Step 3 - Dock every ligand against every protein with smina.

Expected layout (fpocket must have been run on each PDB first):

    proteins/
        pdb3sk3_ButyrateKinase.pdb
        pdb3sk3_ButyrateKinase_out/          <- created by `fpocket -f <pdb>`
            pockets/pocket1_atm.pdb
            pockets/pocket1_vert.pqr
    ligands/
        Benzoic_acid.pdbqt ...

Output: docking_results.csv -> Chemical,Protein,Score(kcal/mol),Status
        poses/<protein>__<chemical>.pdbqt

See docs/CHANGES.md for what differs from legacy/step3.py, and docs/CAVEATS.md
before you interpret any of these numbers.
"""

from __future__ import annotations

import argparse
import csv
import glob
import os
import re
import shutil
import subprocess
import sys

# smina/vina mode table:  "   1       -7.5      0.000      0.000"
MODE_LINE = re.compile(r"^\s*(\d+)\s+(-?\d+\.\d+)\s")

# Non-ligand heteroatoms that should never be part of a rigid receptor.
DROP_HET = {"HOH", "WAT", "DOD", "SO4", "PO4", "GOL", "EDO", "PEG", "MPD",
            "ACT", "CL", "NA", "K", "MG", "CA", "ZN", "IOD", "BME", "TRS",
            "DMS", "FMT", "NO3", "CIT", "TLA", "PGE", "1PE", "BOG"}


# --------------------------------------------------------------------------
# Receptor preparation
# --------------------------------------------------------------------------

def clean_pdb(pdb_in: str, pdb_out: str, keep_het: bool = False) -> int:
    """Keep model 1, drop waters/cryo-additives, keep altloc A. Returns #atoms."""
    kept = 0
    in_model = True
    with open(pdb_in, errors="ignore") as src, open(pdb_out, "w") as dst:
        for line in src:
            rec = line[:6]
            if rec == "MODEL ":
                in_model = line.split()[1] == "1" if len(line.split()) > 1 else True
                continue
            if rec == "ENDMDL":
                in_model = False
                continue
            if not in_model:
                continue
            if rec.startswith("ATOM") or rec.startswith("HETATM"):
                resname = line[17:20].strip().upper()
                altloc = line[16]
                element = line[76:78].strip().upper()
                if altloc not in (" ", "A"):
                    continue
                if element == "H":
                    continue  # re-added by Open Babel at the requested pH
                if rec.startswith("HETATM"):
                    if resname in DROP_HET or not keep_het:
                        continue
                dst.write(line)
                kept += 1
        dst.write("END\n")
    return kept


def prepare_receptor(pdb_file: str, workdir: str, ph: float, keep_het: bool,
                     overwrite: bool) -> str | None:
    """Clean the PDB and convert it to a rigid PDBQT receptor."""
    os.makedirs(workdir, exist_ok=True)
    stem = os.path.splitext(os.path.basename(pdb_file))[0]
    clean = os.path.join(workdir, f"{stem}_clean.pdb")
    pdbqt = os.path.join(workdir, f"{stem}.pdbqt")

    if os.path.exists(pdbqt) and os.path.getsize(pdbqt) > 0 and not overwrite:
        return pdbqt

    n = clean_pdb(pdb_file, clean, keep_het=keep_het)
    if n == 0:
        print(f"    [ERROR] no atoms left after cleaning {pdb_file}")
        return None

    proc = subprocess.run(
        ["obabel", clean, "-O", pdbqt, "-xr", "-p", str(ph)],
        capture_output=True, text=True,
    )
    if not os.path.exists(pdbqt) or os.path.getsize(pdbqt) == 0:
        print(f"    [ERROR] obabel failed on {clean}: "
              f"{(proc.stderr or '').strip()[:120]}")
        return None
    return pdbqt


# --------------------------------------------------------------------------
# Pocket geometry
# --------------------------------------------------------------------------

def _read_xyz(path: str, prefixes: tuple[str, ...]) -> list[tuple[float, float, float]]:
    """Coordinates from a PDB or PQR file.

    Tries fixed PDB columns first, then whitespace splitting - some tools write
    free-format PQR where the columns no longer line up.
    """
    coords = []
    with open(path, errors="ignore") as fh:
        for line in fh:
            if not line.startswith(prefixes):
                continue
            try:
                coords.append((float(line[30:38]), float(line[38:46]),
                               float(line[46:54])))
                continue
            except ValueError:
                pass
            fields = line.split()
            for start in (5, 6):  # with / without chain ID
                try:
                    coords.append((float(fields[start]), float(fields[start + 1]),
                                   float(fields[start + 2])))
                    break
                except (ValueError, IndexError):
                    continue
    return coords


def pocket_geometry(pdb_file: str, rank: int, padding: float,
                    box_min: float, box_max: float):
    """
    Center + box size for fpocket pocket `rank`.

    Prefers pocket<N>_vert.pqr (alpha-sphere centers - these describe the
    cavity itself) and falls back to pocket<N>_atm.pdb (the lining residues,
    which bias the center toward bulky side chains).
    """
    root = os.path.splitext(os.path.basename(pdb_file))[0]
    pocket_dir = os.path.join(os.path.dirname(pdb_file) or ".",
                              f"{root}_out", "pockets")

    vert = os.path.join(pocket_dir, f"pocket{rank}_vert.pqr")
    atm = os.path.join(pocket_dir, f"pocket{rank}_atm.pdb")

    if os.path.exists(vert):
        coords, source = _read_xyz(vert, ("ATOM", "HETATM")), "alpha_spheres"
    elif os.path.exists(atm):
        coords, source = _read_xyz(atm, ("ATOM", "HETATM")), "pocket_atoms"
    else:
        return None

    if len(coords) < 3:
        return None

    xs, ys, zs = zip(*coords)
    center = (sum(xs) / len(xs), sum(ys) / len(ys), sum(zs) / len(zs))
    size = tuple(
        min(box_max, max(box_min, (max(v) - min(v)) + 2 * padding))
        for v in (xs, ys, zs)
    )
    return {"center": center, "size": size, "source": source, "n": len(coords)}


# --------------------------------------------------------------------------
# Docking
# --------------------------------------------------------------------------

def parse_affinity(stdout: str) -> float | None:
    """Best-mode affinity from smina/vina stdout. None if not parseable."""
    for line in stdout.splitlines():
        m = MODE_LINE.match(line)
        if m and m.group(1) == "1":
            return float(m.group(2))
    return None


def dock(receptor: str, ligand: str, geom, out_pose: str, exhaustiveness: int,
         cpu: int, seed: int, timeout: int) -> tuple[float | None, str]:
    cx, cy, cz = geom["center"]
    sx, sy, sz = geom["size"]
    cmd = [
        "smina", "-r", receptor, "-l", ligand,
        "--center_x", f"{cx:.3f}", "--center_y", f"{cy:.3f}", "--center_z", f"{cz:.3f}",
        "--size_x", f"{sx:.1f}", "--size_y", f"{sy:.1f}", "--size_z", f"{sz:.1f}",
        "--exhaustiveness", str(exhaustiveness),
        "--num_modes", "1", "--cpu", str(cpu), "--seed", str(seed),
        "--out", out_pose,
    ]
    try:
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        return None, "timeout"

    score = parse_affinity(proc.stdout)
    if score is None:
        err = (proc.stderr or proc.stdout or "").strip().splitlines()
        return None, f"no_score: {err[-1][:70] if err else 'empty output'}"
    return score, "ok"


# --------------------------------------------------------------------------
# Driver
# --------------------------------------------------------------------------

def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("-p", "--proteins", default="proteins",
                   help="folder of .pdb files, each with its fpocket *_out folder "
                        "(default: %(default)s)")
    p.add_argument("-l", "--ligands", default="ligands",
                   help="folder of .pdbqt ligands from step 2 (default: %(default)s)")
    p.add_argument("-o", "--output", default="docking_results.csv",
                   help="results table (default: %(default)s)")
    p.add_argument("--poses", default="poses", help="folder for output poses")
    p.add_argument("--receptor-dir", default="receptors_prepared",
                   help="where cleaned/converted receptors are cached")
    p.add_argument("--pocket-rank", type=int, default=1,
                   help="which fpocket pocket to use; 1 = highest score "
                        "(default: %(default)s)")
    p.add_argument("--box-padding", type=float, default=4.0,
                   help="angstrom added on each side of the pocket extent "
                        "(default: %(default)s)")
    p.add_argument("--box-min", type=float, default=18.0, help="minimum box edge (A)")
    p.add_argument("--box-max", type=float, default=30.0, help="maximum box edge (A)")
    p.add_argument("--exhaustiveness", type=int, default=8,
                   help="smina sampling effort; 16 is a reasonable production value "
                        "(default: %(default)s)")
    p.add_argument("--cpu", type=int, default=4, help="CPUs per docking run")
    p.add_argument("--seed", type=int, default=42,
                   help="fixed RNG seed so runs are reproducible (default: %(default)s)")
    p.add_argument("--timeout", type=int, default=900, help="seconds per docking run")
    p.add_argument("--ph", type=float, default=6.5, help="pH for receptor protonation")
    p.add_argument("--keep-het", action="store_true",
                   help="keep non-water HETATM records (cofactors) in the receptor")
    p.add_argument("--overwrite", action="store_true",
                   help="ignore existing results and prepared receptors")
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)

    for tool in ("smina", "obabel"):
        if shutil.which(tool) is None:
            sys.exit(f"'{tool}' not found on PATH. See README.md > Installation.")

    proteins = sorted(glob.glob(os.path.join(args.proteins, "*.pdb")))
    proteins = [p for p in proteins if not p.endswith("_clean.pdb")]
    ligands = sorted(glob.glob(os.path.join(args.ligands, "*.pdbqt")))

    if not proteins:
        sys.exit(f"No .pdb files in {args.proteins}/")
    if not ligands:
        sys.exit(f"No .pdbqt files in {args.ligands}/ - run step 2 first.")

    print(f"{len(proteins)} proteins x {len(ligands)} ligands "
          f"= {len(proteins) * len(ligands)} docking runs")

    os.makedirs(args.poses, exist_ok=True)

    done = set()
    write_header = True
    if os.path.exists(args.output) and not args.overwrite:
        with open(args.output, newline="") as fh:
            for row in csv.DictReader(fh):
                done.add((row["Chemical"], row["Protein"]))
        write_header = False
        if done:
            print(f"Resuming: {len(done)} pairs already in {args.output}")

    mode = "w" if write_header or args.overwrite else "a"
    out = open(args.output, mode, newline="")
    writer = csv.writer(out)
    if mode == "w":
        writer.writerow(["Chemical", "Protein", "Score(kcal/mol)", "Status"])
        out.flush()

    for prot_pdb in proteins:
        prot_name = os.path.splitext(os.path.basename(prot_pdb))[0]
        print(f"\n=== {prot_name} ===")

        geom = pocket_geometry(prot_pdb, args.pocket_rank, args.box_padding,
                               args.box_min, args.box_max)
        if geom is None:
            print(f"  SKIP - no fpocket output found. Run:  fpocket -f {prot_pdb}")
            continue
        cx, cy, cz = geom["center"]
        sx, sy, sz = geom["size"]
        print(f"  pocket{args.pocket_rank} from {geom['source']} ({geom['n']} points)")
        print(f"  center ({cx:.2f}, {cy:.2f}, {cz:.2f})  box {sx:.0f}x{sy:.0f}x{sz:.0f} A")

        receptor = prepare_receptor(prot_pdb, args.receptor_dir, args.ph,
                                    args.keep_het, args.overwrite)
        if receptor is None:
            continue

        for lig in ligands:
            lig_name = os.path.splitext(os.path.basename(lig))[0]
            if (lig_name, prot_name) in done:
                continue
            pose = os.path.join(args.poses, f"{prot_name}__{lig_name}.pdbqt")
            score, status = dock(receptor, lig, geom, pose, args.exhaustiveness,
                                 args.cpu, args.seed, args.timeout)
            print(f"  {lig_name[:38]:<38} "
                  f"{score if score is not None else status}")
            writer.writerow([lig_name, prot_name,
                             "" if score is None else f"{score:.2f}", status])
            out.flush()

    out.close()
    print(f"\nDone -> {args.output}")


if __name__ == "__main__":
    main()
