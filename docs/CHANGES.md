# What changed vs. the original scripts

The originals are preserved untouched in [`../legacy/`](../legacy/). This file
lists every substantive change, why it was made, and whether it can change your
numbers.

Legend: **[crash]** the original could fail outright · **[silent]** the original
could produce wrong numbers without telling you · **[result]** the fix changes
results, so old and new runs are not directly comparable · **[qol]** convenience
only.

---

## step1.py → `scripts/step1_fetch_smiles.py`

| # | Issue in the original | Fix |
|---|---|---|
| 1 | **[crash]** `compounds[0].isomeric_smiles`. PubChem retired the `IsomericSMILES` PUG-REST property in favour of `SMILES` / `ConnectivitySMILES`. PubChemPy 1.0.5 keeps `isomeric_smiles` only as a deprecated alias of `.smiles`; on 1.0.4 and earlier it still requests the retired property. So the line either warns or returns nothing depending on which version is installed. | Try `.smiles` → `.isomeric_smiles` → `.connectivity_smiles` → `.canonical_smiles`, then fall back to a direct `pcp.get_properties()` call. Works on both versions. |
| 2 | **[silent] [result]** Desalting kept fragments with `len(part) > 3`. String length is not molecular size: `[Na+]` is 5 characters and survives; ethanol `CCO` is 3 and is thrown away. | RDKit largest-fragment by heavy-atom count, plus an explicit counterion list. Heuristic heavy-atom estimate when RDKit is absent. |
| 3 | **[result]** Polymer filter was `len(smiles) > 350` characters. SMILES length is inflated by branching, ring closures and stereo tags, so the cutoff does not correspond to any physical size. | Molecular-weight cutoff (`--max-mw`, default 900 Da) with a heavy-atom fallback. **This changes which compounds are skipped.** |
| 4 | **[silent]** A bare `except Exception: return None` made a PubChem outage indistinguishable from "compound not in PubChem". | Network errors are separated from misses and retried with exponential backoff (`--retries`). |
| 5 | **[silent]** Failures were printed to the terminal and then lost. A run that resolved 60 of 200 compounds left no record of the other 140. | Everything that fails is written to `chemicals_failed.csv` with a machine-readable reason. |
| 6 | **[qol]** Output kept only `Name,SMILES`, discarding CAS. Nothing downstream can be traced back to a registry number. | Output carries `Name,CAS,CID,SMILES,MolecularWeight,Source`. |
| 7 | **[qol]** `time.sleep(0.2)` was applied per *compound*, but up to two PubChem requests were issued per compound. | Pause applied per request (`--pause`); PubChem asks for ≤5 requests/s. |
| 8 | **[qol]** Not restartable — an interrupted run started from zero. | Resumes from an existing output file unless `--overwrite`. |

**Not fixed, by design:** the script still takes `compounds[0]`, i.e. PubChem's
first hit. For trivial names that is a real ambiguity. Querying CAS first
mitigates it, and the `CID` column now lets you spot-check. Review the CIDs
before you trust a screen.

---

## step2.py → `scripts/step2_build_ligands.py`

| # | Issue in the original | Fix |
|---|---|---|
| 1 | **[crash]** The SMILES string was interpolated into a shell command string with `shell=True`. Any SMILES containing `$`, a backtick, a quote or a backslash breaks the command — and in principle executes whatever follows. | `subprocess.run()` with an argument list. No shell involved. |
| 2 | **[crash]** Relied on the GNU `timeout` binary, which is not present on macOS by default. | Python's own `subprocess` timeout. Portable. |
| 3 | **[silent]** Exit code 0 was treated as success. Open Babel can exit 0 after writing an empty or single-atom `.pdbqt` when 3D embedding fails; that file then goes into docking and produces meaningless scores. | The output file is validated (≥3 `ATOM`/`HETATM` records) before a molecule counts as built. |
| 4 | **[silent]** Two chemicals whose names sanitize to the same stem overwrote each other's file, so one compound was silently replaced by the other. | Collisions get `_2`, `_3` … suffixes. |
| 5 | **[silent]** No mapping from the sanitized file name back to the compound. Steps 3 and 4 only ever see the file stem, so a hit could not be traced to a CAS or SMILES. | `ligands/ligand_index.csv` records `SafeName ↔ Name, CAS, CID, SMILES, MW, Status`. |
| 6 | **[qol]** `obabel` was never checked for; a missing binary produced 127 for every molecule. | Checked up front with a clear install message. |
| 7 | **[qol]** Resume logic was commented out. | `--overwrite` flag; reuse is the default. |
| 8 | **[qol]** `--gen3d` quality and force-field relaxation were not reachable. | `--gen3d`, `--minimize FF`, `--steps`. |

---

## step3.py → `scripts/step3_dock.py`

This is where the most consequential bugs were.

| # | Issue in the original | Fix |
|---|---|---|
| 1 | **[silent]** Score parsing: `if line.strip().startswith("1")`. This matches any output line beginning with the digit 1, and because the loop never breaks, the **last** match wins rather than the mode-1 row. When nothing parsed, the score defaulted to `100.0` — a sentinel written into the results file as if it were a measurement. Step 4 then removed it with its `< 0` filter, so failed runs disappeared without ever being counted. | Anchored regex on the mode table, first match only, plus an explicit `Status` column. A failed run writes an empty score and its reason, and stays visible. |
| 2 | **[result]** Receptor preparation was `obabel prot.pdb -O prot.pdbqt -xr`, which converts *everything* in the file: crystallographic waters, glycerol, sulfate, cryoprotectants, and any co-crystallised ligand. That ligand very often sits in pocket 1 — the pocket you are docking into — and blocks it. | The PDB is cleaned first: model 1 only, alternate conformers other than A dropped, existing hydrogens stripped before re-protonating at the requested pH, and HETATM removed (`--keep-het` to retain cofactors). **This changes results, usually for the better.** |
| 3 | **[result]** Pocket center = mean of the atoms in `pocket1_atm.pdb`, i.e. the residues lining the cavity. That mean is pulled toward whichever wall has more atoms and toward bulky side chains, not toward the cavity. | Prefers `pocket1_vert.pqr` (fpocket's alpha-sphere centers, which describe the cavity volume itself) and falls back to `_atm.pdb`. The script prints which one it used. |
| 4 | **[result]** Fixed 20 × 20 × 20 Å box for every pocket. Too small for a large site (it truncates the accessible space and biases the pose) and too large for a small one (sampling is wasted and the ligand can drift into a neighbouring surface groove). | Box derived from the pocket extent plus `--box-padding`, clamped to `--box-min`/`--box-max` (18–30 Å). |
| 5 | **[result]** `--exhaustiveness` was never set, so smina's default of 8 was used, and no `--seed` was given, so no run was reproducible. | Both exposed; seed fixed at 42 by default. Use `--exhaustiveness 16` or higher for anything you intend to publish. |
| 6 | **[qol]** `prot_pdb.replace(".pdb", ".pdbqt")` wrote prepared receptors into the input folder. | Cached in `receptors_prepared/`, path built with `os.path.splitext`. |
| 7 | **[qol]** No resume. A crash 9 000 pairs into a 10 000-pair run meant starting over. | Existing `(Chemical, Protein)` pairs are read from the results file and skipped. |
| 8 | **[qol]** Pose files named `{protein}_{ligand}`, ambiguous when either name contains `_`. | `{protein}__{ligand}` in `poses/`. |
| 9 | **[qol]** No check for `smina` / `obabel`; no timeout on a hung docking run. | Both added. |

---

## step4.py → `scripts/step4_analyze.py`

| # | Issue in the original | Fix |
|---|---|---|
| 1 | **[silent]** `std()` uses `ddof=1`, so a chemical docked against a single protein gets `NaN`. `NaN > 1.5` and `NaN <= -1.5` are both `False`, so that chemical was dropped with no message. More generally, a z-score computed over 2–3 proteins is not interpretable. | `--min-proteins` (default 5), and a printed audit showing how many rows and chemicals were removed at each stage. |
| 2 | **[silent]** The score column was never coerced to numeric, so a mixed column (numbers plus blanks or `100.0` sentinels) compared in ways that are hard to predict. | `pd.to_numeric(errors="coerce")`, then explicit dropping. |
| 3 | **[silent]** Re-running step 3 appends, so duplicate `(chemical, protein)` pairs were counted more than once and distorted the mean and standard deviation. | Deduplicated, keeping the best score per pair. |
| 4 | **[result]** Only the per-chemical (row) z-score was computed. A large hydrophobic pocket scores *every* ligand well, so it looks like a specific hit for all of them. | Per-protein (column) z-score added as `Z_Prot`; optionally filter on it with `--protein-z-threshold`. |
| 5 | **[result]** The `≤ −8.0 kcal/mol` pre-filter is size-biased. Vina-type scores scale roughly with heavy-atom count, so this filter mostly selects large molecules and structurally excludes small ones — benzoic acid (11 heavy atoms) will essentially never reach −8, no matter how well it fits. | The filter is now optional (`--no-lead-filter`, `--score-lead-min`), and ligand efficiency (score ÷ heavy atoms) is reported when RDKit and `ligand_index.csv` are available. See [CAVEATS.md](CAVEATS.md). |
| 6 | **[qol]** No wide-form output. | `docking_matrix.csv` (chemical × protein) for heatmaps and supplements. |
| 7 | **[qol]** Thresholds hard-coded at the top of the file. | `argparse`, with the same values as defaults. |

---

## Things deliberately **not** changed

* The overall four-step design, the file formats, and the default threshold
  values (`z ≤ −1.5`, `std > 1.5`, `≤ −8.0 kcal/mol`, pH 6.5) are unchanged, so
  a rerun is comparable to your previous work except where a **[result]** row
  above says otherwise.
* No decoy set, promiscuity control, rescoring or MD refinement has been added.
  Those are study-design decisions, not bugs. [CAVEATS.md](CAVEATS.md) explains
  why they matter before these numbers go into a manuscript.
