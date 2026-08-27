# legacy/

The original four scripts, unmodified, exactly as they were first written.

They are kept here for provenance: if you have results produced before the
rewrite, this is the code that produced them, and this is what you cite or
re-run to reproduce them.

**Do not use these for new work.** They contain several bugs that either crash
the run or silently produce wrong numbers - most importantly the score parser in
`step3.py`, which can record `100.0` for a docking run that actually succeeded.
Every issue is listed in [`../docs/CHANGES.md`](../docs/CHANGES.md), with the
fix that replaced it in `../scripts/`.

| legacy file | replaced by |
|---|---|
| `step1.py` | `../scripts/step1_fetch_smiles.py` |
| `step2.py` | `../scripts/step2_build_ligands.py` |
| `step3.py` | `../scripts/step3_dock.py` |
| `step4.py` | `../scripts/step4_analyze.py` |
| `step1_input_file_example.csv` | `../data/example_chemicals.csv` |
