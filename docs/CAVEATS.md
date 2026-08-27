# What this pipeline can and cannot tell you

Read this before a number from `druggable.csv` goes into a figure, a grant, or a
manuscript. Nothing here is a bug — these are the limits of the method itself.

**The one-line summary:** this pipeline generates ranked hypotheses. It does not
provide evidence that a chemical binds a protein.

---

## 1. Docking scores scale with molecular size

Vina-family scoring functions (smina's default is Vina) sum favourable atomic
contacts, so the score grows with heavy-atom count almost automatically. Across
a chemically diverse library, a large flexible molecule will out-score a small
rigid one even when the small one is the better, more specific binder.

Consequences for this design:

* Any absolute cutoff (the `≤ −8.0 kcal/mol` lead filter) is a **size filter
  wearing an affinity filter's clothes**. Benzoic acid (11 heavy atoms) cannot
  reach −8 kcal/mol in practice; a maltodextrin fragment reaches it against
  almost anything.
* The per-chemical z-score is the right instinct: normalising each chemical
  against its own mean cancels the size term, because size is constant along a
  row. Keep that as the primary ranking.
* Report **ligand efficiency** (score ÷ heavy atoms) alongside the raw score.
  `step4_analyze.py` computes it when RDKit and `ligands/ligand_index.csv` are
  available. Typical guidance: LE better (more negative) than about
  −0.30 kcal/mol per heavy atom is interesting; worse than −0.20 usually is not.

## 2. Polymers cannot be docked

Gum arabic, carrageenan, carob bean gum, tragacanth and karaya gum are
heterogeneous polysaccharides of 10⁴–10⁶ Da. Their PubChem records are a
representative repeat unit or one component sugar. Docking that fragment tells
you nothing about the polymer: the entropy, the avidity, the conformational
ensemble and the actual chemistry present in the food matrix are all different.

Either remove them from the chemical list up front, or state explicitly in the
methods that what was screened was a monomer surrogate, and do not extrapolate
the result to the additive.

## 3. Rigid receptor, single ligand protonation state

* The receptor is rigid. No induced fit, no side-chain rearrangement. Pockets
  that only open in the holo form will be missed entirely.
* All structural waters are removed. Bridging waters are catalytically or
  energetically essential in many active sites; removing them can both create
  false positives (an artificially open pocket) and false negatives.
* Open Babel's `-p 6.5` assigns one protonation state by simple rules. It does
  not enumerate tautomers or alternative ionisation states. For carboxylic
  acids this matters: benzoic acid (pKa ≈ 4.2) is essentially fully
  deprotonated at pH 6.5, and the charged and neutral forms dock very
  differently. If ionisation is central to your argument, prepare ligands with
  a dedicated tool (Dimorphite-DL, Epik, or `obabel` at several pH values) and
  dock each state.
* Open Babel-generated PDBQT files are not identical to AutoDockTools /
  Meeko output (charge model, rotatable-bond perception, rigid rings). smina
  accepts them, but if you compare against published AutoDock numbers, prepare
  ligands the same way they did.

## 4. fpocket pocket 1 is not necessarily the active site

`pocket1` is the pocket fpocket ranks highest by its own druggability score.
That is often, but by no means always, the functional site. For every protein
in the panel you should confirm the chosen pocket either overlaps the known
catalytic residues, or coincides with a co-crystallised ligand, or at minimum
report which pocket rank was used. `step3_dock.py` accepts `--pocket-rank`;
consider docking the top 2–3 pockets and reporting all of them.

## 5. The z-score has no null distribution

An extreme z-score means "unusual within this particular matrix", not
"significant". With C chemicals × P proteins you are looking at C·P values; some
will reach z ≤ −1.5 by chance alone. To make a specificity claim you need a
null:

* **Property-matched decoys** (DUD-E style) docked into the same pockets — the
  standard way to show the score is doing something beyond matching physical
  properties.
* **Scrambled control:** permute the chemical–protein assignment and re-run the
  filter; how many "hits" survive?
* **Positive control:** re-dock a known co-crystal ligand into its own structure
  and confirm the pose is recovered (heavy-atom RMSD < 2 Å). A pipeline that
  cannot re-dock a known ligand cannot be trusted on unknown ones.

Also note that the two filters are not independent: `Chem_Std > 1.5` and
`Z_Chem ≤ −1.5` are both functions of the same row.

## 6. Docking is stochastic; report the spread

smina's search is a Monte-Carlo procedure. `step3_dock.py` fixes `--seed 42` so
a run is reproducible, but reproducible is not the same as stable. Re-run the
shortlist with 3–5 different seeds and at higher `--exhaustiveness` (16–32). If
the run-to-run standard deviation is comparable to the differences you are
ranking on, the ranking is noise.

## 7. Score accuracy

Published benchmarks put Vina-type scoring at roughly 2–3 kcal/mol RMS error
against experimental affinities, with rank correlations that are modest on
diverse targets. Practical rule: **differences smaller than ~1.5 kcal/mol are
not interpretable**, and a numeric score should never be converted to a Kd in a
figure.

---

## A defensible minimal protocol

If these results are going into a manuscript, the reviewer-proof version of this
work is roughly:

1. Drop polymers; screen only well-defined small molecules.
2. Verify each pocket against known site information; report the pocket rank.
3. Re-dock a native ligand per structure as a positive control (RMSD < 2 Å).
4. Dock a property-matched decoy set through the identical pipeline.
5. Rank on per-chemical z-score **and** ligand efficiency, not raw score.
6. Repeat the shortlist with ≥3 seeds at exhaustiveness ≥16; report mean ± SD.
7. Take the top few hits to an orthogonal assay — ITC, MST, SPR, thermal shift,
   or an enzyme activity assay. Nothing above substitutes for this.
