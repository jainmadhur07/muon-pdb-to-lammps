# MUON.py — Implementation Notes
---

## 1. Pipeline overview

`MUON.py` runs five stages, matching the spec's step numbering:

| Stage | Spec step(s) | Module |
|---|---|---|
| Parse CLI args | Section 2 | `muon/cli.py` |
| Read structure (PDB/CIF), renumber atoms, relabel chains | Steps 1, 4a, 4c | `muon/structure_reader.py` |
| Parse `ffnonbonded.itp` → nonbonded matrix | Step 2 | `muon/forcefield_reader.py` |
| Resolve atom type + charge per atom (`.rtp`/`.r2b` lookup) | Step 4b | `muon/forcefield_reader.py` |
| Parse `rigid_definitions.dat`, assign rigid body IDs | Step 4d | `muon/rigid_parser.py` |
| Write `metadata.dat`, `mass.dat`, `pairs.dat`, `atom.dat`; concatenate into `data_file` | Steps 1–3, 4e, 5 | `muon/writers.py` |

Internally, the "atoms matrix" is a list of per-atom dicts rather than a NumPy array; each stage adds columns (keys) to it rather than replacing it, which keeps the data flow linear and easy to log at every step.

---

## 2. Behaviors not specified in the original note

Four things had to be reverse-engineered by diffing output against the reference `data_file`, since the spec doesn't mention them. **All four have since been confirmed correct by Prof. Jignesh** (see commit history / prior email thread); summarized here for anyone else reading the repo.

### 2.1 Terminal residues require different `.rtp` blocks

The first and last residue of every chain use different atom names than mid-chain residues (e.g. N-terminal alanine has `H1/H2/H3` instead of a single `H`), corresponding to different `.rtp` blocks (`NALA`, `CALA`, etc.) rather than the plain residue name. The mapping from plain residue name to its N-ter/C-ter block name is read from the accompanying `.r2b` file.

Implementation: `structure_reader.relabel_chains_and_atoms()` tags each atom `is_nterm`/`is_cterm` based on whether it belongs to the lowest/highest-numbered residue in its chain. `forcefield_reader.ResidueDatabase._resolve_block_name()` uses those flags to pick the correct block name before doing the atom-name lookup.

**Confirmed correct.**

### 2.2 Histidine protonation state is resolved per-residue, not globally

The input structure contains histidines in more than one protonation state (some with `HD1` present, others with `HE2` present). Each histidine residue's actual hydrogen placement — not a fixed assumption — determines whether it resolves to `HID`, `HIE`, or `HIP`.

Implementation: `forcefield_reader.ResidueDatabase.resolve_histidine_block()` tries `HID` first, then `HIE`, then `HIP` (via their GROMACS-convention names `HISD`/`HISE`/`HISH`, mapped through `aminoacids.r2b`), accepting the first candidate block whose atom list is a superset of the residue's actual atom names. Terminal histidines are handled the same way, using the N-ter/C-ter variants of each candidate.

**Confirmed correct** (try order: HID → HIE → HIP, as specified).

### 2.3 Atom coordinates are recentered on the structure's own centroid

Output coordinates are the input coordinates minus the mean (x, y, z) across all atoms — not a shift by half the box length, which only approximately matched. Confirmed as expected behavior for structures prepared in GROMACS.

Implementation: computed once in `MUON.py` from the full atom list, passed into `writers.write_atoms()` and subtracted per-atom.

**Confirmed correct.**

### 2.4 Zero-mass force field entries are written as `0.0001`

A small number of `ffnonbonded.itp` entries (virtual/dummy sites, e.g. `MW`, `MCH3`, `MNH3`, `EP`) have mass `0.0000`, which LAMMPS rejects. These are written to `mass.dat` as `0.0001` to match the reference file and satisfy LAMMPS's positive-mass requirement.

Implementation: `writers.write_mass()`.

**Confirmed correct.**

---

## 3. Other clarifications received

- **Atom type count:** the reference `data_file` (88 types) was generated with an older version of the force field folder; the current `ffnonbonded.itp` (90 types, including `S4`/`KQ`) should be used in full, not filtered down to match the older reference.
- **Reference file provenance:** the reference `data_file` was generated from `new1.pdb`, not `new2.pdb`. `new2.pdb` is intended for a later stage.

---

## 4. Validation

Running the pipeline on `sample_data/new1.pdb` reproduces `sample_data/reference_data_file` **exactly**, byte-for-byte, in the `Masses`, `Pair Coeffs`, and `Atoms` sections:

```bash
cd sample_data
python3 ../MUON.py -pdb new1.pdb -rigid rigid_definitions.dat -ff ff \
    -box 150 150 150 -out generated_data_file -log out.log
diff generated_data_file reference_data_file
```

No diff output (aside from the atom-type-count line, expected per §3).

---

## 5. Known gaps / not yet exercised

- **DNA/RNA structures.** `ResidueDatabase` falls through to `dna.rtp`/`rna.rtp` when a residue isn't found in `aminoacids.rtp`, per the spec's fallback instruction, but this path hasn't been tested against an actual nucleic acid structure yet.
- **Rigid types 1/2.** `rigid_parser.apply_rigid_ids()` passes `rigid_type` through unchanged; only type `0` appears in the current `rigid_definitions.dat`, and no logic exists yet for the flexible-model types (by design — not needed until that feature exists).
- **CIF input.** `structure_reader.read_cif()` handles the standard `_atom_site` loop and prefers `auth_asym_id`/`auth_seq_id` over `label_*` fields so that CIF-derived residue numbering stays consistent with `rigid_definitions.dat`'s conventions. Not yet run against a real CIF file — only PDB inputs have been tested so far.

---

## 6. Project layout

```
MUON.py                launch script
muon/
  cli.py                CLI argument parsing
  structure_reader.py    PDB/CIF -> atoms list, chain relabeling, terminal-residue tagging
  forcefield_reader.py   ffnonbonded.itp -> nonbonded matrix; .rtp/.r2b residue lookup
  rigid_parser.py        rigid_definitions.dat -> rigid body IDs
  writers.py              metadata.dat / mass.dat / pairs.dat / atom.dat + final assembly
sample_data/            reference input/output files used as a regression test
```
