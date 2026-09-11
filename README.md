# MUON.py

Converts a protein structure (PDB/CIF) + a rigid-body definition file +
a GROMACS-style force field directory into a LAMMPS `data_file`,

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

## 4. Validation

Reproduce the check:
```bash
cd sample_data
python3 ../MUON.py -pdb new1.pdb -rigid rigid_definitions.dat -ff ff \
    -box 150 150 150 -out generated_data_file -log out.log
diff generated_data_file reference_data_file
```
## Usage

```bash
python3 MUON.py -pdb new1.pdb -rigid rigid_definitions.dat \
    -ff /path/to/forcefield_dir -box 150 150 150 \
    -out data_file -log out.log
```

This produces four intermediate files (`metadata.dat`, `mass.dat`,
`pairs.dat`, `atom.dat`) in the current directory plus the final merged
`data_file`, and writes a log describing every step to the `-log` path.


## Project layout

```
MUON.py                     launch script -- wires everything together
muon/
  cli.py                    command-line argument parsing
  structure_reader.py       PDB/CIF -> atoms matrix, chain relabeling
  forcefield_reader.py      ffnonbonded.itp -> nonbonded matrix,
                             .rtp/.r2b lookup for atom type + charge,
                             histidine protonation state resolution
  rigid_parser.py           rigid_definitions.dat -> rigid body ids
  writers.py                metadata.dat / mass.dat / pairs.dat /
                             atom.dat writers + final assembly
sample_data/                files Prof. Jignesh provided, used both as a
                             working example and as a regression test
```
