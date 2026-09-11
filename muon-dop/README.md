# MUON.py

Converts a protein structure (PDB/CIF) + a rigid-body definition file +
a GROMACS-style force field directory into a LAMMPS `data_file`, per the
spec in Prof. Jignesh Prajapati's "Flow of the code" document.

## Usage

```bash
python3 MUON.py -pdb new1.pdb -rigid rigid_definitions.dat \
    -ff /path/to/forcefield_dir -box 150 150 150 \
    -out data_file -log out.log
```

This produces four intermediate files (`metadata.dat`, `mass.dat`,
`pairs.dat`, `atom.dat`) in the current directory plus the final merged
`data_file`, and writes a log describing every step to the `-log` path.

## Validated against the reference output

Reproduce the check:
```bash
cd sample_data
python3 ../MUON.py -pdb new1.pdb -rigid rigid_definitions.dat -ff ff \
    -box 150 150 150 -out generated_data_file -log out.log
diff generated_data_file reference_data_file
```

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
