# MUON.py

Converts a protein structure (PDB/CIF) + a rigid-body definition file +
a GROMACS-style force field directory into a LAMMPS `data_file`, per the
spec in Prof. Jignesh Prajapati's "Flow of the code" documents.

## Usage

```bash
python3 MUON.py -pdb new1.pdb -rigid rigid_definitions.dat \
    -ff /path/to/forcefield_dir -box 150 150 150 \
    -out data_file -log out.log [-ions]
```

`-ions` is optional (see "Ion inclusion" below); everything else is
required. This produces four intermediate files (`metadata.dat`,
`mass.dat`, `pairs.dat`, `atom.dat`) in the current directory plus the
final merged `data_file`, and writes a log describing every step to the
`-log` path.

## Validated against the reference output

Running the command above (without `-ions`) on `sample_data/new1.pdb`
reproduces `sample_data/reference_data_file` exactly, atom for atom, in
the `Atoms`/`Masses`/`Pair Coeffs` sections -- every atom's position,
charge, rigid body id, and atom type id matches. The only difference in
the full file is the atom type count (90 here vs. 88 in the reference),
which is expected: this `ffnonbonded.itp` has two additional entries
(`S4`, `KQ`) that were added after the reference file was generated
(confirmed with Prof. Jignesh -- both should be kept).

Reproduce the check:
```bash
cd sample_data
python3 ../MUON.py -pdb new1.pdb -rigid rigid_definitions.dat -ff ff \
    -box 150 150 150 -out generated_data_file -log out.log
diff generated_data_file reference_data_file
```

**Ion inclusion has not been validated against a reference file yet** --
Prof. Jignesh hasn't shared an ion-containing example yet. It has been
tested against small hand-built structures (see "Ion inclusion" below)
and matches the one worked example in his spec (a K+ ion's van der Waals
radius comes out to 1.519 A, exactly as in his document).

## Project layout

```
MUON.py                     launch script -- wires everything together
muon/
  cli.py                    command-line argument parsing
  structure_reader.py       PDB/CIF -> atoms matrix, chain relabeling,
                             ion padding (-ions)
  forcefield_reader.py      ffnonbonded.itp -> nonbonded matrix,
                             .rtp/.r2b lookup for atom type + charge,
                             histidine protonation state resolution,
                             ions.itp -> ion list + van der Waals radii
  rigid_parser.py           rigid_definitions.dat -> rigid body ids
                             (supports both a fixed id and an id range)
  writers.py                metadata.dat / mass.dat / pairs.dat /
                             atom.dat writers + final assembly
sample_data/                files Prof. Jignesh provided, used both as a
                             working example and as a regression test
```

## Design notes (clarified with Prof. Jignesh)

A few details weren't fully spelled out in the original "Flow of the
code" note. These were resolved by comparing output against the
reference `data_file` and then confirmed directly with him:

1. **Terminal residues use different rtp blocks.** The first and last
   residue of each chain have different atom names (e.g. `H1/H2/H3`
   instead of `H` at the N-terminus, `OC1/OC2` instead of `O` at the
   C-terminus), looked up via the N-ter/C-ter columns of the `.r2b`
   files (e.g. `NALA`, `CSER`) rather than the plain residue name.

2. **Histidine protonation state is resolved per residue.** A plain PDB
   labels every histidine `HIS` regardless of its actual protonation
   state. For each histidine, the code tries matching it as `HID` first;
   if any of its atoms don't fit that block, it tries `HIE`, then `HIP`
   (terminal histidines use the N-ter/C-ter variant of each, e.g.
   `NHID`/`CHID`).

3. **Atom coordinates are recentered using the structure's own
   centroid** (mean x, y, z across all atoms), which is the expected
   convention for a structure prepared in GROMACS.

4. **Zero-mass atom types are written as `0.0001`, not `0.0000`.**
   A few force-field entries are massless virtual/dummy sites (`MW`,
   `MCH3`, `MNH3`, `EP`); LAMMPS requires a strictly positive mass.

## Ion inclusion (`-ions`, from Flow of the code_v2)

LAMMPS requires every rigid body to have at least 3 atoms. A lone ion
(e.g. a single K+) is only 1 atom, so when `-ions` is passed, every atom
belonging to a known ion residue gets two extra "padding" atoms inserted
right after it in the atoms list:

- Both copy the ion's chain, residue id, y, and z.
- One is shifted `-radius` along x, the other `+radius` along x, where
  `radius` is that specific ion's van der Waals radius.
- Both get atom name `MW` **and** residue name `MW`, atom type `MW` (an
  existing zero-mass, zero-charge dummy type already in
  `ffnonbonded.itp`), and charge `0` -- this bypasses the normal rtp
  lookup entirely for these two atoms.

**Where the van der Waals radius comes from**, per ion:
1. The list of ion residue names (`IB+, CA, CL, NA, MG, K, RB, CS, LI,
   ZN`) comes from `ions.itp`.
2. Each ion's force field atom type comes from its own block in
   `aminoacids.rtp` (ions are defined there, not in a separate file).
   Note: the atom's name inside that block doesn't always match the
   residue name -- e.g. residue `IB+`'s one atom is named `IB`.
3. That atom type's `sigma` value comes from the nonbonded matrix
   (built from `ffnonbonded.itp`, same as everywhere else in the
   pipeline).
4. `radius = (sigma_nm * 10) / 2` (nm -> Angstrom, diameter -> radius).

This step runs *before* chain relabeling and atom renumbering, so the
new dummy atoms get included in both like any other atom, and *before*
the rtp/charge lookup loop, which checks for residue name `MW` first
and assigns its type/charge directly rather than searching the rtp
files (there's no rtp entry for a residue named MW).

**Clarified with Prof. Jignesh:** his original spec said to check
whether "the residue name is MW", but his own worked example only
changed the *atom* name to `MW`, not the residue name. He confirmed
(2026-10-12) that the residue name should be changed to `MW` as well,
to match his original instruction -- both are now set.

## Rigid body id ranges (`rigid_definitions.dat`, general, not ion-specific)

The third column of `rigid_definitions.dat` can now be either a single
number (the whole residue range shares one rigid body, as before) or a
range, e.g.:
```
C1   110-114   2-6   0
```
assigns residue 110 -> rigid id 2, 111 -> 3, ..., 114 -> 6 (one id per
residue, in order). This is essential for ions: a run of many individual
ions, each needing its own separate rigid body, would otherwise require
one `rigid_definitions.dat` row per ion. `rigid_parser.py` raises an
error if the residue count in the range doesn't match the id count, to
catch a mismatched row rather than silently assigning the wrong ids.

## Open items

- **CIF testing**: `structure_reader.py` reads both PDB and CIF files
  into the same atoms-matrix format, but this has only been exercised
  against PDB input so far.
- **Ion validation**: no reference `data_file` with ions exists yet to
  diff against (see "Validated against the reference output" above).

## Extending

- More CLI flags: add them in `muon/cli.py`.
- Flexible rigid types (1, 2): `rigid_parser.apply_rigid_ids` already
  passes `rigid_type` through unchanged for when that logic is needed.
- DNA/RNA structures: `forcefield_reader.ResidueDatabase` already tries
  `dna.rtp` / `rna.rtp` as fallbacks, but hasn't been tested against a
  nucleic acid structure yet.
