#!/usr/bin/env python3
"""MUON.py -- launch script.

Usage:
    python3 MUON.py -pdb new1.pdb -rigid rigid_definitions.dat \
        -ff /path/to/forcefield_dir -box 150 150 150 \
        -out data_file -log out.log
"""
import sys
import time

from muon import cli
from muon.structure_reader import (
    read_structure, relabel_chains_and_atoms, group_residues,
)
from muon.forcefield_reader import read_nonbonded_matrix, ResidueDatabase
from muon.rigid_parser import read_rigid_definitions, apply_rigid_ids
from muon.writers import write_metadata, write_mass, write_pairs, write_atoms, assemble_data_file


def log(fh, message):
    timestamp = time.strftime("%Y-%m-%d %H:%M:%S")
    line = f"[{timestamp}] {message}"
    print(line)
    fh.write(line + "\n")
    fh.flush()


def main(argv=None):
    args = cli.parse_args(argv)

    with open(args.log_file, "w") as logfh:
        log(logfh, "Starting MUON.py")
        log(logfh, f"structure={args.structure_file} rigid={args.rigid_file} "
                    f"ff_dir={args.ff_dir} box={args.box} out={args.out_file}")

        atoms = read_structure(args.structure_file)
        log(logfh, f"Read {len(atoms)} atoms from {args.structure_file}")

        atoms = relabel_chains_and_atoms(atoms)
        chains_seen = sorted(set(a["new_chain"] for a in atoms))
        log(logfh, f"Relabeled chains: {chains_seen}")

        nonbonded = read_nonbonded_matrix(args.ff_dir)
        log(logfh, f"Read {len(nonbonded)} atom types from ffnonbonded.itp")
        fftype_to_id = {row["name"]: row["id"] for row in nonbonded}

        resdb = ResidueDatabase(args.ff_dir)

        # Histidine: try HID first, then HIE, then HIP (per residue).
        residues = group_residues(atoms)
        histidine_blocks = {}  # (new_chain, resid) -> (block_name, family)
        for key, residue_atoms in residues.items():
            if residue_atoms[0]["resname"] != "HIS":
                continue
            atom_names = {a["name"] for a in residue_atoms}
            is_nterm = residue_atoms[0]["is_nterm"]
            is_cterm = residue_atoms[0]["is_cterm"]
            try:
                block_name, family = resdb.resolve_histidine_block(
                    atom_names, is_nterm, is_cterm
                )
            except LookupError as e:
                log(logfh, f"ERROR: {e}")
                sys.exit(1)
            histidine_blocks[key] = (block_name, family)
            log(logfh, f"Histidine at {key} resolved to rtp block '{block_name}'")

        # Look up atom type + charge for every atom.
        not_found = []
        for atom in atoms:
            key = (atom["new_chain"], atom["resid"])
            try:
                if key in histidine_blocks:
                    block_name, family = histidine_blocks[key]
                    fftype, charge = resdb.rtp[family][block_name][atom["name"]]
                else:
                    fftype, charge, family = resdb.lookup(
                        atom["resname"], atom["name"], atom["is_nterm"], atom["is_cterm"]
                    )
            except (LookupError, KeyError):
                not_found.append(
                    f"Could not find atom '{atom['name']}' of residue "
                    f"'{atom['resname']}' at {key}"
                )
                continue
            atom["fftype"] = fftype
            atom["charge"] = charge
            if fftype not in fftype_to_id:
                not_found.append(
                    f"Atom type '{fftype}' (residue {atom['resname']}, atom "
                    f"{atom['name']}) not found in ffnonbonded.itp"
                )
                continue
            atom["atom_type_id"] = fftype_to_id[fftype]

        if not_found:
            for msg in not_found:
                log(logfh, f"ERROR: {msg}")
            log(logfh, f"Aborting: {len(not_found)} atom(s) could not be resolved.")
            sys.exit(1)

        log(logfh, "All atoms matched to a force field type and charge.")

        rigid_defs = read_rigid_definitions(args.rigid_file)
        apply_rigid_ids(atoms, rigid_defs)
        log(logfh, f"Applied {len(rigid_defs)} rigid body definition(s).")

        n_atom_types = len(nonbonded)
        write_metadata(
            "metadata.dat",
            header_comment="# LAMMPS data file for rigid bodies",
            n_atoms=len(atoms),
            n_atom_types=n_atom_types,
            box=args.box,
        )
        log(logfh, "Wrote metadata.dat")

        write_mass("mass.dat", nonbonded)
        write_pairs("pairs.dat", nonbonded)
        log(logfh, "Wrote mass.dat and pairs.dat")

        # Recenter the structure at the origin using its own centroid.
        n = len(atoms)
        centroid = (
            sum(a["x"] for a in atoms) / n,
            sum(a["y"] for a in atoms) / n,
            sum(a["z"] for a in atoms) / n,
        )
        log(logfh, f"Structure centroid (subtracted from all coordinates): {centroid}")
        write_atoms("atom.dat", atoms, box_center_shift=centroid)
        log(logfh, "Wrote atom.dat")

        assemble_data_file(args.out_file, ["metadata.dat", "mass.dat", "pairs.dat", "atom.dat"])
        log(logfh, f"Wrote final data file: {args.out_file}")
        log(logfh, "Done.")


if __name__ == "__main__":
    main()
