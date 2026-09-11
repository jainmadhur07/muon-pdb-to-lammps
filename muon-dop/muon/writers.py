"""Write the intermediate .dat files and merge them into the final
LAMMPS data_file, matching the reference file's format exactly."""


def write_metadata(path, header_comment, n_atoms, n_atom_types, box):
    """metadata.dat: header comment, atom count, atom type count, box bounds."""
    bx, by, bz = box
    with open(path, "w") as fh:
        if header_comment:
            fh.write(f"{header_comment}\n\n")
        fh.write(f"{n_atoms} atoms\n\n")
        fh.write(f"{n_atom_types}  atom types\n\n")
        fh.write(f"0 {bx:g} xlo xhi\n")
        fh.write(f"0 {by:g} ylo yhi\n")
        fh.write(f"0 {bz:g} zlo zhi\n\n")


def write_mass(path, nonbonded_matrix):
    """mass.dat: 'Masses' header, blank line, then 'id mass' per line.
    Zero-mass dummy atom types (e.g. MW, EP) are written as 0.0001, since
    LAMMPS requires a strictly positive mass."""
    with open(path, "w") as fh:
        fh.write("Masses\n\n")
        for row in nonbonded_matrix:
            mass_value = row["mass"]
            if float(mass_value) == 0.0:
                mass_value = "0.0001"
            fh.write(f"{row['id']} {mass_value}\n")
        fh.write("\n")


def write_pairs(path, nonbonded_matrix):
    """pairs.dat: 'Pair Coeffs' header, then 'id  epsilon/4.184  sigma*10'
    per line (kJ/mol -> kcal/mol, nm -> Angstrom)."""
    with open(path, "w") as fh:
        fh.write("Pair Coeffs\n\n")
        for row in nonbonded_matrix:
            epsilon_kcal = row["epsilon_kjmol"] / 4.184
            sigma_ang = row["sigma_nm"] * 10
            fh.write(f"{row['id']} {epsilon_kcal:g} {sigma_ang:g}\n")
        fh.write("\n")


def write_atoms(path, atoms, box_center_shift=None):
    """atom.dat: 'Atoms' header, then 7 columns per atom:
    atom_number  rigid_body_id  atom_type_id  charge  x  y  z

    box_center_shift is a (dx, dy, dz) tuple subtracted from every
    coordinate -- pass the structure's centroid to center it at the origin.
    """
    dx, dy, dz = box_center_shift if box_center_shift else (0, 0, 0)
    with open(path, "w") as fh:
        fh.write("Atoms\n\n")
        for atom in atoms:
            x, y, z = atom["x"] - dx, atom["y"] - dy, atom["z"] - dz
            fh.write(
                f"{atom['new_serial']} {atom['rigid_id']} {atom['atom_type_id']} "
                f"{atom['charge']:g}\t{x:8.3f}{y:9.3f}{z:9.3f}\n"
            )


def assemble_data_file(out_path, part_paths):
    """Concatenate metadata.dat, mass.dat, pairs.dat, atom.dat (in order)
    into the final output file, with a trailing blank line."""
    with open(out_path, "w") as out_fh:
        for part in part_paths:
            with open(part) as part_fh:
                out_fh.write(part_fh.read())
        out_fh.write("\n")
