"""Read the force field directory: ffnonbonded.itp -> nonbonded matrix,
.rtp/.r2b files -> per-atom charge and force field type lookup, and
ions.itp -> the list of ion residue names (for -ions padding)."""
import os


def read_nonbonded_matrix(ff_dir):
    """Parse ffnonbonded.itp -> list of dicts (id, name, mass, sigma_nm,
    epsilon_kjmol), one per atom type, numbered 1..N in file order."""
    path = os.path.join(ff_dir, "ffnonbonded.itp")
    matrix = []
    next_id = 1
    with open(path) as fh:
        for line in fh:
            stripped = line.split(";", 1)[0].strip()
            if not stripped or stripped.startswith("["):
                continue
            tokens = stripped.split()
            # atomtypes data line: name  at.num  mass  charge  ptype  sigma  epsilon
            if len(tokens) < 7:
                continue
            try:
                float(tokens[2])  # validate it's a real data line, not a header
                sigma = float(tokens[-2])
                epsilon = float(tokens[-1])
            except ValueError:
                continue
            matrix.append({
                "id": next_id,
                "name": tokens[0],
                "mass": tokens[2],  # kept as original string for exact formatting
                "sigma_nm": sigma,
                "epsilon_kjmol": epsilon,
            })
            next_id += 1
    return matrix


def read_ion_resnames(ff_dir):
    """Parse ions.itp -> list of ion residue names (e.g. IB+, CA, CL, NA,
    MG, K, RB, CS, LI, ZN), in file order. Each ion is a GROMACS
    [ moleculetype ] block; its name is the molname on the first real
    data line after that header.
    """
    path = os.path.join(ff_dir, "ions.itp")
    names = []
    with open(path) as fh:
        lines = fh.readlines()
    for i, raw in enumerate(lines):
        stripped = raw.split(";", 1)[0].strip()
        if stripped != "[ moleculetype ]":
            continue
        for j in range(i + 1, len(lines)):
            data_line = lines[j].split(";", 1)[0].strip()
            if not data_line:
                continue
            names.append(data_line.split()[0])
            break
    return names


def _read_r2b(path):
    """Parse a .r2b file into {gmx_resname: {"main":.., "nter":.., "cter":.., "ter2":..}}."""
    mapping = {}
    if not os.path.exists(path):
        return mapping
    with open(path) as fh:
        for line in fh:
            stripped = line.strip()
            if not stripped or stripped.startswith(";"):
                continue
            tokens = stripped.split()
            if len(tokens) < 4:
                continue
            gmx, main, nter, cter = tokens[0], tokens[1], tokens[2], tokens[3]
            ter2 = tokens[4] if len(tokens) > 4 else "-"
            mapping[gmx] = {"main": main, "nter": nter, "cter": cter, "ter2": ter2}
    return mapping


def _read_rtp_blocks(path):
    """Parse an .rtp file into {block_name: {atom_name: (fftype, charge)}}."""
    blocks = {}
    if not os.path.exists(path):
        return blocks
    current_block = None
    in_atoms_section = False
    with open(path) as fh:
        for raw in fh:
            line = raw.split(";", 1)[0].strip()  # strip inline comments too
            if not line:
                continue
            if line.startswith("[") and line.endswith("]"):
                name = line[1:-1].strip()
                if name == "atoms":
                    in_atoms_section = True
                    continue
                in_atoms_section = False
                if name == "bondedtypes":
                    current_block = None
                    continue
                current_block = name
                blocks[current_block] = {}
                continue
            if in_atoms_section and current_block is not None:
                tokens = line.split()
                if len(tokens) >= 3:
                    atom_name, fftype, charge = tokens[0], tokens[1], tokens[2]
                    blocks[current_block][atom_name] = (fftype, float(charge))
    return blocks


class ResidueDatabase:
    """Loads aminoacids/dna/rna .rtp + .r2b files and resolves
    (resname, atomname, is_nterm, is_cterm) -> (fftype, charge)."""

    FAMILIES = ["aminoacids", "dna", "rna"]
    HISTIDINE_TRY_ORDER = ["HISD", "HISE", "HISH"]  # try HID, then HIE, then HIP

    def __init__(self, ff_dir):
        self.ff_dir = ff_dir
        self.r2b = {}
        self.rtp = {}
        for family in self.FAMILIES:
            self.r2b[family] = _read_r2b(os.path.join(ff_dir, f"{family}.r2b"))
            self.rtp[family] = _read_rtp_blocks(os.path.join(ff_dir, f"{family}.rtp"))

    def _resolve_block_name(self, family, resname, is_nterm, is_cterm):
        entry = self.r2b[family].get(resname)
        if entry is None:
            return resname  # e.g. ions/water with no N-ter/C-ter variants
        if is_nterm and entry["nter"] != "-":
            return entry["nter"]
        if is_cterm and entry["cter"] != "-":
            return entry["cter"]
        return entry["main"] if entry["main"] != "-" else resname

    def lookup(self, resname, atomname, is_nterm, is_cterm):
        """Returns (fftype, charge, family) or raises LookupError."""
        for family in self.FAMILIES:
            block_name = self._resolve_block_name(family, resname, is_nterm, is_cterm)
            block = self.rtp[family].get(block_name)
            if block is None:
                continue
            if atomname in block:
                fftype, charge = block[atomname]
                return fftype, charge, family
        raise LookupError(
            f"Could not find atom '{atomname}' of residue '{resname}' "
            f"(is_nterm={is_nterm}, is_cterm={is_cterm}) in any of "
            f"{self.FAMILIES}'s .rtp files under {self.ff_dir}"
        )

    def resolve_histidine_block(self, atom_names, is_nterm, is_cterm):
        """For a histidine residue, try HID first; if its atoms don't all
        fit, try HIE; then HIP (terminal residues use the N-ter/C-ter
        variant of each, e.g. NHID/CHID). Returns (block_name, family)."""
        family = "aminoacids"
        tried = []
        for gmx_name in self.HISTIDINE_TRY_ORDER:
            block_name = self._resolve_block_name(family, gmx_name, is_nterm, is_cterm)
            block = self.rtp[family].get(block_name)
            tried.append(block_name)
            if block is None:
                continue
            if atom_names.issubset(block.keys()):
                return block_name, family
        raise LookupError(
            f"Histidine atom set {sorted(atom_names)} did not match any of "
            f"the three candidate blocks tried (in order): {tried}. "
            f"(is_nterm={is_nterm}, is_cterm={is_cterm})"
        )

    def ion_fftype(self, ion_resname):
        """Returns the force field atom type (e.g. 'K', 'Na', 'C0') for an
        ion residue, by looking up its one atom in aminoacids.rtp (ions are
        defined there, not in a separate .rtp file). Each ion's rtp block
        has exactly one atom -- usually named the same as the residue
        (e.g. [ K ] contains atom K), but not always (e.g. [ IB+ ]
        contains atom IB), so this takes whichever single atom is in the
        block rather than assuming the atom name matches the residue name.
        """
        block = self.rtp["aminoacids"].get(ion_resname)
        if not block:
            raise LookupError(
                f"Ion residue '{ion_resname}' not found as its own block "
                f"in aminoacids.rtp"
            )
        if len(block) != 1:
            raise LookupError(
                f"Expected ion residue '{ion_resname}' to have exactly one "
                f"atom in aminoacids.rtp, found {len(block)}: {sorted(block)}"
            )
        (fftype, _charge), = block.values()
        return fftype


def build_ion_vdw_radii(ion_resnames, resdb, nonbonded_matrix):
    """For each ion residue name, find its force field atom type (via
    ResidueDatabase.ion_fftype), then that type's sigma in the nonbonded
    matrix, and compute its van der Waals radius = (sigma_nm * 10) / 2
    (nm -> Angstrom, diameter -> radius). Returns {ion_resname: radius_ang}.
    """
    sigma_by_fftype = {row["name"]: row["sigma_nm"] for row in nonbonded_matrix}
    radii = {}
    for ion_resname in ion_resnames:
        fftype = resdb.ion_fftype(ion_resname)
        if fftype not in sigma_by_fftype:
            raise LookupError(
                f"Ion '{ion_resname}' has force field type '{fftype}', "
                f"which was not found in ffnonbonded.itp"
            )
        sigma_nm = sigma_by_fftype[fftype]
        radii[ion_resname] = (sigma_nm * 10) / 2
    return radii
