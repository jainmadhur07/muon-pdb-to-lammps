"""Read the force field directory: ffnonbonded.itp -> nonbonded matrix, and
.rtp/.r2b files -> per-atom charge and force field type lookup."""
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
                float(tokens[2])
                sigma = float(tokens[-2])
                epsilon = float(tokens[-1])
            except ValueError:
                continue
            matrix.append({
                "id": next_id,
                "name": tokens[0],
                "mass": tokens[2], 
                "sigma_nm": sigma,
                "epsilon_kjmol": epsilon,
            })
            next_id += 1
    return matrix


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
            line = raw.split(";", 1)[0].strip()
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
            return resname
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
