"""Read a PDB/CIF file into a uniform "atoms matrix" (list of atom dicts)."""
from collections import defaultdict


def read_pdb(path):
    """Parse a standard fixed-column PDB file into a list of atom dicts."""
    atoms = []
    with open(path) as fh:
        for line in fh:
            record = line[0:6].strip()
            if record not in ("ATOM", "HETATM"):
                continue
            atoms.append({
                "serial": int(line[6:11]),
                "name": line[12:16].strip(),
                "resname": line[17:20].strip(),
                "chain": line[21:22].strip(),
                "resid": int(line[22:26]),
                "x": float(line[30:38]),
                "y": float(line[38:46]),
                "z": float(line[46:54]),
            })
    return atoms


def read_cif(path):
    """Parse the ATOM/HETATM rows of a minimal mmCIF file."""
    atoms = []
    with open(path) as fh:
        lines = fh.readlines()

    header = []
    in_loop = False
    data_start = None
    for i, line in enumerate(lines):
        stripped = line.strip()
        if stripped.startswith("_atom_site."):
            in_loop = True
            header.append(stripped.split(".", 1)[1])
            continue
        if in_loop and not stripped.startswith("_") and stripped:
            data_start = i
            break

    if data_start is None:
        raise ValueError(f"No _atom_site loop found in CIF file: {path}")

    col = {name: idx for idx, name in enumerate(header)}

    def get(tokens, *names, cast=str, default=None):
        for n in names:
            if n in col:
                return cast(tokens[col[n]])
        if default is not None:
            return default
        raise KeyError(f"None of {names} found in CIF _atom_site header")

    for line in lines[data_start:]:
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        if stripped.startswith("_") or stripped.startswith("loop_"):
            break
        tokens = stripped.split()
        if tokens[0] not in ("ATOM", "HETATM"):
            continue
        atoms.append({
            "serial": get(tokens, "id", cast=int),
            "name": get(tokens, "label_atom_id", cast=lambda v: v.strip('"')),
            "resname": get(tokens, "label_comp_id"),
            "chain": get(tokens, "auth_asym_id", "label_asym_id"),
            "resid": get(tokens, "auth_seq_id", "label_seq_id", cast=int),
            "x": get(tokens, "Cartn_x", cast=float),
            "y": get(tokens, "Cartn_y", cast=float),
            "z": get(tokens, "Cartn_z", cast=float),
        })
    return atoms


def read_structure(path):
    """Dispatch to read_pdb or read_cif based on file extension."""
    lower = path.lower()
    if lower.endswith(".cif"):
        return read_cif(path)
    if lower.endswith(".pdb"):
        return read_pdb(path)
    raise ValueError(f"Unrecognized structure file extension: {path}")


def group_residues(atoms):
    """Group atoms into residues keyed by (new_chain, resid), in order of
    first appearance. Call after relabel_chains_and_atoms."""
    residues = {}
    for atom in atoms:
        key = (atom["new_chain"], atom["resid"])
        residues.setdefault(key, []).append(atom)
    return residues


def relabel_chains_and_atoms(atoms):
    """Renumber atoms 1..N and relabel chains with an incremental numeric
    suffix (A, B, A, B, A, C -> A1, B1, A2, B2, A3, C1). Also tags each atom
    with is_nterm / is_cterm (first/last residue of its own chain run)."""
    run_counts = defaultdict(int)
    prev_chain = object()  # sentinel, never equal to a real chain id
    current_new_chain = None

    for atom in atoms:
        if atom["chain"] != prev_chain:
            run_counts[atom["chain"]] += 1
            current_new_chain = f"{atom['chain']}{run_counts[atom['chain']]}"
            prev_chain = atom["chain"]
        atom["new_chain"] = current_new_chain

    for i, atom in enumerate(atoms, start=1):
        atom["new_serial"] = i

    resid_by_new_chain = defaultdict(set)
    for atom in atoms:
        resid_by_new_chain[atom["new_chain"]].add(atom["resid"])

    bounds = {
        chain: (min(resids), max(resids))
        for chain, resids in resid_by_new_chain.items()
    }

    for atom in atoms:
        lo, hi = bounds[atom["new_chain"]]
        atom["is_nterm"] = (atom["resid"] == lo)
        atom["is_cterm"] = (atom["resid"] == hi)

    return atoms
