"""Parse rigid_definitions.dat (chain, resid range, rigid id, rigid type)
and apply rigid body IDs to atoms."""


def read_rigid_definitions(path):
    """Returns a list of dicts: chain, resid_lo, resid_hi, rigid_id, rigid_type."""
    defs = []
    with open(path) as fh:
        for line in fh:
            stripped = line.strip()
            if not stripped or stripped.startswith(";") or stripped.startswith("#"):
                continue
            tokens = stripped.split()
            if len(tokens) < 4:
                continue
            chain, resid_range, rigid_id, rigid_type = tokens[0], tokens[1], tokens[2], tokens[3]
            lo_str, hi_str = resid_range.split("-")
            defs.append({
                "chain": chain,
                "resid_lo": int(lo_str),
                "resid_hi": int(hi_str),
                "rigid_id": int(rigid_id),
                "rigid_type": int(rigid_type),
            })
    return defs


def apply_rigid_ids(atoms, rigid_defs):
    """Sets atom["rigid_id"] = 0 by default, then overwrites it for any
    atom whose (new_chain, resid) falls inside a rigid_definitions.dat row.
    Mutates and returns the same list.
    """
    for atom in atoms:
        atom["rigid_id"] = 0

    for d in rigid_defs:
        for atom in atoms:
            if (atom["new_chain"] == d["chain"]
                    and d["resid_lo"] <= atom["resid"] <= d["resid_hi"]):
                atom["rigid_id"] = d["rigid_id"]
    return atoms
