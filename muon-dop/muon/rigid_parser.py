"""Parse rigid_definitions.dat (chain, resid range, rigid id, rigid type)
and apply rigid body IDs to atoms. The rigid id column can be a single
number (whole range shares one rigid body) or a range like "2-6"
(sequential id per residue, e.g. for a run of individual ions)."""


def _parse_range(text):
    """'1-110' -> (1, 110). A bare number '1' -> (1, 1)."""
    if "-" in text:
        lo_str, hi_str = text.split("-")
        return int(lo_str), int(hi_str)
    value = int(text)
    return value, value


def read_rigid_definitions(path):
    """Returns a list of dicts: chain, resid_lo, resid_hi, rigid_id_lo,
    rigid_id_hi, rigid_type. rigid_id_lo == rigid_id_hi for a plain
    (non-range) rigid id.
    """
    defs = []
    with open(path) as fh:
        for line in fh:
            stripped = line.strip()
            if not stripped or stripped.startswith(";") or stripped.startswith("#"):
                continue
            tokens = stripped.split()
            if len(tokens) < 4:
                continue
            chain, resid_range, rigid_id_field, rigid_type = tokens[0], tokens[1], tokens[2], tokens[3]
            resid_lo, resid_hi = _parse_range(resid_range)
            rigid_id_lo, rigid_id_hi = _parse_range(rigid_id_field)
            defs.append({
                "chain": chain,
                "resid_lo": resid_lo,
                "resid_hi": resid_hi,
                "rigid_id_lo": rigid_id_lo,
                "rigid_id_hi": rigid_id_hi,
                "rigid_type": int(rigid_type),
            })
    return defs


def apply_rigid_ids(atoms, rigid_defs):
    """Sets atom["rigid_id"] for every atom matching a rigid_definitions.dat
    row (default 0). Mutates and returns the same list."""
    for atom in atoms:
        atom["rigid_id"] = 0

    for d in rigid_defs:
        if d["rigid_id_lo"] == d["rigid_id_hi"]:
            # Plain case: one fixed id for the whole residue range.
            for atom in atoms:
                if (atom["new_chain"] == d["chain"]
                        and d["resid_lo"] <= atom["resid"] <= d["resid_hi"]):
                    atom["rigid_id"] = d["rigid_id_lo"]
            continue

        # Range case: one rigid id per residue, assigned in order.
        residues_in_range = sorted({
            atom["resid"] for atom in atoms
            if atom["new_chain"] == d["chain"]
            and d["resid_lo"] <= atom["resid"] <= d["resid_hi"]
        })
        expected_count = d["rigid_id_hi"] - d["rigid_id_lo"] + 1
        if len(residues_in_range) != expected_count:
            raise ValueError(
                f"rigid_definitions.dat row for chain {d['chain']} residues "
                f"{d['resid_lo']}-{d['resid_hi']} expects {expected_count} "
                f"residues (rigid id range {d['rigid_id_lo']}-{d['rigid_id_hi']}) "
                f"but found {len(residues_in_range)} residues in that range."
            )
        rigid_id_for_resid = {
            resid: d["rigid_id_lo"] + i
            for i, resid in enumerate(residues_in_range)
        }
        for atom in atoms:
            if atom["new_chain"] == d["chain"] and atom["resid"] in rigid_id_for_resid:
                atom["rigid_id"] = rigid_id_for_resid[atom["resid"]]

    return atoms
