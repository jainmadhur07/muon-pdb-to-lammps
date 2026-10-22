"""Read ffbonded.itp -> bonds / angles / impropers / dihedrals matrices (lists
of dicts), scaled for the -temperature flexible model."""
import math
import os

SECTIONS = {
    "bondtypes": ("bonds", 2),
    "angletypes": ("angles", 3),
    "dihedraltypes": ("dihedrals", 4),
}


def _temperature_factor(temperature):
    return (1 + math.log(temperature / 300)) ** 2


def _read_sections(ff_dir):
    """Returns {"bonds", "angles", "impropers", "dihedrals"} -> raw
    (atom_types, numeric_tokens) rows."""
    data = {"bonds": [], "angles": [], "impropers": [], "dihedrals": []}
    current = None
    with open(os.path.join(ff_dir, "ffbonded.itp")) as fh:
        for raw in fh:
            stripped = raw.split(";", 1)[0].strip()
            if raw.lstrip().startswith("["):
                name = stripped.strip("[] ").strip()
                if name == "dihedraltypes" and "improper" in raw:
                    current = ("impropers", 4)
                else:
                    current = SECTIONS.get(name)
                continue
            if current is None or not stripped:
                continue
            key, n_atoms = current
            tokens = stripped.split()
            data[key].append((tokens[:n_atoms], tokens[n_atoms + 1:]))  # skip func
    return data


def read_bonded_matrices(ff_dir, temperature):
    """Returns (bonds, angles, impropers, dihedrals): lists of dicts, IDs 1..N."""
    raw = _read_sections(ff_dir)
    factor = _temperature_factor(temperature)

    bonds = []
    for i, (types, vals) in enumerate(raw["bonds"], start=1):
        b0, kb = float(vals[0]), float(vals[1])
        bonds.append({"id": i, "types": types,
                      "kb": (kb / 836.8) * factor, "b0": b0 * 10})

    angles = []
    for i, (types, vals) in enumerate(raw["angles"], start=1):
        th0, cth = float(vals[0]), float(vals[1])
        angles.append({"id": i, "types": types,
                       "cth": (cth / 8.368) * factor, "th0": th0})

    impropers = []
    for i, (types, vals) in enumerate(raw["impropers"], start=1):
        d, k, n = float(vals[0]), float(vals[1]), int(vals[2])
        if d == 0:
            d = 1
        elif d == 180:
            d = -1
        else:
            raise ValueError(
                f"unknown improper type is observed, set the parameters: {types} phase={d}")
        k = factor if k == 0 else (k / 4.184) * factor
        impropers.append({"id": i, "types": types, "k": k, "d": d, "n": n})

    dihedrals = []
    for i, (types, vals) in enumerate(raw["dihedrals"], start=1):
        d, k, n = float(vals[0]), float(vals[1]), int(vals[2])
        k = factor if k == 0 else (k / 4.184) * factor
        dihedrals.append({"id": i, "types": types, "k": k, "n": n,
                          "d": int(math.floor(d + 0.5)), "weight": 1})

    return bonds, angles, impropers, dihedrals
