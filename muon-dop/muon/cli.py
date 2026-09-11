"""Command-line argument parsing for MUON.py."""
import argparse


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        prog="MUON.py",
        description="Convert a PDB/CIF structure + force field into a LAMMPS data file.",
    )
    parser.add_argument(
        "-pdb", required=True, dest="structure_file",
        help="Input structure file in PDB or CIF format (e.g. new2.pdb)",
    )
    parser.add_argument(
        "-rigid", required=True, dest="rigid_file",
        help="Rigid body definitions file (e.g. rigid_definitions.dat)",
    )
    parser.add_argument(
        "-ff", required=True, dest="ff_dir",
        help="Path to the force field directory (contains .rtp / .itp files)",
    )
    parser.add_argument(
        "-box", required=True, nargs=3, type=float, metavar=("X", "Y", "Z"),
        dest="box",
        help="Simulation box dimensions, e.g. -box 150 150 150",
    )
    parser.add_argument(
        "-out", default="data_file", dest="out_file",
        help="Name of the final LAMMPS data file to write (default: data_file)",
    )
    parser.add_argument(
        "-log", default="out.log", dest="log_file",
        help="Name of the log file to write (default: out.log)",
    )
    return parser.parse_args(argv)
