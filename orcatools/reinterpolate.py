#!/usr/bin/env python3
import argparse
import math
import pathlib
import sys
from io import StringIO

import ase.io
import numpy as np
from ase import Atoms
from ase.io import read
from geodesic_interpolate.coord_utils import align_geom
from geodesic_interpolate.geodesic import Geodesic
from geodesic_interpolate.interpolation import redistribute


# Mathe
def gcd(a, b):
    '''Return the greatest common divisor using Euclid's Algorithm.'''
    while b:
        a, b = b, a % b
    return a

def lcm(a, b):
    '''Return lowest common multiple.'''
    return a * b // gcd(a, b)


def read_allxyz(filename):
    with open(filename, "r") as f:
        file = f.read()
    xyz = file.replace(">\n", "")
    file = StringIO(xyz)
    read_atoms = ase.io.read(file, index=':', format='xyz')
    return read_atoms

def center_all(structures):
    new_structures = []
    for structure in structures:
        structure.center()
        new_structures.append(structure)
    return new_structures

def write_allxyz(trajectory, filename):
    with StringIO() as file, open(filename, "w") as allxyz_file:
        file.name = "file.xyz"
        ase.io.write(file, trajectory, format="xyz")
        file.seek(0)
        first_frame = True
        while True:
            # Read number of atoms
            num_atoms_line = file.readline()
            if not num_atoms_line:
                break  # End of file

            # Add '>' separator if it's not the first frame
            if not first_frame:
                allxyz_file.write(">\n")
            first_frame = False  # Only the first frame skips the separator

            num_atoms = int(num_atoms_line.strip())

            # Read comment line (or metadata)
            comment_line = file.readline().strip()

            # Write the frame to the .allxyz file
            allxyz_file.write(f"{num_atoms}\n")
            allxyz_file.write(f"{comment_line}\n")

            # Read each atom's data and write it to the .allxyz file
            for _ in range(num_atoms):
                atom_data = file.readline()
                allxyz_file.write(atom_data)


def ase_geodesic_interpolate(initial_mol, final_mol, n_images=20, friction=0.01, dist_cutoff=3, scaling=1.7, sweep=None,
                             tol=0.002, maxiter=15, microiter=20):
    atom_string = initial_mol.symbols
    atoms = list(atom_string)
    initial_pos = [initial_mol.positions]
    final_pos = [final_mol.positions]
    total_pos = initial_pos + final_pos

    # First redistribute number of images.  Perform interpolation if too few and subsampling if too many
    # images are given
    raw = redistribute(atoms, total_pos, n_images, tol=tol * 5)

    # Perform smoothing by minimizing distance in Cartesian coordinates with redundant internal metric
    # to find the appropriate geodesic curve on the hyperspace.
    smoother = Geodesic(atoms, raw, scaling, threshold=dist_cutoff, friction=friction)

    if sweep is None:
        sweep = len(atoms) > 35
    if sweep:
        smoother.sweep(tol=tol, max_iter=maxiter, micro_iter=microiter)
    else:
        smoother.smooth(tol=tol, max_iter=maxiter)
    all_mols = []
    for pos in smoother.path:
        mol = Atoms(atom_string, pos)
        all_mols.append(mol)
    return all_mols


def rmsd(a, b):
    return np.sqrt(np.mean(np.square(a.positions - b.positions)))


def moved(atoms, positions):
    atoms = atoms.copy()
    atoms.positions = positions
    return atoms


def build_path(traj, n_new):
    """Join the structures in traj with n_new[i] geodesic frames between traj[i] and traj[i + 1].

    All frames share one frame of reference: traj[0] is centered, every later structure is rotated and
    translated (internal geometry unchanged) to continue the path without a rigid jump.
    """
    current = moved(traj[0], traj[0].positions - np.mean(traj[0].positions, axis=0))
    path = [current]
    for i, n in enumerate(n_new):
        if n > 0:
            # The geodesic path contains both end points, aligned frame by frame onto the start
            frames = ase_geodesic_interpolate(current, traj[i + 1], n_images=n + 2)
            if len(frames) != n + 2:
                raise RuntimeError(f"Geodesic interpolation returned {len(frames)} frames instead of {n + 2}.")
            if rmsd(frames[0], current) > 1e-6:
                raise RuntimeError("Geodesic interpolation moved the start structure.")
            path.extend(frames[1:-1])
            current = moved(traj[i + 1], frames[-1].positions)
        else:
            current = moved(traj[i + 1], align_geom(current.positions, traj[i + 1].positions)[1])
        path.append(current)
    return path


def select_by_distance(path, nimages):
    """Pick nimages frames from path at equal cumulative RMSD, plus both endpoints."""
    if len(path) < nimages + 2:
        raise RuntimeError(f"Path has {len(path)} frames, cannot select {nimages + 2}.")
    s = np.concatenate([[0.0], np.cumsum([rmsd(a, b) for a, b in zip(path, path[1:])])])
    if s[-1] == 0:
        raise ValueError("All structures are identical, the path has zero length.")
    indices = [0]
    for k in range(1, nimages + 1):
        nearest = int(np.argmin(np.abs(s - k * s[-1] / (nimages + 1))))
        # strictly increasing, and leave one frame for each remaining target and the end point
        indices.append(min(max(nearest, indices[-1] + 1), len(path) - 2 - (nimages - k)))
    indices.append(len(path) - 1)
    return [path[i] for i in indices]


def interpolation_core(traj, nimages=14, use_distance=False, verbose=False, oversample=8):
    """Return nimages + 2 frames from traj[0] to traj[-1] along a geodesic path through all structures in traj.

    Image mode: the path is split so that the output frames are separated by the same number of path frames;
    input structures are kept wherever the numbers allow it.
    Distance mode: a path with about oversample times more frames than the output is built, and the output
    frames are taken from it at equal cumulative RMSD. Interior input structures are not kept.
    All output frames share one frame of reference (see build_path).
    """
    if len(traj) < 2:
        raise ValueError(f"Interpolation needs at least two structures, got {len(traj)}.")
    num_segments = len(traj) - 1
    if use_distance:
        lengths = [align_geom(a.positions, b.positions)[0] for a, b in zip(traj, traj[1:])]
        total = sum(lengths)
        if total == 0:
            raise ValueError("All structures are identical, the path has zero length.")
        # intervals per segment proportional to its length, oversample intervals per output step
        n_new = [max(0, math.ceil(oversample * (nimages + 1) * length / total) - 1) for length in lengths]
        if verbose:
            print("Division by distance, new frames per segment:", n_new)
        return select_by_distance(build_path(traj, n_new), nimages)
    # Smallest number of path intervals that both the segments and the output steps divide
    n_intervals = lcm(num_segments, nimages + 1)
    per_segment = n_intervals // num_segments
    step = n_intervals // (nimages + 1)
    if verbose:
        print(f"Division by images, {per_segment - 1} new frames per segment, every {step}. frame is used")
    path = build_path(traj, [per_segment - 1] * num_segments)
    return path[::step]

def reinterpolate(images=14, output="traj.xyz", structures=None, trajectory=None, distance=False, verbose=False):
    if verbose:
        print("Reinterpolation tool V0.0.3 - 02.11.2024")
    # Validate that there is one input
    if structures is None and trajectory is None:
        print("Please specify at least one structure or trajectory.")
        sys.exit(1)
    if structures is not None and trajectory is not None:
        print("Cannot specify both --structures and --trajectory.")
        sys.exit(1)
    outpath = pathlib.Path(output)
    if outpath.suffix not in (".xyz", ".allxyz"):
        raise ValueError(f"Output file format {outpath.suffix} is not supported, use .xyz or .allxyz.")
    # Convert the input into a ase trajectory
    if structures is not None:
        if verbose:
            print("Creating Trajectory from input structures...")
        traj = []
        for structure in structures:
            atoms = read(structure)
            traj.append(atoms)
    else:
        if pathlib.Path(trajectory).suffix not in (".xyz", ".allxyz"):
            raise ValueError(f"Trajectory file format {pathlib.Path(trajectory).suffix} is not supported, "
                             "use .xyz or .allxyz.")
        if pathlib.Path(trajectory).suffix == ".xyz":
            if verbose:
                print("Reading Trajectory from input xyz file")
            traj = read(trajectory, index=":")
        else:
            if verbose:
                print("Reading Trajectory from input allxyz file")
            traj = read_allxyz(trajectory)
    # Perform Interpolation
    traj = center_all(traj)
    # No second center_all: it would shift every frame by its own bounding-box center. The frames are
    # already aligned with a common center.
    output_traj = interpolation_core(traj, nimages=images, verbose=verbose, use_distance=distance)
    print("n_output", len(output_traj))
    # Format output file
    if outpath.suffix == ".xyz":
        if verbose:
            print("Writing Trajectory to output xyz file")
        ase.io.write(outpath, output_traj, format="xyz")
    else:
        if verbose:
            print("Writing Trajectory to output allxyz file")
        write_allxyz(output_traj, outpath)
    return output_traj

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="reinterpolate.py",
        description="Re-interpolate trajectories or between structures using geodesic interpolation.",
        epilog=''
    )
    parser.add_argument("-n", "--images", type=int, help="Number of images to generate (excluding start and end)", required=True)
    parser.add_argument("-o", "--output", type=str, help="Output file, format can be xyz or allxyz", required=True)
    parser.add_argument("-s", "--structures", type=str, nargs='+', help="Structures to be included in the interpolation ordered from beginning to end", required=False, default=None)
    parser.add_argument("-t", "--trajectory", type=str, help="Trajectory to be re-interpolated, formatted as xyz or allxyz", required=False, default=None)
    parser.add_argument("-d", "--distance", action="store_true", help="Distribute images evenly using RMSD distance. This is recommended for NEB trajectories but will break NEB-CI trajectories", required=False, default=False)
    parser.add_argument("-v", "--verbose", action='store_true', help="Verbose output", required=False)
    args = parser.parse_args()
    reinterpolate(images=args.images, output=args.output, structures=args.structures, trajectory=args.trajectory, distance=args.distance, verbose=args.verbose)


