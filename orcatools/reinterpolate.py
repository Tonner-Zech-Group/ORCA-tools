#!/home/franzthiemann/ThC/bin/vnev_arch/bin/python
import math
import os
import pathlib
import tempfile
import numpy as np
from ase.io import read
import ase.io
from ase import Atoms
from geodesic_interpolate.interpolation import redistribute
from geodesic_interpolate.geodesic import Geodesic
import argparse

# Mathe
def gcd(a, b):
    '''Return the greatest common divisor using Euclid's Algorithm.'''
    while b:
        a, b = b, a % b
    return a

def lcm(a, b):
    '''Return lowest common multiple.'''
    return a * b / gcd(a, b)


def read_allxyz(filename):
    with open(filename, "r") as f:
        file = f.read()
    xyz = file.replace(">\n", "")
    tmp_dir = tempfile.gettempdir()
    file = os.path.join(tmp_dir, "structure.xyz")
    f = open(file, "w")
    f.write(xyz)
    f.close()
    read_atoms = ase.io.read(file, index=':')
    return read_atoms

def center_all(structures):
    new_structures = []
    for structure in structures:
        structure.center()
        new_structures.append(structure)
    return new_structures

def convert_xyz(source, destination):
    with open(source, "r") as xyz_file, open(destination, "w") as allxyz_file:
        first_frame = True
        while True:
            # Read number of atoms
            num_atoms_line = xyz_file.readline()
            if not num_atoms_line:
                break  # End of file

            # Add '>' separator if it's not the first frame
            if not first_frame:
                allxyz_file.write(">\n")
            first_frame = False  # Only the first frame skips the separator

            num_atoms = int(num_atoms_line.strip())

            # Read comment line (or metadata)
            comment_line = xyz_file.readline().strip()

            # Write the frame to the .allxyz file
            allxyz_file.write(f"{num_atoms}\n")
            allxyz_file.write(f"{comment_line}\n")

            # Read each atom's data and write it to the .allxyz file
            for _ in range(num_atoms):
                atom_data = xyz_file.readline()
                allxyz_file.write(atom_data)


def write_allxyz(trajectory, filename):
    tmp_dir = tempfile.gettempdir()
    file = os.path.join(tmp_dir, "output.xyz")
    ase.io.write(file, trajectory, format="xyz")
    convert_xyz(file, filename)


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
    try:
        if sweep:
            smoother.sweep(tol=tol, max_iter=maxiter, micro_iter=microiter)
        else:
            smoother.smooth(tol=tol, max_iter=maxiter)
    finally:
        all_mols = []
        for pos in smoother.path:
            mol = Atoms(atom_string, pos)
            all_mols.append(mol)
        return all_mols


def interpolation_core(traj, nimages=14, use_distance=False, verbose=False):
    # Basic idea: Interpolate between each image, combine the path
    # Path is then split into equidistant images or into images seperated by the same number of invisible images
    nimages_tot = nimages + 2
    if nimages_tot == len(traj):
        return traj
    num_segments = len(traj) - 1
    # Calculate smallest multiplier between nimages and segments to get total target images
    # Efficient
    n_target = lcm(num_segments, (nimages_tot - 1)) + 1
    # Dumb
    # n_target = ((nimages_tot - 1) * num_segments) + 1
    if verbose:
        print("n_target", n_target) ## Debug

    n_generated = n_target - 2 - (len(traj) - 2)
    if verbose:
        print("n_generated", n_generated)
    n_per_segment = n_generated / num_segments
    if verbose:
        print("n_per_segment", n_per_segment)
    total_trj = [traj[0]]
    for i in range(num_segments):
        if verbose:
            print(f"Segment: {i+ 1} of {num_segments}")
        geo_frames = ase_geodesic_interpolate(traj[i], traj[i + 1], n_images=n_per_segment + 1)
        # print("geo_frames", len(geo_frames))
        total_trj = total_trj + geo_frames[1:] + [traj[i + 1]]
    print("Expected: ", n_target, len(total_trj))
    if not use_distance:
        if verbose:
            print("Devision by images")
        ret_images = []
        ret_images.append(total_trj[0])
        counter = 0
        if verbose:
            print((n_target - 1), "/", (nimages + 1))
        step = int((n_target - 1) / (nimages + 1))
        if verbose:
            print("step", step)
        for i in range(nimages):
            counter += step
            if verbose:
                print("counter", counter)
            ret_images.append(total_trj[counter])
        ret_images.append(total_trj[-1])
        return ret_images
    else:
        # calculate RMSD distances between all the images
        if verbose:
            print("Devision by distance")
        rmsds = [0]
        srmsds = [0]
        sum = 0
        for i in range(int(n_target) - 1):
            step = total_trj[i].positions - total_trj[i + 1].positions
            rms = np.sqrt(np.mean(np.square(step)))
            rmsds.append(rms)
            sum += rms
            srmsds.append(sum)
        rmsd_per_frame = sum / (nimages + 1)
        ret_images = []
        if verbose:
            print("rmsds", srmsds)
        if verbose:
            print("length", len(rmsds))
        target_rmsd = []
        for i in range(int(n_target) - 1):
            target_rmsd.append(rmsd_per_frame * i)
        if verbose:
            print("target_rmsd", target_rmsd)
        if verbose:
            print("length", len(srmsds), len(total_trj))
        i = 0
        target_rmsd.append(math.inf)
        while i in range(int(n_target)):
            image = total_trj[i]
            if srmsds[i] >= target_rmsd[0]:
                ret_images.append(image)
                target_rmsd = target_rmsd[1:]
            else:
                i = i + 1
        ret_images.append(total_trj[-1])
        return ret_images

def reinterpolate(images=14, output="traj.xyz", structures=None, trajectory=None, distance=False, verbose=False):
    if verbose:
        print("Reinterpolation tool V0.0.2 - 25.10.2024")
    # Validate that there is one input
    if structures is None and trajectory is None:
        print("Please specify at least one structure or trajectory.")
        exit(1)
    if structures is not None and trajectory is not None:
        print("Cannot specify both --structures and --trajectory.")
        exit(1)
    # Convert the input into a ase trajectory
    if structures is not None:
        if verbose:
            print("Creating Trajectory from input structures...")
        traj = []
        for structure in structures:
            atoms = read(structure)
            traj.append(atoms)
    else:
        if os.path.splitext(pathlib.Path(trajectory))[1] == ".xyz":
            if verbose:
                print("Reading Trajectory from input xyz file")
            traj = read(trajectory, index=":")
        elif os.path.splitext(pathlib.Path(trajectory))[1] == ".allxyz":
            if verbose:
                print("Reading Trajectory from input allxyz file")
            traj = read_allxyz(trajectory)
        else:
            print(f"File format {os.path.splitext(pathlib.Path(trajectory))[1]} is not supported.")
    # Perform Interpolation
    traj = center_all(traj)
    output_traj = interpolation_core(traj, nimages=images, verbose=verbose, use_distance=distance)
    output_traj = center_all(output_traj)
    print("n_output", len(output_traj))
    # Format output file
    outpath = pathlib.Path(output)
    if os.path.splitext(outpath)[1] == ".xyz":
        if verbose:
            print("Writing Trajectory to output xyz file")
        ase.io.write(outpath, output_traj, format="xyz")
        return output_traj
    elif os.path.splitext(outpath)[1] == ".allxyz":
        if verbose:
            print("Writing Trajectory to output allxyz file")
        write_allxyz(output_traj, outpath)
        return output_traj
    else:
        print(f"File format {os.path.splitext(outpath)[1]} is not supported.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog="reinterpolate.py",
        description="Re-interpolate trajectories or between structures using geodesic interpolation.",
        epilog=''
    )
    parser.add_argument("-n", "--images", type=int, help="Number of images to genereate (excluding start and end)", required=True)
    parser.add_argument("-o", "--output", type=str, help="Output file, format can be xyz or allxyz", required=True)
    parser.add_argument("-s", "--structures", type=str, nargs='+', help="Structurs to be included in the interpolation ordered from beginning to end", required=False, default=None)
    parser.add_argument("-t", "--trajectory", type=str, help="Trajectory to be re-interpolated, formated as xyz or allxyz", required=False, default=None)
    parser.add_argument("-d", "--distance", action="store_true", help="Distribute images evenly using RMSD distance. This is recommended for NEB trajectories but will break NEB-CI trajectories", required=False, default=False)
    parser.add_argument("-v", "--verbose", action='store_true', help="Verbose output", required=False)
    args = parser.parse_args()
    reinterpolate(images=args.images, output=args.output, structures=args.structures, trajectory=args.trajectory, distance=args.distance, verbose=args.verbose)


