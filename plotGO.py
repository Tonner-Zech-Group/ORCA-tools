#!/usr/bin/env python3
#
# Script to plot ORCA Geometry Convergence Output
# by Patrick Melix
# 2023/05
#
# You can import the module and then call .main() or use it as a script
# Needs grep and tail.
import argparse, os, subprocess, sys, glob
import numpy as np
from matplotlib.ticker import MaxNLocator
from matplotlib import pyplot as plt
from ase.units import create_units

def plot(xAxis, iterations, labels, filename, lw=3, s=0):
    msbig = 9
    for i_label, label in enumerate(labels):
        data = np.array([ d[i_label][0] for d in iterations ])
        data_converged = np.array([ d[i_label][-1] for d in iterations ])
        if True in data_converged:
            not_converged_max_idx = min(np.where(data_converged == True)[0])
        else:
            not_converged_max_idx = len(data)
        ax = plt.figure().gca()
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        plt.xlabel('Step Number')
        plt.ylabel(label)
        plt.gca().yaxis.grid(True)
        plt.plot(xAxis, data, color='black', ls=':', lw=lw)
        plt.scatter(xAxis[0:not_converged_max_idx], data[0:not_converged_max_idx], marker='x', color='red', s=(msbig+s)**2)
        plt.scatter(xAxis[not_converged_max_idx:], data[not_converged_max_idx:], marker='x', color='green', s=(msbig+s)**2)
        #plt.legend()
        plt.tight_layout()
        plt.savefig("{}_{}.png".format(filename,label.replace(' ', '_')))
        #plt.show()
        plt.close()


def main(filename='convergence.png', presentation=False):
    input_file = [ f for f in glob.glob('*.inp') ]
    if len(input_file) > 1:
        raise ValueError("More than one input file found.")
    elif len(input_file) == 0:
        raise ValueError("No inp file found.")

    # get output file
    output_file = input_file[0].replace('.inp', '.out')
    assert os.path.isfile(output_file), "Output file {:} does not exist.".format(output_file)
    print("Reading output file {:}".format(output_file))
    iterations = []
    n_iterations = 0
    labels = []
    with open(output_file, 'r') as f:
        line = f.readline()
        while line:
            if '|Geometry convergence|' in line:
                iterations.append([])
                n_iterations += 1
                line = f.readline() #skip header
                line = f.readline() #skip --- line
                line = f.readline() #read first data line
                while line:
                    tmp = line.split()
                    if n_iterations == 1:
                        labels.append(" ".join(reversed(tmp[-4::-1])))
                    iterations[-1].append([])
                    iterations[-1][-1].append(float(tmp[-3])) #value
                    iterations[-1][-1].append(float(tmp[-2])) #tolerance
                    iterations[-1][-1].append(bool(tmp[-1] == "YES")) #converged
                    line = f.readline().strip()
                    if line == len(line) * line[0]: #only dots
                        break
            line = f.readline()
    print("Found the following labels: {}".format(", ".join(labels)))

    if presentation:
        lw = 5
        s = 3
        plt.rcParams.update({'font.size': 22})
        plt.rcParams.update({'legend.fontsize': 22})
    else:
        lw = 3
        s = 0

    plot(np.arange(1, n_iterations + 1), iterations, labels, filename=filename, lw=lw, s=s)


if __name__ == "__main__":
    f = "/home/patrickm/git/Python4ChemistryTools/mpl-settings.py"
    if os.path.isfile(f):
        exec(open(f).read())
    parser = argparse.ArgumentParser(description='Plot ORCA Geometry Optimization Convergence')
    parser.add_argument('--file', help='Plot Filename Beginning, will be appended with _<label>.png', default='convergence')
    parser.add_argument('--presentation', help='Presentation Mode (i.e. thicker lines)', action='store_true')
    args = parser.parse_args()
    main(args.file, args.presentation)
