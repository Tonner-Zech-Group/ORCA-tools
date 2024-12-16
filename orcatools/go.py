import argparse
import os
import glob
import re

import numpy as np
from matplotlib.ticker import MaxNLocator
from matplotlib import pyplot as plt
import ase.units as units

def make_go_plot(xAxis, iterations, labels, filename, filelabels=None, lw=3, s=0, show=False):
    msbig = 9
    if filelabels is None:
        filelabels = labels
    for i_label, label in enumerate(labels):
        data = np.array([ d[i_label][0] for d in iterations ])
        data_converged = np.array([ d[i_label][-1] for d in iterations ])
        converged_indices = np.where(data_converged)[0]
        not_converged_indices = np.where(~data_converged)[0]
        ax = plt.figure().gca()
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        plt.xlabel('Step Number')
        plt.ylabel(label)
        plt.gca().yaxis.grid(True)
        plt.plot(xAxis, data, color='black', ls=':', lw=lw)
        plt.scatter(xAxis[not_converged_indices], data[not_converged_indices], marker='x', color='red', s=(msbig+s)**2)
        plt.scatter(xAxis[converged_indices], data[converged_indices], marker='x', color='green', s=(msbig+s)**2)
        plt.tight_layout()
        plt.savefig("{}_{}.png".format(filename,filelabels[i_label].replace(' ', '_')))
        if show:
            plt.show()
        plt.close()


def read_output(output_file):
    with open(output_file, 'r') as f:
        file = f.read()
        get_lines = r"[-]+\|Geometry convergence\|[-]+\n(?:.*\n)*?[-|\s]+\n((?:\s*\S+\s\S+\s+-?\d+\.\d+\s+\d+\.\d+\s+(?:YES|NO)\s*\n)+)"
        get_values = r"\s*(\S+\s\S+)\s+(-?\d+\.\d+)\s+(\d+\.\d+)\s+(YES|NO)"
        sections = re.findall(get_lines, file)
        iterations = []
        labels = []
        # Very inefficient dual pass method, but I do not know a safe way to avoid this
        # Pass 1 get Labels
        for section in sections:
            values = re.findall(get_values, section)
            for value in values:
                if value[0] not in labels:
                    labels.append(value[0])
        # second pass to fill the data array with values
        for section in sections:
            frame = []
            values = re.findall(get_values, section)
            for label in labels:
                found = False
                for value in values:
                    if label == value[0]:
                        frame.append([float(value[1]), float(value[2]), bool(value[3] == "YES")])
                        found = True
                        break
                if not found:
                    frame.append([float("NaN"), float("NaN"), False])
            iterations.append(frame)
    return iterations, labels

def convert_units(data, labels, target="au"):
    unit_labels = labels.copy()
    # print(data)
    if target == "au":
        for i in range(len(unit_labels)):
            if unit_labels[i] == "Energy change":
                unit_labels[i] = "Energy change / $E_h$"
            if "gradient" in unit_labels[i]:
                unit_labels[i] = unit_labels[i] + " / $E_h\\, a_0^{-1}$"
            if "step" in unit_labels[i]:
                unit_labels[i] = unit_labels[i] + " / $a_0$"
    if target == "si":
        for i in range(len(unit_labels)):
            if unit_labels[i] == "Energy change":
                unit_labels[i] = "Energy change / $\\mathrm{kJ}\\, \\mathrm{mol}^{-1}$"
                data = convert_data(data, i, units.Hartree * units.mol / units.kJ)
            if "gradient" in unit_labels[i]:
                unit_labels[i] = unit_labels[i] + " / $\\mathrm{kJ}\\, \\mathrm{mol}^{-1}\\, \\mathrm{nm}^{-1}$"
                data = convert_data(data, i, (units.Hartree * units.mol / units.kJ) * (units.nm / units.Bohr))
            if "step" in unit_labels[i]:
                unit_labels[i] = unit_labels[i] + " / $\\mathrm{nm}$"
                data = convert_data(data, i, units.Bohr / units.nm)
    return data, unit_labels

def convert_data(data, index, factor):
    for i in range(len(data)):
        pass
        data[i][index][0] = data[i][index][0] * factor
    return data


def plot_orca_go(filename='convergence.png', presentation=False, path='.', show=False, unit="si"):
    # check for numerical subfolders
    subfolders = [ int(f) for f in glob.glob(os.path.join(path,'*')) if os.path.isdir(f) and f.isdigit() ]
    subfolders.sort()
    subfolders = [ str(d) for d in subfolders ]
    print("Found {} numerical subfolders, iterating over them and root:".format(len(subfolders)), end=' ')
    subfolders.append(os.path.join(path, '.'))
    iterations = []
    labels = []
    for dir in subfolders:
        if dir == '.':
            print("root", end=', ')
        else:
            print(dir, end=', ')
        input_file = [ f for f in glob.glob(os.path.join(dir,'*.inp')) if 'scfhess.inp' not in f ]
        if len(input_file) > 1:
            raise ValueError("More than one input file found.")
        elif len(input_file) == 0:
            print("\n No .inp file found in {:}, skipping.".format(dir))
            continue
        # get output file
        output_file = input_file[0].replace('.inp', '.out')
        assert os.path.isfile(output_file), "Output file {:} does not exist.".format(output_file)
        #print("Reading output file {:}".format(output_file))
        tmp, tmp_labels = read_output(output_file)
        # Use add the documented units to the labels
        tmp, unit_labels = convert_units(tmp, tmp_labels, target=unit)
        iterations.extend(tmp)
        labels.append(tmp_labels)
    print("Done!")
    n_iterations = len(iterations)
    assert len(set([" ".join(la) for la in labels])) == 1, "Found different labels in different subfolders!"
    labels = labels[0]
    print("Found the following labels: {}".format(", ".join(labels)))

    if presentation:
        lw = 5
        s = 3
        plt.rcParams.update({'font.size': 22})
        plt.rcParams.update({'legend.fontsize': 22})
    else:
        lw = 3
        s = 0

    fn = os.path.join(path, filename)
    make_go_plot(np.arange(1, n_iterations + 1), iterations, unit_labels, filelabels=labels, filename=fn, lw=lw, s=s, show=show)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plot ORCA Geometry Optimization Convergence')
    parser.add_argument('--file', help='Plot Filename Beginning, will be appended with _<label>.png', default='convergence')
    parser.add_argument('--presentation', help='Presentation Mode (i.e. thicker lines)', action='store_true')
    parser.add_argument("--unit", "-u", help="Unit System to use for plotting, choose between 'au' or 'si'", default='si')
    args = parser.parse_args()
    plot_orca_go(args.file, presentation=args.presentation, unit=args.unit)
