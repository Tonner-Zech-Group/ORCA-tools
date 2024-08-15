import argparse
import os
import glob
import numpy as np
from matplotlib.ticker import MaxNLocator
from matplotlib import pyplot as plt

def make_go_plot(xAxis, iterations, labels, filename, lw=3, s=0, show=False):
    msbig = 9
    for i_label, label in enumerate(labels):
        data = np.array([ d[i_label][0] for d in iterations ])
        data_converged = np.array([ d[i_label][-1] for d in iterations ])
        converged_indices = np.where(data_converged is True)[0]
        not_converged_indices = np.where(data_converged is False)[0]
        ax = plt.figure().gca()
        ax.xaxis.set_major_locator(MaxNLocator(integer=True))
        plt.xlabel('Step Number')
        plt.ylabel(label)
        plt.gca().yaxis.grid(True)
        plt.plot(xAxis, data, color='black', ls=':', lw=lw)
        plt.scatter(xAxis[not_converged_indices], data[not_converged_indices], marker='x', color='red', s=(msbig+s)**2)
        plt.scatter(xAxis[converged_indices], data[converged_indices], marker='x', color='green', s=(msbig+s)**2)
        #plt.legend()
        plt.tight_layout()
        plt.savefig("{}_{}.png".format(filename,label.replace(' ', '_')))
        if show:
            plt.show()
        plt.close()


def read_output(output_file):
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
    return iterations, labels

def plot_orca_go(filename='convergence.png', presentation=False, path='.', show=False):
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
    make_go_plot(np.arange(1, n_iterations + 1), iterations, labels, filename=fn, lw=lw, s=s, show=show)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Plot ORCA Geometry Optimization Convergence')
    parser.add_argument('--file', help='Plot Filename Beginning, will be appended with _<label>.png', default='convergence')
    parser.add_argument('--presentation', help='Presentation Mode (i.e. thicker lines)', action='store_true')
    args = parser.parse_args()
    plot_orca_go(args.file, args.presentation)
