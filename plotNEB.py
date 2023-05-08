#!/usr/bin/env python3
#
# Script to plot ORCA NEB calculation results
# by Patrick Melix
# 2023/04
#
# You can import the module and then call .main() or use it as a script
# Needs grep and tail.
import argparse, os, subprocess, sys, glob
import numpy as np
from matplotlib.ticker import MaxNLocator
from matplotlib import pyplot as plt
from ase.units import create_units


def plot(reactionCoord, reactionCoordImageAxis, energies, energySpline, filename, forces=None, lw=3, s=0, highlight=None, dispersion=None, unit_label='kJ/mol'):
    msbig = 9
    ax = plt.figure().gca()
    ax.xaxis.set_major_locator(MaxNLocator(integer=True))
    plt.xlabel('Reaction Coordinate')# [Å]
    ax.set_xticklabels([]) #no numbers on x
    plt.ylabel(r'$\Delta E$ [{}]'.format(unit_label))
    #plt.ylim([-10**exp,10**exp])
    #plt.yscale('symlog')
    #plt.gca().yaxis.grid(True)
    plt.plot(reactionCoord, energySpline, color='black', ls=':', label='Cubic Spline', lw=lw)
    plt.scatter(reactionCoordImageAxis, energies, marker='P', color='red', s=(msbig+s)**2, label='NEB Energy')
    dScale = 0.02
    maxX = max(reactionCoordImageAxis)
    delta = dScale*maxX
    yRange = max(energySpline) - min(energySpline)
    for i,x in enumerate(reactionCoordImageAxis):
        if forces:
            tangentX = [x-delta, x+delta]
            tangentY = [energies[i]+(delta*forces[i]),energies[i]-(delta*forces[i])] #invert sign of forces from neb output
            if i == 0:
                label = 'NEB Force'
            else:
                label = None
            # limit y range of tangent
            tangentYRange = max(tangentY) - min(tangentY)
            factor = delta * 0.1
            n = 1
            while tangentYRange > 0.1 * yRange:
                delta -= factor
                tangentX = [x-delta, x+delta]
                tangentY = [energies[i]+(delta*forces[i]),energies[i]-(delta*forces[i])] #invert sign of forces from neb output
                tangentYRange = max(tangentY) - min(tangentY)
                n += 1
                if n >= 11:
                    raise ValueError("Tangent Problem")
            plt.plot(tangentX, tangentY, color='green', ls='-', lw=lw, label=label)

    if highlight is not None:
        plt.scatter(reactionCoordImageAxis[highlight], energies[highlight], marker='o', s=(msbig+s+30)**2, facecolors='none', edgecolors='orange', lw=lw+2, clip_on=False)

    if dispersion is not None:
        plt.scatter(reactionCoordImageAxis[1:], dispersion[1:], color='brown', marker='o', label='Dispersion', s=(msbig+s)**2)
    #plt.xticks(x, printDirs[:], rotation=90)
    plt.legend()
    plt.tight_layout()
    plt.savefig(filename)
    #plt.show()
    plt.close()


def main(filename='NEB.png', presentation=False, highlight=None, plot_all=False, plot_dispersion=False, unit='kJ/mol'):
    unitDict = create_units('2014')
    bohr2ang = unitDict['Bohr'] / unitDict['Angstrom']
    conv = unitDict['Hartree'] #ORCA uses Hartree
    if '/' in unit:
        tmp = unit.split('/')
        conv /= unitDict[tmp[0]]
        for u in tmp[1:]:
            conv *= unitDict[u]
    else:
        conv /= unitDict[unit]
    print("Conversion factor from Hartree to {}: {:}".format(unit, conv))
    print("Conversion factor from Bohr to Å: {:}".format(bohr2ang))

    interp_file = [ f for f in glob.glob('*.interp') if not ".final." in f ]
    if len(interp_file) > 1:
        raise ValueError("More than one interp file found.")
    elif len(interp_file) == 0:
        raise ValueError("No interp file found.")

    iterations = []
    n_iterations = 0
    n_images = None
    with open(interp_file[0]) as f:
        line = f.readline()
        while line:
            if 'Iteration' in line:
                iterations.append({'interpolation': {'x': [], 'dist': [], 'energy': []}, 'images': {'x': [], 'dist': [], 'energy': []}})
                n_iterations += 1
                line = f.readline() #skip header
                line = f.readline()
                while line:
                    tmp = line.split()
                    iterations[-1]['images']['x'].append(float(tmp[0]))
                    iterations[-1]['images']['dist'].append(float(tmp[1]))
                    iterations[-1]['images']['energy'].append(float(tmp[2]))
                    line = f.readline().strip()
                    if line == '':
                        if n_images is None:
                            n_images = len(iterations[-1]['images']['x'])
                        else:
                            assert len(iterations[-1]['images']['x']) == n_images, "Number of images in iteration {:} does not match previous iterations.".format(n_iterations)
                        break
            elif 'Interp.' in line:
                line = f.readline()
                while line:
                    if line == '\n':
                        break
                    else:
                        tmp = line.split()
                        iterations[-1]['interpolation']['x'].append(float(tmp[0]))
                        iterations[-1]['interpolation']['dist'].append(float(tmp[1]))
                        iterations[-1]['interpolation']['energy'].append(float(tmp[2]))
                    line = f.readline().strip()
            line = f.readline()
    print("Found {:} iterations with {} images.".format(n_iterations, n_images))


    spline = np.array([iterations[-1]['interpolation']['dist'], iterations[-1]['interpolation']['energy']])
    nebData = np.array([iterations[-1]['images']['dist'], iterations[-1]['images']['energy']])
    print("Energies loaded.")

    reactionCoord = spline[0][:] * bohr2ang
    reactionCoordImageAxis = nebData[0][:] * bohr2ang
    energies = nebData[1][:] * conv #convert to kJ/mol
    energySpline = spline[1][:] * conv #convert to kJ/mol
    forces = None
    #forces = # needs to be loaded from *.log

    if presentation:
        lw = 5
        s = 3
        plt.rcParams.update({'font.size': 22})
        plt.rcParams.update({'legend.fontsize': 22})
    else:
        lw = 3
        s = 0

    dispersion = None
    if plot_dispersion:
        raise RuntimeError("Dispersion plotting not yet implemented.")
        # print("Collecting dispersion energies from image outputs.")
        # dispersion = []
        # for i in range(nImages):
        #     path = "{:02d}".format(i)
        #     assert os.path.isdir(path), "Could not find dir {}".format(path)
        #     outcar = os.path.join(path,'OUTCAR')
        #     assert os.path.isfile(outcar), "Could not find file {}".format(outcar)
        #     child = subprocess.Popen(["grep 'Edisp (eV)' {:} | tail -1".format(outcar)], stdout=subprocess.PIPE, shell=True)
        #     dispE = float(child.communicate()[0].strip().split()[-1])
        #     dispersion.append(dispE)
        # dispersion = np.array(dispersion)
        # dispersion -= dispersion[0]
        # dispersion /= conv

    plot(reactionCoord, reactionCoordImageAxis, energies, energySpline, filename, forces, lw=lw, s=s, highlight=highlight, dispersion=dispersion, unit_label=unit)

    if plot_all:
        #plot the main image and then one with every point highlighted
        filename = filename.split('.')
        filename[-2] += "-{:02d}"
        filename = ".".join(filename)
        for i in range(n_images):
            plot(reactionCoord, reactionCoordImageAxis, energies, energySpline, filename.format(i), forces, lw=lw, s=s, highlight=i, dispersion=dispersion, unit_label=unit)


if __name__ == "__main__":
    exec(open("/home/patrickm/git/Python4ChemistryTools/mpl-settings.py").read())
    parser = argparse.ArgumentParser(description='Plot orca NEB results')
    parser.add_argument('--file', help='Plot Filename', default='NEB.png')
    parser.add_argument('--presentation', help='Presentation Mode (i.e. thicker lines)', action='store_true')
    parser.add_argument('--highlight', help='Circle Point N', type=int, default=None)
    parser.add_argument('--plotall', help='Create main plot and each highlighted plot.', action='store_true')
    parser.add_argument('--plotdispersion', help='Include dispersion contributions in plot.', action='store_true')
    parser.add_argument('--unit', help='Set the unit used to plot, must be ase compatible.', default='kJ/mol')
    args = parser.parse_args()
    main(args.file, args.presentation, args.highlight, args.plotall, args.plotdispersion, args.unit)
