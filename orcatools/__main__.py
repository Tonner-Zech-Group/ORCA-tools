
import argparse
parser = argparse.ArgumentParser(description='Plot ORCA Outputs')
subparser = parser.add_subparsers(help="Which Type of Operation to perform", dest='subcommand')
subparser.required = True

parser.add_argument('--path', type=str, help='Path to the ORCA calculation', default='.')
parser.add_argument('--presentation', help='Presentation Mode (i.e. thicker lines)', action='store_true')


neb = subparser.add_parser('neb', help='Plot NEB')
neb.add_argument('--highlight', help='Circle Point N', type=int, default=None)
neb.add_argument('--plotall', help='Create main plot and each highlighted plot.', action='store_true')
neb.add_argument('--plotdispersion', help='Include dispersion contributions in plot.', action='store_true')
neb.add_argument('--unit', help='Set the unit used to plot, must be ase compatible.', default='kJ/mol')
neb.add_argument('--file', help='PNG Filename', default='NEB')

go = subparser.add_parser('go-convergence', help='Plot Geometry Optimization Convergence')
go.add_argument('--file', help='Plot Filename Beginning, will be appended with _<label>.png', default='convergence')

interpolate = subparser.add_parser('interpolate', help='Interpolation between structures or reinterpolation of trajectories')
interpolate.add_argument("-n", "--images", type=int, help="Number of images to genereate (excluding start and end)", required=True)
interpolate.add_argument("-o", "--output", type=str, help="Output file, format can be xyz or allxyz", required=True)
interpolate.add_argument("-s", "--structures", type=str, nargs='+', help="Structurs to be included in the interpolation ordered from beginning to end", required=False, default=None)
interpolate.add_argument("-t", "--trajectory", type=str, help="Trajectory to be re-interpolated, formated as xyz or allxyz", required=False, default=None)
interpolate.add_argument("-d", "--distance", action="store_true", help="Distribute images evenly using RMSD distance. This is recommended for NEB trajectories but will break NEB-CI trajectories", required=False, default=False)
interpolate.add_argument("-v", "--verbose", action='store_true', help="Verbose output", required=False)

args = parser.parse_args()

if args.subcommand == 'neb':
    from orcatools.neb import plot_orca_neb
    plot_orca_neb(args.file, args.presentation, args.highlight, args.plotall, args.plotdispersion, args.unit, args.path)
elif args.subcommand == 'go-convergence':
    from orcatools.go import plot_orca_go
    plot_orca_go(args.file, args.presentation, args.path)
elif args.subcommand == 'interpolate':
    from orcatools.reinterpolate import reinterpolate
    reinterpolate(images=args.images, output=args.output, structures=args.structures, trajectory=args.trajectory,
                  distance=args.distance, verbose=args.verbose)