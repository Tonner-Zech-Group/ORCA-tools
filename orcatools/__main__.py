
import argparse
parser = argparse.ArgumentParser(description='Plot ORCA Outputs')
subparser = parser.add_subparsers(help="Which Type of Calculation to Plot", dest='subcommand')
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

args = parser.parse_args()

if args.subcommand == 'neb':
    from orcatools.neb import plot_orca_neb
    plot_orca_neb(args.file, args.presentation, args.highlight, args.plotall, args.plotdispersion, args.unit, args.path)
elif args.subcommand == 'go-convergence':
    from orcatools.go import plot_orca_go
    plot_orca_go(args.file, args.presentation, args.path)