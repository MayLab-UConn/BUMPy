''' Main script for memcurv project'''

from argparse import ArgumentParser
import shapes2
from molecules import Molecules
import nonrigid_coordinate_transformations as nrb
import rigid_body_transforms as rb
import numpy as np

__version__ = '0.5'



def parse_command_lines():
    ''' Parses command line for parameters, returns parsed arguments '''
    prog_name =  ''
    prog_description = ''
    geometry_description = 'temp_geometry'
    optional_description = 'temp_optional'
    output_description   = 'temp_output'

    parser = ArgumentParser(prog=prog_name,description=prog_description,
                            add_help=False)
    # groups
    required_inputs     = parser.add_argument_group('required inputs')
    geometric_arguments = parser.add_argument_group('geometric arguments',
                          geometry_description)
    optional_arguments  = parser.add_argument_group('optional arguments',
                          optional_description)
    output_arguments    = parser.add_argument_group('output arguments',
                          output_description)
    # mandatory input
    required_inputs.add_argument('-s','--shape',required=True)
    required_inputs.add_argument('-f','--bilayer',required=True)
    required_inputs.add_argument('-g','--geometry',nargs='*',required=True)
    required_inputs.add_argument('-z','--zo',required=True, type=float)

    # optional arguments
    optional_arguments.add_argument('-h','--help', action='help',
                             help='show this help message and exit')
    optional_arguments.add_argument('--outer_leaflet',default='top')
    optional_arguments.add_argument('-uapl','--upper_area_per_lipid')
    optional_arguments.add_argument('-lapl','--lower_area_per_lipid')
    # output files
    output_arguments.add_argument('-o') # mandatory
    output_arguments.add_argument('-p') # optional
    output_arguments.add_argument('-n') # optional
    return parser.parse_args()

def display_parameters(cl_args):
    '''Displays selected parameters upon command line execution'''
    pass


# -----------------------------------------------------------------------------
# main commands, turn into function main() at end
# -----------------------------------------------------------------------------
print('Parsing command line arguments\n')
args = parse_command_lines()
display_parameters(args)  # show user what they selected

# parse arguments
geometric_args = {garg.split(':')[0]:float(garg.split(':')[1]) for garg in args.geometry } # use a comprehension
zo = args.zo
#shape_tobuild = getattr(shapes2.shapes,args.shape) This will be correct when compiled together
shape_tobuild = getattr(shapes2,args.shape)

# adjust size of template bilayer
template_bilayer = Molecules(infile=args.bilayer)
mult_factor = np.ceil( shape_tobuild.dimension_requirements(**geometric_args)/template_bilayer.boxdims[0:2]).astype(int)
template_bilayer.duplicate_laterally(*mult_factor)

# construct the shape
shape = shape_tobuild.gen_shape(template_bilayer,zo,**geometric_args)
shape.boxdims = shape_tobuild.final_dimensions(**geometric_args)
# file output
shape.write_pdb(args.o)
if args.p:
    shape.write_topology(args.p)
if args.n:
    shape.write_index(args.n)
