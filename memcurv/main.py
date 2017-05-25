''' Main script for memcurv project'''

from argparse import ArgumentParser
import shapes2
from molecules import Molecules
import nonrigid_coordinate_transformations as nrb
import rigid_body_transforms as rb


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
    required_inputs.add_argument('-s','--shape')
    required_inputs.add_argument('-f','--bilayer')
    # geometric parameters, which are used depends on shape input
    geometric_arguments.add_argument('--r_sphere',type=float)
    geometric_arguments.add_argument('--r_cylinder',type=float)
    geometric_arguments.add_argument('--r_junction',type=float)
    geometric_arguments.add_argument('--cylinder_length',type=float)
    geometric_arguments.add_argument('--completeness',type=float,default=1.0)
    geometric_arguments.add_argument('--flat_area',default='match_average')
    # optional arguments
    optional_arguments.add_argument('-h','--help', action='help',
                             help='show this help message and exit')
    optional_arguments.add_argument('--outer_leaflet',default='top')
    optional_arguments.add_argument('--thickness',type=float)
    optional_arguments.add_argument('--area_matching_method',default='scaling')
    optional_arguments.add_argument('-uapl','--upper_area_per_lipid')
    optional_arguments.add_argument('-lapl','--lower_area_per_lipid')
    # output files
    output_arguments.add_argument('-o') # mandatory
    output_arguments.add_argument('-p') # optional
    output_arguments.add_argument('-n') # optional
    return parser.parse_args()

def display_parameters(cl_args):
    '''Displays selected parameters upon command line execution'''
    pass # develop this

def check_argument_sanity(cl_args):
    '''Checks arguments to make sure necessary arguments are present and with
       correct file types.

       Matching geometric parameters to shapes will occur in each shape class,
       as different shapes require different parameters
    '''
    shape_options = ['sphere','cylinder','torus','spherocylinder',
                    'cylinder-pierced-sphere','semisphere-bilayer',
                    'semicylinder-bilayer','cylinder-bilayer']
    # list which parameters are applicable to shapes
    if cl_args.shape is None:
        print('Error: You must specify a shape using -s')
        return False
    if cl_args.bilayer[-4:] != '.pdb':
        print('Error: input file from -f does not have a .pdb extension\n')
        return False
    return True


def duplicate_flat_bilayer(bilayer_obj):
    '''Duplicates a flat bilayer 3 times, appends to make a bilayer 4 times
       the original size with same dimension ratios
    '''
    dimensions = bilayer_obj.get_current_dims()
    x_plus  = bilayer_obj.slice_pdb(np.arange(0,len(bilayer_obj.atomno)))
    y_plus  = bilayer_obj.slice_pdb(np.arange(0,len(bilayer_obj.atomno)))
    xy_plus = bilayer_obj.slice_pdb(np.arange(0,len(bilayer_obj.atomno)))
    x_plus.coords = x_plus.coords + [dimensions[0],0,0]
    y_plus.coords = y_plus.coords + [0,dimensions[1],0]
    xy_plus.coords = xy_plus_coords + [dimensions[0],dimensions[1],0]
    return bilayer_obj.append_pdb(x_plus).append_pdb(y_plus).append_pdb(xy_plus)

# -----------------------------------------------------------------------------
# main commands, turn into function main() at end
# -----------------------------------------------------------------------------
print('Parsing command line arguments\n')
args = parse_command_lines()
display_parameters(args)  # show user what they selected
#if not check_argument_sanity(args): # exit early if there's an issue#
#    exit
# determine which shape_assembly to call
#functions = {
#    'semisphere'             : shapes.semisphere,
#    'sphere'                 : shapes.sphere,
#    'cylinder'               : shapes.cylinder,
#    'torus'                  : shapes.torus,
#    'spherocylinder'         : shapes.spherocylinder,
#    'cylinder_pierced_sphere': shapes.cylinder_pierced_sphere,
#    'semisphere_bilayer'     : shapes.semisphere_bilayer,
#    'semicylinder_bilayer'   : shapes.semicylinder_bilayer,
#    'cylinder_bilayer'       : shapes.cylinder_bilayer
 #}
template_bilayer = Molecules(infile=args.bilayer)
shape_tobuild = shapes2.shapes.gen_shape(args,template_bilayer)
shape_tobuild.write_pdb(args.o)
shape_tobuild.write_topology(args.p)
shape_tobuild.write_index(args.n)
