
import numpy as np
from .utils import arg_parsing
import .shape_classes as shp

''' Flowthrough
    # 1. Parse inputs
    # 2. Check connections
    # 3. Calculate areas
    # 4. Generate/Manipulate PDBs
    # 5. Merge
'''

def attempt_connection(shape1,shape2,mode='adjacent'):
    ''' Compares two shapes, evaluates compatibility for connection based on
        the following criteria:
            connecting_dir: are acyl chain connecting directions lined up
                            parallel or orthoganally
            symmetry:       is there a discrepency between symmetries? If so
                            some shapes can merge with connector, others can't
            dimensions:     For radially symmetric shapes, radius is important
                            for laterally symmetric shapes, compare dimensions
        Produces the following outputs:
            possible_merge: True or false, will be true if any possible
                            connection can be made
            rotation:       Is rotation necessary? [x y z] angle
            connector:      Is connecting segment necessary? If so, lateral
                            or radial? ('no','lateral',radial')
        Note that if rotations are necessary, suggested angles will be
        calculated for the SECOND object'''






# main function
def main():
    pass




# run main function
if __name__ = 'main':
    main()
