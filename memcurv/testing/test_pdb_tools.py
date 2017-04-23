''' Hopefully comprehensive testing of the functions of the pdb_tools class

    test in IPYTHON!!!
'''
%load_ext autoreload
%autoreload 2

import numpy as np
import utils.pdb_tools as pdb

# Allows for on the fly modification of code


flat_bilayer_file = '/Users/kevinboyd/Projects/tim23/Lipid_Curvature/folders_waiting/PC_80_CL1_20/charmm-gui/gromacs/step5_extended.pdb'

start_bilayer = pdb.bilayer_pdb()
start_bilayer.read_pdb(flat_bilayer_file)

#test get_current_dims

#test get_current_minmax
