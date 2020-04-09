# Runs commands to generate reference data that is compared to in the unit tests
# Every file regenerated should be checked for correctness

bumpy=../../bumpy.py
# topology
python3 $bumpy -f reference_files/input/input_asymm.gro -s sphere -z 10 -o temp.gro -g r_sphere:50 -p reference_files/test_topology/topol_complex.top
python3 $bumpy -f reference_files/input/input_asymm.gro -s sphere -z 10 -o temp.gro -g r_sphere:50 -p reference_files/test_topology/topol_with_dummy.top --gen_dummy_particles --dummy_grid_thickness 50

# index
python3 $bumpy -f reference_files/input/input_asymm.gro -s sphere -z 10 -o temp.gro -g r_sphere:50 -n reference_files/test_index/index_complex.ndx
python3 $bumpy -f reference_files/input/input_asymm.gro -s sphere -z 10 -o temp.gro -g r_sphere:50 -n reference_files/test_index/index_with_dummy.ndx --gen_dummy_particles --dummy_grid_thickness 50

# basic shapes
common_args="python3 $bumpy -f reference_files/input/input_asymm.gro -z 10"
$common_args -o reference_files/test_shapes/test_buckle.pdb                  -s buckle                  -g r_buckle:50 l_buckle:60
$common_args -o reference_files/test_shapes/test_capped_cylinder.pdb         -s capped_cylinder         -g r_cylinder:50 l_cylinder:60
$common_args -o reference_files/test_shapes/test_cylinder.pdb                -s cylinder                -g r_cylinder:50 l_cylinder:60
$common_args -o reference_files/test_shapes/test_double_bilayer_cylinder.pdb -s double_bilayer_cylinder -g r_cylinder:50 l_cylinder:50 r_junction:40 l_flat:250
$common_args -o reference_files/test_shapes/test_flat_bilayer.pdb            -s flat_bilayer            -g x_dimension:50 y_dimension:40
$common_args -o reference_files/test_shapes/test_inner_quarter_torus.pdb     -s inner_quarter_torus     -g r_torus:75 r_tube:50
$common_args -o reference_files/test_shapes/test_outer_quarter_torus.pdb     -s outer_quarter_torus     -g r_torus:75 r_tube:50
$common_args -o reference_files/test_shapes/test_semicylinder_plane.pdb      -s semicylinder_plane      -g r_cylinder:50 l_cylinder:30 r_junction:40 l_flat:50
$common_args -o reference_files/test_shapes/test_semisphere_plane.pdb        -s semisphere_plane        -g r_sphere:50 r_junction:30 l_flat:200
$common_args -o reference_files/test_shapes/test_semisphere.pdb              -s semisphere              -g r_sphere:50
$common_args -o reference_files/test_shapes/test_torus.pdb                   -s  torus                  -g r_torus:75 r_tube:50
