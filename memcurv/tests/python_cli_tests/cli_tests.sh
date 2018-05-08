#!/bin/bash
executable=../../main.py
input=../asymm.pdb
zo=10.5
# ------------------------------------------------------------------------------------
# runs through basic tests of every shape in repository to make sure nothing is broken
# ------------------------------------------------------------------------------------
outfolder="shape_test_files"

while read line; do
  shape=`echo $line | cut -d ' ' -f 1`
  gargs=`echo $line | cut -d ' ' -f 2-`
  echo "testing $shape, using arguments $gargs"
  python3 $executable -f $input -z $zo -o $outfolder/test_$shape.pdb -s $shape -g $gargs
  if [ $? -gt 0 ]; then
    echo "ERROR OMG"
    break
  fi
done < shapes.txt


# ------------------------------------------------------------------------------------
# runs tests on dummy parameter inputs
# ------------------------------------------------------------------------------------
outfolder="dummy_test_files"
python3 $executable --gen_dummy_particles -f $input -z $zo -o $outfolder/test_basic_dummy.pdb  -s semicylinder_plane -g r_cylinder:100 r_junction:50 l_flat:100 l_cylinder:100 -p $outfolder/test_basic_dummy.top -n $outfolder/test_basic_dummy.ndx --dummy_grid_thickness 50
python3 $executable --gen_dummy_particles -f $input -z $zo -o $outfolder/test_dummy_name.pdb   -s semicylinder_plane -g r_cylinder:100 r_junction:50 l_flat:100 l_cylinder:100 -p $outfolder/test_dummy_name.top  -n $outfolder/test_dummy_name.ndx  --dummy_name TEST --dummy_grid_thickness 50
python3 $executable --gen_dummy_particles -f $input -z $zo -o $outfolder/test_dummy_grid.pdb   -s semicylinder_plane -g r_cylinder:100 r_junction:50 l_flat:100 l_cylinder:100  --dummy_grid_spacing 3 --dummy_grid_thickness 50


# ------------------------------------------------------------------------------------
# test residue exclusion command line options
# ------------------------------------------------------------------------------------
outfolder="exclusion_test_files"
python3 $executable -f $outfolder/asymm_solvated.pdb -z $zo -o $outfolder/test_exclusion_1excl.pdb -n $outfolder/test_exclusion_1excl.ndx -p $outfolder/test_exclusion_1excl.top  --ignore_resnames W        -s sphere -g r_sphere:50
python3 $executable -f $outfolder/asymm_solv_ion.pdb -z $zo -o $outfolder/test_exclusion_2excl.pdb -n $outfolder/test_exclusion_2excl.ndx -p $outfolder/test_exclusion_2excl.top  --ignore_resnames W:SOD    -s sphere -g r_sphere:50
exit
