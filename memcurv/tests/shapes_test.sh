#!/bin/bash

# just runs through basic tests of every shape in repository to make sure nothing is broken
outfolder="shape_test_pdbs"


while read line; do
  shape=`echo $line | cut -d ' ' -f 1`
  gargs=`echo $line | cut -d ' ' -f 2-`
  echo "testing $shape, using arguments $gargs"
  python3 ../main.py -f ../small_flat_bilayer.pdb -z 10.5 -o $outfolder/test_$shape.pdb -s $shape -g $gargs
  if [ $? -gt 0 ]; then
    echo "ERROR OMG"
    break
  fi
done < shapes.txt

exit
