memcurv (temp name)

This program will be used to generate coordinate files for use in molecular
dynamics (MD) simulations of curved lipid bilayers

--------------------------------------------------------------------------------
DEPENDENCIES
--------------------------------------------------------------------------------
numpy
--------------------------------------------------------------------------------
INSTALLATION
--------------------------------------------------------------------------------
I think I can do this with numpy as the only dependency, merge all files into
one .py file, so installation shouldn't even be necessary.
--------------------------------------------------------------------------------
USAGE
--------------------------------------------------------------------------------
PROJNAME can be used in two ways: First, one can build from a number of
predetermined shapes, with user-inputted geometric parameters. This is the
simplest way to do things.

Second, one could build a custom shape by merging various shapes in a grid

--------------------------------------------------------------------------------
CODE SETUP
--------------------------------------------------------------------------------
This being my first real coding project, I'm aware that my setup may not be
ideal (or sensical in any way), am willing to take any suggestions.

Besides the main() bit of code, there are 5 primary sections of code, which will each be developed in more detail below
1. Argument Parsing
2. Molecules
3. Rigid body transformations
4. Non-rigid body transformations
5. Shape assemblies

DETAILS:
1. **Argument Parsing** - Initial release of this program will be command line only. Will use argparse module (built in, I think, don't have to list as a dependency) to parse commands
2. **Molecules** - Main class for storing bilayer properties, including coordinates, atom numbers, residue numbers, and relevant strings for PDB output. Will have the following functionalities:
  -reading PDBs
  -writing PDBs
  -writing simple .top files
  -Determining (and storing) resids in top vs bottom leaflets
  -Determining (and storing) system dimensions
  -renumbering atoms and resids
  -slicing self to create smaller version of same class (ie in prep for
   transformations)
  -appending multiple instances of self (splicing partial shapes back together)
3. **Rigid body transformations** - Library of functions (either as a module or class, not sure yet) which can manipulate coordinates in the following ways:
  -Euclidian 3D translation
    -ie, centering around origin, moving components prior to splicing
  -Euclidian 3D rotation  
4. **Non-rigid body transformations** - Library of functions to perform transformations on coordinates that involve dimension scaling or curvature
  -coordinate scaling
  -spherical transforms
  -cylindrical transforms
  -junction transforms
  -Maybe more? -not sure if any other type is necessary right now.
5. **Shape Assemblies** - Library of commands to generate and assemble components of a desired shape, taking input from command line arguments then passing along commands to Molecules and the transformation libraries. See flow-through for an example. The idea is that if someone wants to take my program and make a new shape with it, they would be able to do so by creating
a new shape assembly function, without modifying anything else besides adding the shape name to argument parsing

Example flowthrough for a spherocylinder (cylinder capped on both ends by a
semi-sphere):

1. *argparse* reads in commands, calls spherocylinder function in *shape      assemblies*
2. *Shape assemblies* creates new *Molecules* instance, populates it by reading in a pdb file of a flat bilayer (specified in *argparse*)
3. *Shape assemblies* checks which leaflet bilayer is "top", directs *Molecules*
instance to invert z axis if necessary
4. *Shape assemblies* directs *Molecules* to make two slices (each an      individual instance of the *Molecules* class):
  a. Pre-cylindrical slice, rectangular, size is (2 * pi * radius) by length,
     with these dimensions known to *Shape assemblies* through *argparse*
  b. Pre-semispherical slice, circular, slice radius of
     (2 * pi * radius(*from shape assembly*) / 4),
5. *Shape assemblies* manipulates coordinates of each slice using the
*non-rigid body transformations* functions. Note that we use spherical and
cylindrical transforms, without needing a separate function for semi-spheres, as the dimensions of the slice combined with the input radius will ensure the correct semi-spherical wrapping.
6. *Shape assemblies* has the semi-spherical slice is duplicated, now there are 3 instances of *Molecules* not counting the original flat bilayer.
7. *Shape assemblies* aligns the three shapes in space using the *rigid body transformations* functions to rotate and translate.
8. Shapes 2 and 3 are appended to shape 1, making one instance containing all
particles in the correct orientation
9. Resid and atom numbers are reset and organized by leaflet for clarity and
   ease of generating topology files
10. *Molecules* is directed to write out pdb and topology   

Note: The trickiest part will be steps 4 and 5. Warping a flat bilayer into a curved region changes the relative areas of the inner and outer leaflets, which
needs to be accounted for.  **MISSING EXPLANATION OF HOW WE DO THIS**

Note:
Shape assemblies can be combined for code reusage. For instance, in above example, gen_cylinder could be used for cylindrical shape, and gen_semisphere
could be used for semi sphere prior to duplication, so shape_assemblies should
build on themselves.


--------------------------------------------------------------------------------
DEVELOPMENT AND FUTURE FEATURES
--------------------------------------------------------------------------------
{list on what this has been tested}

Future Features:
-i/o extension to .gro files
-addition of proteins
-GUI implementation
-Interactive shape building (would be easiest with GUI)
--------------------------------------------------------------------------------
Bugs
--------------------------------------------------------------------------------
