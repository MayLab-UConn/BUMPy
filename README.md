memcurv (temp name)

This program will be used to generate coordinate files for use in molecular
dynamics (MD) simulations of curved lipid bilayers

--------------------------------------------------------------------------------
DEPENDENCIES
--------------------------------------------------------------------------------

numpy
argparse








--------------------------------------------------------------------------------
INSTALLATION
--------------------------------------------------------------------------------

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
1. :Argument Parsing::
Initial release of this program will be command line only. Will use argparse
module (built in, I think, don't have to list as a dependency) to parse commands
2. :Molecules::
Main class for storing bilayer properties, including coordinates, atom numbers,
residue numbers, and relevant strings for PDB output. Will have the following
functionalities
  -reading PDBs
  -writing PDBs
  -Determining (and storing) resids in top vs bottom leaflets
  -Determining (and storing) system dimensions
  -renumbering atoms and resids
  -slicing self to create smaller version of same class (ie in prep for
   transformations)
  -appending multiple instances of self (splicing partial shapes back together)
3. :Rigid body transformations::
Library of functions (either as a module or class, not sure yet) which can
manipulate coordinates in the following ways:
  -Euclidian 3D translation
    -ie, centering around origin, moving components prior to splicing
  -Euclidian 3D rotation  
4. :Non-rigid body transformations::
Library of functions to perform transformations on coordinates that involve
dimension scaling or curvature
  -coordinate scaling
  -spherical transforms
  -cylindrical transforms
  -junction transforms
  -Maybe more? -not sure if any other type is necessary right now.
5. :Shape Assemblies::
Library of commands to generate and assemble components of a desired shape, taking input from command line arguments then passing along commands to Molecules and the transformation libraries. See flow-through for an example.






--------------------------------------------------------------------------------
DEVELOPMENT AND FUTURE FEATURES
--------------------------------------------------------------------------------

--------------------------------------------------------------------------------
Bugs
--------------------------------------------------------------------------------
