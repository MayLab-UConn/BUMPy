			memcurv   v 1.0

This program will be used to generate coordinate files for use in molecular
dynamics (MD) simulations of curved lipid bilayers, or just visualization

Please read and cite 
REFERENCE

--------------------------------------------------------------------------------
INSTALLATION
--------------------------------------------------------------------------------
No installation is necessary as long as you have a python (v3) interpreter

--------------------------------------------------------------------------------
DEPENDENCIES
--------------------------------------------------------------------------------
just numpy!


--------------------------------------------------------------------------------
USAGE
--------------------------------------------------------------------------------
memcurv is designed to be used at the command line, with something like the
command "python memcurv [ options ]" 

For a description of command line options, use the -h option at the command line
For a list of supported curved shapes, use the --list option at the command line


--------------------------------------------------------------------------------
DEVELOPMENT
--------------------------------------------------------------------------------
memcurv is maintained on github (link: https://github.com/scal444/memcurv)
I'll happily accept code contributions to improve the tool or add shapes to 
the repository 


--------------------------------------------------------------------------------
Bugs
--------------------------------------------------------------------------------

To report bugs, email Kevin Boyd at kevin.boyd@uconn.edu


--------------------------------------------------------------------------------
tips
--------------------------------------------------------------------------------

Input bilayer
-Templates used as starting structures for memcurv need to contain WHOLE lipid
 molecules (ie - not broken across periodic boundaries)
-Template bilayers do NOT need to be large enough to wrap the flat bilayer into
 the desired shape - memcurv will multiply the input laterally to create a flat
 bilayer of sufficient size. this does mean that the original box dimensions
 need to be correctly set 

Shapes
-

