#			memcurv 

This program is  used to generate coordinate files for use in molecular dynamics (MD) simulations of curved lipid bilayers
Please read and cite 
REFERENCE

### INSTALLATION
No installation is necessary as long as you have a python (v3) interpreter

### DEPENDENCIES
Just numpy!

### USAGE
Memcurv is designed to be used at the command line, with something like the command "python memcurv [ options ]" 

For a description of command line options, use the -h option at the command line
For a list of supported curved shapes, use the --list option at the command line

### DEVELOPMENT
* Memcurv is maintained on github (link: https://github.com/scal444/memcurv)
* I'll happily accept code contributions to improve the tool or add shapes to the repository 
* If there's a shape you want to simulate, I can likely make a template for it very quickly. Email me with a description of what you want built and I'll add it to the repository.
### Bugs
To report bugs, email Kevin Boyd at kevin.boyd@uconn.edu

### Tips

#### Input bilayer
* Templates used as starting structures for memcurv need to contain WHOLE lipid molecules (ie - not broken across periodic boundaries)
* Template bilayers do NOT need to be large enough to wrap the flat bilayer into the desired shape; memcurv will multiply the input laterally to create a flat bilayer of sufficient size. This does mean that the original box dimensions need to be correctly set 

#### Shapes
* All shapes are contained in the shapes class. The requirements for adding a shape to the repository are described in the shapes documentation
* Descriptions of supported shapes can be found in the shapes.pdf document. Note that for some shapes the terms "inner" and "outer" leaflets are ambiguous, so for those shapes a convention of "inner" and "outer" was chosen and is listed in the document. This convention does not change any of the pivotal-plane calculations. 



