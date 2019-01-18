#			BUMPy

This program is used to generate coordinate files for use in molecular dynamics (MD) simulations of curved lipid bilayers

### CITATION
If you use our tool, please read and cite the BUMPy publication:

Kevin Boyd and Eric May. BUMPy: A Model-Independent Tool for Constructing Lipid Bilayers of Varying Curvature and Composition. 
Journal of Chemical Theory and Computation, 2018, 14(12), pp 6642-6652.

### AUTHORS
BUMPy is written and maintained by Kevin Boyd (kevin.boyd@uconn.edu), in the lab of Eric May at the University of Connecticut.

### INSTALLATION
No installation is necessary as long as you have a python (v3) interpreter

### DEPENDENCIES
Just numpy and scipy!

### USAGE
#### Command line tool
BUMPY is designed to be used at the command line, with something like the command "python bumpy.py [ options ]". All you
need is the bumpy.py file! You can also take the file out of the BUMPy directory and use it wherever you want.
* For a description of command line options, use the -h option at the command line.
* For a list of supported curved shapes, use the --list option at the command line, or see the shape repository pdf.
* Example commands with BUMPy are available in the examples directory
#### Accurate area estimation
A novel feature of BUMPy is quantitative estimation of areas in curved bilayers using the monolayer pivotal plane - a surface
within the monolayer that does not undergo area changes upon curvature deformations. On the command line, this option is
controlled with the -z flag. The pivotal plane location is composition-dependent.
#### Does having an accurate pivotal plane estimate matter?
It depends on what you're using the system for. Having an inaccurate pivotal plane estimate when building these systems
leads to area mismatch. We quantify some of the effects of such area mismatch in our publication. The summary is,
the effects of this area mismatch are quantifiable for a number of observables such as lipid splay and diffusion, but
generally quite small if your estimated pivotal plane location is within a few Angstroms of the true value. This is good
news, as most of the lipids we've calculated pivotal plane locations for only vary in location by a few angstroms. 10 A
is a good starting guess. The key is to know what properties you want to measure, and if the subtle differences due to
area mismatch are on the same order as the degree of precision you want in your measurements. Do take a look at our publication
to see what we're talking about!
#### Suggested values of pivotal planes
* If you don't care overly much about area matching, and just want to build a shape, feel free to use 10A as the -z input,
and you should have reasonable inner and outer leaflet areas
* If you're working with the Martini forcefield, we've calculated pivotal plane locations for over a dozen lipids, found
in the pivotal_planes folder. We hope to add to this repository over time
* If you want a more quantitative estimate of the pivotal plane location, you can infer a location from some flat bilayer
properties such as thickness or the lateral pressure profile, though we haven't managed to find an exact relationship between
z values and flat bilayer properties (and we've tried). See our publication for details.
* If you absolutely need rigorous area matching, and the z value of the lipid you want to simulate hasn't been calculated yet,
you can do it yourself! See our publication and cite the authors who came up with the z measurement- (Wang and Deserno, J. Chem. Phys. 16, 164109 (2015)). If you do calculate your own pivotal planes, please let us know and we'll update our repository with your reported values!

### Equilibration of BUMPy systems
By mixing and matching different building block shapes, the potential arises for clashes at shape interfaces, which can
(and typically does) lead to non-finite forces during energy minimization. To allow minimization to proceed, we use soft-core
potentials to scale down short-range nonbonded interactions. If you are using Gromacs, the following .mdp snippet (taken from the CHARMM-GUI's
suggested minimization scheme) can be used - just paste it into a typical minimization script, and minimization should work.
```
free-energy              = yes
init-lambda              = 0.01
sc-alpha                 = 4
sc-power                 = 2
sc-coul                  = yes
nstdhdl                  = 0
couple-moltype           = system
; we are changing both the vdw and the charge. In the initial state, both are on
couple-lambda0           = vdw-q
; in the final state, both are off.
couple-lambda1           = none
couple-intramol          = yes
```
Please note that the use of soft-core potentials slows down minimization by about an order of magnitude. We therefore suggest
a brief (~50 step) minimization using soft-core potentials, followed by a typical minimization without soft-core potentials. We
have found that every system we've created in BUMPy can be successfully minimized with these techniques, so please do let us know
if you come across a usage case where soft-core potentials are not sufficient!

### DEVELOPMENT
* BUMPy is maintained on github (link: https://github.com/MayLab-UConn/BUMPy)
* I'll happily accept code contributions to improve the tool or add shapes to the repository.
* Work will be done on the dev branch, and merged into master at official release points. This is so that
  one can pinpoint a release version if a science-affecting bug is found - but hopefully it won't matter!
* If there's a shape you want to simulate, I can likely make a template for it very quickly. Email me with a description of what you want built and I'll add it to the repository.
### Bugs
To report bugs, email Kevin Boyd at kevin.boyd@uconn.edu

Known bugs:
* .gro reading is unreliable. Suggest using pdb files as input
### Tips
* Units are in Angstroms - be careful with your inputs!
#### Input bilayer
* Templates used as starting structures for BUMPy need to contain WHOLE lipid molecules (ie - not broken across periodic boundaries)
* Template bilayers do NOT need to be large enough to wrap the flat bilayer into the desired shape; BUMPy will multiply the input laterally to create a flat bilayer of sufficient size. This does mean that the original box dimensions need to be correctly set
#### Shapes
* All shapes are contained in the shapes class. The requirements for adding a shape to the repository are described in the shapes documentation
* Descriptions of supported shapes can be found in the shapes.pdf document. Note that for some shapes the terms "inner" and "outer" leaflets are ambiguous, so for those shapes a convention of "inner" and "outer" was chosen and is listed in the document. This convention does not change any of the pivotal-plane calculations.
