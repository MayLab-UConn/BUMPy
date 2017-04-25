'''Module (or class?) of cookie cutter commands for generating shapes

   Things that should be done here:
    -organizational flow
    -flipping bilayer when relevant
    -slicing bilayers appropriately/monolayers I guess
'''
import numpy as np
# -----------------------------------------------------------------------------
class simple_shape(object):
    '''Parent class for all shapes, containing parameters each one will have'''
    def __init__(self,args):
        # things to initialize
        self.bilayer_thickness     = args.thickness
        self.area_matching_method  = args.area_matching_method
        self.outer_leaflet         = args.outer_leaflet
        self.upper_apl             = args.upper_area_per_lipid
        self.lower_apl             = args.lower_area_per_lipid
        self.template_bilayer      = molecules.Molecules(file=args.bilayer)
        self.curved_bilayer        = None
        # things to do
        if args.outer_leaflet == 'bot':   # flip upside down if necessary
            self.bilayer.coords = rb.center_coordinates_3D(self.bilayer_coords)
            self.template_bilayer.invert_z_axis()
        self.autodetect_params()
    def autodetect_params(self):
        '''Runs through optional parameters, calculates automatically from
           input bilayer if set to None
        '''
        if self.bilayer_thickness is None:
            pass
        if self.upper_area_per_lipid is None:
            pass
        if self.lower_area_per_lipid is None:
            pass
# -----------------------------------------------------------------------------

class semisphere(simple_shape):
    def __init__(self,args):
        super(simple_shape,self).__init__(args)
        self.r_sphere = args.r_sphere
    def check_params(self):
        '''Make sure correct parameters are present, warn if extraneous inputs
           are present'''
        pass
    def gen_shape(self):
        '''Main assembly function, returns a Molecules object'''
        # radius is quarter circle
        top_slice_radius = pi * (self.r_sphere + (self.thickness/2) / 2)
        bot_slice_radius = pi * (self.r_sphere - (self.thickness/2) / 2)
        # origin of circular slice
        slice_origin    = self.template_bilayer.gen_random_slice_point(
                        top_slice_radius)
        # determining slicing indices
        in_outer_circle = self.template_bilayer.circular_slice(
                        top_slice_radius,slice_origin)
        in_inner_circle = self.template_bilayer.circular_slice(
                        bot_slice_radius, slice_origin)
        top_leaflet_ind = self.template_bilayer.leaflets[ equal to 1]# write in
        bot_leaflet_ind = self.template_bilayer.leaflets[ equal to 0]
        # make slices
        top_leaflet = self.template_bilayer.slice_pdb(
                    np.intersect(in_outer_circle,top_leaflet))
        bot_leaflet = self.template_bilayer.slice_pdb(
                    np.intersect(in_inner_circle,bot_leaflet))
       # scale bottom leaflet
        bot_leaflet.coords = nrb.scale_coordinates_radially(
                           self.template_bilayer.coords,
                           (top_slice_radius / bot_slice_radius))
        # merge leaflets
        curved_bilayer = top_leaflet.append_pdb(bot_leaflet)
        curved_bilayer.coords = nrb.spherical_transform(curved_bilayer.coords,
                                                        self.r_sphere)
        return curved_bilayer

class sphere(simple_shape):
    def __init__(self,args):
        super(simple_shape,self).__init__(args)


class semi_cylinder(simple_shape):
    def __init__(self,args):
        super(simple_shape,self).__init__(args)
        self.radial_symmetry = True
        self.lateral_symmetry = True

class cylinder(simple_shape):
    def __init__(self,args):
        super(simple_shape,self).__init__(args)
        self.radial_symmetry = True
        self.lateral_symmetry = False
class torus(simple_shape):
    def __init__(self,args):
        super(simple_shape,self).__init__(args)
        self.radial_symmetry = True
        self.lateral_symmetry = False
class spherocylinder(simple_shape):
    def __init__(self,args):
        super(simple_shape,self).__init__(args)
        self.radial_symmetry = True
        self.lateral_symmetry = True
class cylinder_pierced_sphere(simple_shape):
    def __init__(self,args):
        super(simple_shape,self).__init__(args)
        self.radial_symmetry = True
        self.lateral_symmetry = False
class semicylinder_bilayer(simple_shape):
    def __init__(self,args):
        super(simple_shape,self).__init__(args)
        self.radial_symmetry = True
        self.lateral_symmetry = False
class semisphere_bilayer(simple_shape):
    def __init__(self,args):
        super(simple_shape,self).__init__(args)
        self.radial_symmetry = True
        self.lateral_symmetry = False
class cylinder_bilayer(simple_shape):
    def __init__(self,args):
        super(simple_shape,self).__init__(args)
        self.radial_symmetry = True
        self.lateral_symmetry = False
