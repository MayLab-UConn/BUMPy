'''Module (or class?) of cookie cutter commands for generating shapes

   Things that should be done here:
    -organizational flow
    -flipping bilayer when relevant
    -slicing bilayers appropriately/monolayers I guess
'''
import numpy as np
import rigid_body_transforms as rb
import nonrigid_coordinate_transformations as nrb
from copy import deepcopy,copy
# -----------------------------------------------------------------------------
class simple_shape:
    '''Parent class for all shapes, containing parameters each one will have'''
    def __init__(self,args,template_bilayer):
        # things to initialize
        self.thickness             = args.thickness
        self.area_matching_method  = args.area_matching_method
        self.outer_leaflet         = args.outer_leaflet
        self.upper_apl             = args.upper_area_per_lipid
        self.lower_apl             = args.lower_area_per_lipid
        self.template_bilayer      = template_bilayer
        self.curved_bilayer        = None
        # things to do
        self.template_bilayer.coords = rb.center_coordinates_3D(
                                       self.template_bilayer.coords)
        if args.outer_leaflet == 'bot':   # flip upside down if necessary
            self.template_bilayer.invert_z_axis()
        self.autodetect_params()
        self.thickness = float(self.thickness)
    def autodetect_params(self):
        '''Runs through optional parameters, calculates automatically from
           input bilayer if set to None
        '''
        if self.thickness is None:
            pass
        if self.upper_apl is None:
            pass
        if self.lower_apl is None:
            pass
# -----------------------------------------------------------------------------

class semisphere(simple_shape):
    def __init__(self,args,template_bilayer):
        super().__init__(args,template_bilayer)
        self.r_sphere = float(args.r_sphere)
    def check_params(self):
        '''Make sure correct parameters are present, warn if extraneous inputs
           are present'''
        pass
    def gen_shape(self):
        '''Main assembly function, returns a Molecules object'''
        # radius is quarter circle
        print('Calculating slices')
        slice_radius = np.pi * self.r_sphere / 2
        top_slice_radius = np.pi * (self.r_sphere + (self.thickness/2)) / 2
        bot_slice_radius = np.pi * (self.r_sphere - (self.thickness/2)) / 2
        # origin of circular slice
        slice_origin    = self.template_bilayer.gen_random_slice_point(
                        top_slice_radius)
        # determining slicing indices
        in_outer_circle = self.template_bilayer.circular_slice(
                        slice_origin,top_slice_radius)
        in_inner_circle = self.template_bilayer.circular_slice(
                        slice_origin,bot_slice_radius)
        top_leaflet_ind = np.where(self.template_bilayer.leaflets == 1)[0]
        bot_leaflet_ind = np.where(self.template_bilayer.leaflets == 0)[0]
        # make slices

        top_leaflet = self.template_bilayer.slice_pdb(
                    np.intersect1d(in_outer_circle,top_leaflet_ind))
        bot_leaflet = self.template_bilayer.slice_pdb(
                    np.intersect1d(in_inner_circle,bot_leaflet_ind))
       # scale bottom leaflet
        top_leaflet.write_pdb('test_top.pdb')
        bot_leaflet.write_pdb('test_bot.pdb')
        print('Scaling and merging slices')
        bot_leaflet.coords = nrb.scale_coordinates_radial(bot_leaflet.coords,
                           (slice_radius / bot_slice_radius))
        top_leaflet.coords = nrb.scale_coordinates_radial(
                           top_leaflet.coords,
                           (slice_radius /  top_slice_radius))
        #import pdb; pdb.set_trace()

        # merge leaflets
        top_leaflet.append_pdb(bot_leaflet) # now has both layers
        top_leaflet.write_pdb('test_merge.pdb')
        print('Transforming coordinates')
        top_leaflet.coords = nrb.spherical_transform(top_leaflet.coords,
                                                     self.r_sphere)
        return top_leaflet

class sphere(simple_shape):
    def __init__(self,args):
        super(simple_shape,self).__init__(args,template_bilayer)
        self.r_sphere = float(args.r_sphere)
    def gen_shape(self):
        semisphere_obj = semisphere(args)
        top_half_sphere = semisphere_obj.gen_shape()
        bot_half_sphere = copy(top_half_sphere)
        bot_half_sphere.coords = rb.rotate_coordinates(bot_half_sphere.coords,
                                                       [180,0,0])

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
