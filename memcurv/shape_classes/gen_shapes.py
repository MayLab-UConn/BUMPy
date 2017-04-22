'''Module (or class?) of cookie cutter commands for generating shapes

   Things that should be done here:
    -organizational flow
    -flipping bilayer when relevant
    -slicing bilayers appropriately/monolayers I guess
'''




def gen_sphere(bilayer_pdb_obj,radius,apl=None,outer_leaflet='top'):
    '''The only geometric inputs for sphere are the radius and designation of
       outer vs inner leaflet'''
     # First, invert:
    if outer_leaflet == 'bot':
         bilayer_pdb_obj.invert_z_axis()
    #calculate area differences
    top_area = sphere.area(radius + (thickness/2)
    bot_area = sphere.area(radius - (thickness/2) # thickness needs thought
    # calculate slice_radius
    top_slice_radius = # set surface area equal to a flat circle area
    bot_slice_radius = # same thing
    # find slice point
    sl_center = bilayer_pdb_ojb.gen_circular_slicepoint(slice_radius)
    # top layer is a simple slice
    top_layer = bilayer_pdb_obj.circular_slice(sl_center,top_slice_radius)
    # bot layer needs scaling to top_layer size
    bot_layer = bilayer_pdb_obj.circular_slice(sl_center,bot_slice_radius)
    bot_layer.scale_coordinates_radially(top_slice_radius)
    # merge two together
    merged_bilayer = top_layer.append_pdb(bot_layer)
    # perform transform
    merged_bilayer.coords = nonrigid_coordinate_transforms.spherical_transform(coords)
    return merged_bilayer

def gen_cylinder(bilayer_pdb_obj,radius,length,apl=None,outer_leaflet='top'):
    pass
def gen_semicylinder_plane():
    pass
def gen_semisphere_plane():
    pass
def gen_torus():
    pass
def gen_sphere_cylinder_pbc():
    pass
def gen_spherocylinder():
