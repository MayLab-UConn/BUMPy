''' Coordinate transformations that are to act on EVERY particle in a system.
    Program should be set up so that all exclusions, deletions, etc are
    processed prior to transformation.

    With that in mind,cylindrical and spherical transformation functions can be
    used to create full spheres and cylinders, half, or any completeness in
    between with careful determination of input bilayer size and desired radius.

    In addition, the cylindrical transformation function will be used to
    generate lateral junctions. Radial junctions require a separate function,
    as they depend both on the junction radius and cylindrical radius that is
    being connected.

    These transforms in each case lead to shapes pointing "up" from the bilayer.
    If the option for "outer leaflet" in some cases is set to the bottom bilayer
    leaflet, the flat bilayer will be inverted in the z direction prior to
    transform.

    Geometrically, the original bilayers will be centered around (0,0,0). For
    partial cylinders, spheres, etc, they will center around (0,0,+r) and wrap
    down. For junctions this is (0,0,+r_junction), with a gap at the TOP of
    r_cylinder, and tapering down at z=0 to a flat bilayer.
'''


import numpy as np
import rigid_body_transforms as rb


# ------------------------------------------------------------------------------
# cartesian to polar and back
# ------------------------------------------------------------------------------

def cart2pol(cart_coords):
    ''' Converts cartesian to polar coordinates. Acts on first 2 columns (xy)
        returns tuple of (theta,rho,z)
    '''
    x = cart_coords[:,0]
    y = cart_coords[:,1]
    z = cart_coords[:,2]
    theta    = np.arctan2(y,x)
    rho  = np.sqrt(x**2 + y**2)
    return(theta,rho,z)

def pol2cart(theta,rho,z):
    ''' Converts polar to cartesian coordinates, outputs as nparts * xyz array
    '''
    x = rho * np.cos(theta)
    y = rho * np.sin(theta)
    cart_coords = np.stack((x,y,z),axis=1)
    return cart_coords

# ------------------------------------------------------------------------------
# Coordinate scaling
# ------------------------------------------------------------------------------
def scale_coordinates_radial(coords,ratio):
    '''Coordinate scaling centered on 0'''
    meanvals = np.mean(coords,axis=0)
    coords = coords - [meanvals[0],meanvals[1],0]
    (theta,rho,z) = cart2pol(rb.center_coordinates_3D(coords))
    rho = rho * ratio
    return pol2cart(theta,rho,z) + [meanvals[0],meanvals[1],0]
def scale_coordinates_rectangular(coords,ratio):
    minvals = np.min(coords,axis=0)
    coords = coords - [minvals[0], minvals[1],0] # push to 0,0,0 for minima
    coords[:,0:2] = coords[:,0:2] * ratio
    return coords +  [minvals[0],minvals[1],0]
def scale_coordinates_toroid(coords,current_range,new_range):
    '''Radial coordinate scaling from one range of spaces to another'''
    meanvals = np.mean(coords,axis=0)
    coords = coords - [meanvals[0],meanvals[1],0]
    (theta,rho,z) = cart2pol(rb.center_coordinates_3D(coords))
    curr_range_size = current_range[1] - current_range[0]
    midpoint  = (current_range[0] + current_range[1]) / 2
    new_range_size  = new_range[1] - new_range[0]
    ratio = new_range_size / curr_range_size
    rho = (rho - midpoint)*ratio  + midpoint
    return pol2cart(theta,rho,z) + [meanvals[0],meanvals[1],0]


# ------------------------------------------------------------------------------
# Main curving transformations
# ------------------------------------------------------------------------------

def cylindrical_transform(xyz_coords,r,outer_leaflet='top'):
    ''' Transforms a rectangular segment of coordinates into a cylinder with
        given radius r. Completeness of cylinder will depend on length of y
        dimension, x dimension is long axis of cylinder.

    '''
    xyz_coords = rb.center_coordinates_3D(xyz_coords)
    radii = r + xyz_coords[:,2]
    offset_angle = np.pi / 2
    # calculate arc lengths, this is independent of z
    arc_length = xyz_coords[:,1]
    arc_length_angle = arc_length  /  r
    # cylindrical transform
    y_transform = radii * np.sin(arc_length_angle)
    z_transform = radii * np.cos(arc_length_angle)
    return np.stack((xyz_coords[:,0],y_transform,z_transform),axis=1)

def spherical_transform(xyz_coords,r,outer_leaflet='top'):
    ''' Transforms a circular segment of coordinates into a sphere with
        given radius r. Completeness of sphere will depend on radius of circular
        bilayer patch.

    '''
    xyz_coords = rb.center_coordinates_3D(xyz_coords)
    (theta,rho,z) = cart2pol(xyz_coords)

    radii = r + z
    arc_length_angle = rho / r
    rho_transform = radii * np.sin(arc_length_angle)
    z_transform   = radii * np.cos(arc_length_angle)
    return pol2cart(theta,rho_transform,z_transform)
def toroidal_transform(xyz_coords,r_torus,r_tube):
    xyz_coords = rb.center_coordinates_3D(xyz_coords)
    (theta,rho,z) = cart2pol(xyz_coords)
    arc_length = rho - (r_torus - (np.pi * r_tube / 2))
    arc_length_angle = arc_length / r_tube
    radii = r_tube + z
    z_transform = radii * np.sin(arc_length_angle)
    rho_transform = r_torus + radii * np.sin(arc_length_angle - np.pi/2)

    return pol2cart(theta,rho_transform,z_transform)
def junction_transform(xyz_coords,r_cylinder,r_junction):
    ''' Tranforms a section of bilayer into a hollow junction connecting
        orthogonal bilayer segments. Bilayer should be a hollowed disk.

        Directionality is tricky here, working with two radii. We will use the
        convention of "inner" as the leaflet connecting to the inner leaflet of
        the cylindrical region of the junction. So, if outer_leaflet is set to
        top(in shapes), the BOTTOM leaflet will connect to the
        inner cylindrical species.
    '''
    xyz_coords = rb.center_coordinates_3D(xyz_coords)
    (theta,rho,z) = cart2pol(xyz_coords)

    r_max = r_cylinder + r_junction
    # arc length starts not from 0, but first flat radius corresponding
    # to top of cylinder. If properly set up this is r_max - quarter_circle
    arclength = rho - (r_max - np.pi * r_junction / 2)
    arclength_angle = arclength / r_junction
    # rho tacks on r_cylinder, z_transform is ok with just trig transform
    rho_transform = r_junction * np.sin(arclength_angle) + r_cylinder
    z_transform   = r_junction * np.cos(arclength_angle)
    return pol2cart(theta,rho_transform,z_transform)
