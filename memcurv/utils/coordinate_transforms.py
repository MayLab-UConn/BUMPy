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
# ------------------------------------------------------------------------------
# centering coordinates
# ------------------------------------------------------------------------------
def center_coordinates_3D(coords,
                       x_center_type='mean',
                       y_center_type='mean',
                       z_center_type='range'):
    ''' Centers 3D cartesian coordinate system about the origin. Each dimension
        can be centered based on two different, criteria:
            mean: subtracts average from coordinates
            range:subtracts midpoint between maximum and minimum coordinate.
        "range" should be used for centering z axis, as asymmetric bilayers can
        tilt the average
    '''
        def center_1D(coords1D,type):
            if   type == 'mean':
                return coords1D - np.mean(coords1D)
            elif type == 'range':
                return coords1D - (np.max(coords1D) - np.min(coords1D))/2
        coords[:,0] = center_1D(coords[:,0],x_center_type)
        coords[:,1] = center_1D(coords[:,1],y_center_type)
        coords[:,2] = center_1D(coords[:,2],z_center_type)
        return coords

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
# Main curving transformations
# ------------------------------------------------------------------------------

def cylindrical_transform(xyz_coords,r,outer_leaflet='top'):
    ''' Transforms a rectangular segment of coordinates into a cylinder with
        given radius r. Completeness of cylinder will depend on length of x
        dimension, y dimension is long axis of cylinder.

        Direction of curvature by default is "top" : --(^)--, "bot" : --(_)--.
        While once made it is trivial to rotate the shapes 180 degrees to match
        the direction of curvature, this step determines which leaflet will be
        "inner" and which will be "outer". For ease of future rotations,
        all transformations here will be outer=top, with an upward protrusion
        of the cylinder. If outer_leaflet is set to 'bot', z dimension will be
        inverted about z axis prior to transformation.
    '''
    xyz_coords = center_coordinates_3D(xyz_coords)
    # determine if curving up or down
    # note that curving will be from center of bilayer.
    if outer_leaflet == 'bot':
        xyz_coords[:,2] = (-1) * xyz_coords[:,2]
    radii = r + xyz_coords[:,2]
    offset_angle = np.pi / 2
    # calculate arc lengths, this is independent of z
    arc_length = xyz_coords[:,0]
    arc_length_angle = arc_length  /  r  + offset_angle
    # cylindrical transform
    x_transform = radii * np.cos(arc_length_angle)
    z_transform = radii * np.sin(arc_length_angle)
    return np.stack((x_transform,xyz_coords[:,1],z_transform),axis=1)

def spherical_transform(xyz_coords,r,outer_leaflet='top'):
    ''' Transforms a circular segment of coordinates into a sphere with
        given radius r. Completeness of sphere will depend on radius of circular
        bilayer patch.
        Direction of curvature by default is "top" : --(^)--, "bot" : --(_)--.
        While once made it is trivial to rotate the shapes 180 degrees to match
        the direction of curvature, this step determines which leaflet will be
        "inner" and which will be "outer." For ease of future rotations,
        all transformations here will be outer=top, with an upward protrusion
        of the sphere. If outer_leaflet is set to 'bot', z dimension will be
        inverted about z axis prior to transformation.
    '''
    xyz_coords = center_coordinates_3D(xyz_coords)
    (theta,rho,z) = cart2pol(xyz_coords)

    if outer_leaflet == 'bot':
        xyz_coords[:,2] = (-1) * xyz_coords[:,2]
    radii = r + xyz_coords[:,2]
    offset_angle = np.pi / 2
    arc_length_angle = rho / r
    rho_transform = radii * np.cos(arc_length_angle)
    z_transform   = radii * np.sin(arc_length_angle)
    return pol2cart(theta,rho_transform,z_transform)

def junction_transform(xyz_coords,r_cylinder,r_junction,outer_leaflet='top'):
    ''' Tranforms a section of bilayer into a hollow junction connecting
        orthogonal bilayer segments. Bilayer should be a hollowed disk.

        Directionality is tricky here, working with two radii. We will use the
        convention of "inner" as the leaflet connecting to the inner leaflet of
        the cylindrical region of the junction. So, if outer_leaflet is set to
        top, the BOTTOM leaflet will connect to inner cylindrical species. If
        set to "bot", will invert z axis (flipping bilayer) prior to transform.
    '''
    xyz_coords = center_coordinates_3D(xyz_coords)
    (theta,rho,z) = cart2pol(xyz_coords)

    if outer_leaflet == 'bot':
        xyz_coords[:,2] = (-1) * xyz_coords[:,2]

    r_max = r_cylinder + r_junction
    # arc length starts not from 0, but first flat radius corresponding
    # to top of cylinder. If properly set up this is r_max - quarter_circle
    arclength = rho - (r_max - np.pi * r_junction / 2)
    arclength_angle = arclength / r_junction
    # rho tacks on r_cylinder, z_transform is ok with just trig transform
    rho_transform = r_junction * np.sin(arclength_angle) + r_cylinder
    z_transform   = r_junction * np.cos(arclength_angle)
    return pol2cart(theta,rho_transform,z_transform)
