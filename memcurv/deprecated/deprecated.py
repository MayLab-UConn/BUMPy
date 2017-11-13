

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
