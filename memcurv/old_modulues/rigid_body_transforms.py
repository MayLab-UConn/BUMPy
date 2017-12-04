''' Rigid body transforms include rotation and translation '''

import numpy as np

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
            return coords1D - (np.max(coords1D) + np.min(coords1D))/2
    coords[:,0] = center_1D(coords[:,0],x_center_type)
    coords[:,1] = center_1D(coords[:,1],y_center_type)
    coords[:,2] = center_1D(coords[:,2],z_center_type)
    return coords

def translate_coordinates(coords,directional_vectors):
    '''simple coordinate translation, direction is [x y z]'''
    return coords + directional_vectors

def rotate_coordinates(coords,rotation_angles,com=False,unit='degrees'):
    ''' Rotational transform of coordinate system. Default rotation center is
        about the origin (com=False). If com is set to True, rotation will occur
        about center of mass of coordinates (with no weighting of particles).
        Default input is in degrees, other option is 'radians'
    '''
    if unit == 'degrees': # calculations are in radians
        rotation_angles = np.radians(rotation_angles)
    if com:
        com_coords = np.mean(coords,axis=0)
        coords = coords - com_coords # center around axis, will translate later
    # develop rotation matrices (taken from google, 3D rotations)
    r_x = np.array([[1,                       0,                             0],
                    [0, np.cos(rotation_angles[0]),-np.sin(rotation_angles[0])],
                    [0, np.sin(rotation_angles[0]),np.cos(rotation_angles[0])]])

    r_y = np.array([[np.cos(rotation_angles[1]) ,0, np.sin(rotation_angles[1])],
                    [0                       ,1,                           0  ],
                    [-np.sin(rotation_angles[1]),0,np.cos(rotation_angles[1])]])

    r_z = np.array([[np.cos(rotation_angles[2]),-np.sin(rotation_angles[2]), 0],
                    [np.sin(rotation_angles[2]),  np.cos(rotation_angles[2]),0],
                    [0,                       0,                            1]])
    rotmat = r_x.dot(r_y).dot(r_z)       # net rotation matrix
    rotcoords = np.dot(coords,rotmat)    # do rotation
    if com:
        rotcoords = rotcoords + com_coords # return to original position
    return rotcoords
