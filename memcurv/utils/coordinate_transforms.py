
import numpy as np
import math

def cart2pol(xyz_coords):
    ''' Converts cartesian to polar coordinates. Acts on first 2 columns (xy)
        returns tuple of (theta,rho,z)
    '''
    x = xyz_coords[:,0]
    y = xyz_coords[:,1]
    z = xyz_coords[:,2]
    
    theta    = np.arctan2(y,x)
    rho  = np.sqrt(x**2 + y**2)

    return(theta,rho,z)


def pol2cart(theta,rho,z):
    ''' Converts polar to cartesian coordinates, outputs as nparts * xyz array
    '''
    x = rho * np.cos(theta)
    y = rho * np.sin(theta)
    xyz_coords = np.stack((x,y,z),axis=1)
    return xyz_coords

def cylindrical_transform(xyz_coords,r,curv_dir='up'):
    ''' Transforms a rectangular segment of coordinates into a cylinder with
        given radius r. Completeness of cylinder will depend on length of x
        dimension, y dimension is long axis of cylinder
        curv_center
    '''
    xyz_coords = xyz_coords - np.mean(xyz_coords,axis = 0)
    (theta,rho,z) = cart2pol(xyz_coords)
    radii = r + z
    
    # determine if curving up or down
    if curv_dir=='up':
        offset_angle = 90
    else if curv_dir=='down':
        offset_angle = 270
        
    # calculate arc lengths
    arc_angle = 360 * 
    
    