# checks 1 FRAME for io ratio, used for validating initial shape
# because things are going wrong and it would be nice to
# have a quick summary accessible from the command line

import sys
import mdtraj as md
import numpy as np
infile = sys.argv[1]
traj = md.load(infile)
traj.center_coordinates()
coords = traj.xyz[:,traj.topology.select('name PO4'),:]

def cart2spherical(cart_coords,rho=None):
    # can calculate rho individually, or take in preset
    x = cart_coords[:,:,0]  # flat, just return
    y = cart_coords[:,:,1]
    z = cart_coords[:,:,2]
    if rho is None:
        rho = np.sqrt(x**2 + y**2 + z**2)
    theta = np.arctan2(y,x)
    phi = np.arctan2(np.sqrt(x**2 + y**2),z)
    return rho,theta,phi

def separate_by_local_radii(cartcoords_oneframe,rho,search_radius):
    ''' Attempts to assign inner or outer leaflet on the local level. Given one
        frame of coordinates (cartesian and radial), finds nearest neighbors,
        separates by midpoint between maximal and minimal radial coordinates

        cartcoords should be nparts * xyz, rho should be nparts

        to speed things up, will assign every particle within that radius
    '''
    n_lips = cartcoords_oneframe.shape[0]
    inner_indices = []; outer_indices = []
    already_assigned = np.zeros(n_lips)
    for i in range(n_lips):
        # check if has been assigned
        if not already_assigned[i]:
            dists = np.linalg.norm(cartcoords_oneframe[i,:] - cartcoords_oneframe,axis=1)
            indices = np.where(dists < search_radius)[0]

            radii = rho[indices]
            init_splitter = (np.min(radii) + np.max(radii) ) / 2
            # refine initial split radius. From max/min to mean of leaflets
            refined_min = np.mean(radii[radii < init_splitter])
            refined_max = np.mean(radii[radii > init_splitter])
            splitter = (refined_min + refined_max) / 2
            # assign all indices in range
            for j in indices:
                if not already_assigned[j]:
                    if rho[j] < splitter:
                        inner_indices.append(j)
                    else:
                        outer_indices.append(j)
                    already_assigned[j] = True
    return np.array(inner_indices),np.array(outer_indices)
rho,theta,phi = cart2spherical(coords)
inner,outer = separate_by_local_radii(coords[0,:],rho[0,:],5)
print('inner number: {}\n'.format(inner.size))
print('outer number: {}\n'.format(outer.size))
print('io ratio is: {:.3f}\n'.format(inner.size/outer.size))
print('radius is: {:.2f}\n.'.format((rho[0,inner].mean() + rho[0,outer].mean())/2))
exit
