# checks 1 FRAME for io ratio, used for validating initial shape
# because things are going wrong and it would be nice to
# have a quick summary accessible from the command line

import sys
import mdtraj as md
import numpy as np
infile = sys.argv[1]
shape  = sys.argv[2]
traj = md.load(infile)
traj.center_coordinates()
phoscoords = traj.xyz[0,traj.topology.select('name PO4'),:]
C4coords   = traj.xyz[0,traj.topology.select('name C4A'),:]

def cart2spherical(cart_coords,rho=None):
    x = cart_coords[:,0]
    y = cart_coords[:,1]
    z = cart_coords[:,2]
    rho = np.sqrt(x**2 + y**2 + z**2)
    return rho

def cart2polar(cart_coords):
    y = cart_coords[:,1]
    z = cart_coords[:,2]
    rho = np.sqrt(y**2 + z**2)
    return rho

if shape == 'sphere':
    rho_phos = cart2spherical(phoscoords)
    rho_C    = cart2spherical(C4coords)
elif shape == 'cylinder':
    rho_phos = cart2polar(phoscoords)
    rho_C    = cart2polar(C4coords)
inner = np.where(rho_phos <  rho_C)[0]
outer = np.where(rho_phos >= rho_C)[0]

oi_ratio = outer.size/inner.size
radius = (rho_phos[inner].mean() + rho_phos[outer].mean())/2


if shape == 'sphere':
    zo_nocorr = radius * (np.sqrt(oi_ratio) - 1) / (np.sqrt(oi_ratio) + 1)
    zo_corr = zo_nocorr * ( 1 + (1/4) * (( 2 * zo_nocorr / radius) **2))
elif shape == 'cylinder':
    zo_nocorr = radius * (oi_ratio -1) / (oi_ratio + 1)
    zo_corr = zo_nocorr * (1 + (1/2) * ((zo_nocorr / radius) ** 2))


print('inner number: {}\n'.format(inner.size))
print('outer number: {}\n'.format(outer.size))
print('io ratio is: {:.3f}\n'.format(inner.size/outer.size))
print('radius is: {:.2f}\n.'.format(radius))
print('zo, not corrected: {:.2f}\n'.format(zo_nocorr))
print('zo, corrected: {:.2f}\n'.format(zo_corr))
