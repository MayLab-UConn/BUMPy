from argparse import ArgumentParser
import shapes2
from molecules import Molecules
import nonrigid_coordinate_transformations as nrb
import rigid_body_transforms as rb
import numpy as np
import pandas as pd

template_bilayer = Molecules('/home/kevin/hdd/Projects/Lipid_Diffusion/flat_template_bilayers/POPC.pdb')
outdir = '/home/kevin/hdd/Projects/Lipid_Diffusion/images/'

# shape generation
l_cylinder = 300
r_cylinder = 100
thickness = 2 * 11.6
completeness = 1
# ------------------------------------------------------------------------------
# cylinders
# ------------------------------------------------------------------------------
# typical cylindrical stuff
cylinder_slice_length = 2 * np.pi * r_cylinder * completeness
slice_origin = np.min(template_bilayer.coords,axis=0)[0:2]+40
outer_slice_length = 2 * np.pi * (r_cylinder + (thickness/2)) * completeness
inner_slice_length = 2 * np.pi * (r_cylinder - (thickness/2)) * completeness
xvals = [slice_origin[0],slice_origin[0] + l_cylinder]
yvals_outer = [slice_origin[1],slice_origin[1] + outer_slice_length]
yvals_inner = [slice_origin[1],slice_origin[1] + inner_slice_length]
in_top_slice =  template_bilayer.rectangular_slice(xvals,yvals_outer)
in_bot_slice =  template_bilayer.rectangular_slice(xvals,yvals_inner)
top_leaflet_ind = np.where(template_bilayer.metadata.leaflets == 1)[0]
bot_leaflet_ind = np.where(template_bilayer.metadata.leaflets == 0)[0]

''' Write whole bilayer pdb'''
#template_bilayer.write_pdb(outdir + 'whole_template.pdb',position=False)

top_leaflet = template_bilayer.slice_pdb(np.intersect1d(in_top_slice, top_leaflet_ind))
bot_leaflet = template_bilayer.slice_pdb(np.intersect1d(in_bot_slice, bot_leaflet_ind))

# copy is for making pdb with top and bottom of different sizes
top_copy    = template_bilayer.slice_pdb(np.intersect1d(in_top_slice, top_leaflet_ind))
'''write top and bot separately'''
#top_leaflet.write_pdb(outdir + 'top_initial_slice.pdb',position=False)
#bot_leaflet.write_pdb(outdir + 'bot_initial_slice.pdb',position=False)

top_copy.append_pdb(bot_leaflet)
'''Write top and bottom, different sizes'''
#top_copy.write_pdb(outdir + 'top_and_bot_no_transform.pdb',position='positive')
# now do transforms
top_leaflet.coords = nrb.scale_coordinates_rectangular(top_leaflet.coords,[1,cylinder_slice_length/outer_slice_length])
bot_leaflet.coords = nrb.scale_coordinates_rectangular(bot_leaflet.coords,[1,cylinder_slice_length/inner_slice_length])
top_leaflet.append_pdb(bot_leaflet)
'''Write top and bot in same box'''
#top_leaflet.write_pdb(outdir + 'top_and_bot_same_size.pdb',position='center')
top_leaflet.coords = nrb.cylindrical_transform(rb.center_coordinates_3D(top_leaflet.coords),r_cylinder)
'''Write cylinder'''
top_leaflet.write_pdb(outdir + 'cylindrical_transform.pdb',position='center')


# ------------------------------------------------------------------------------
# sphere
# ------------------------------------------------------------------------------
r_sphere = 100

slice_radius = np.pi * r_sphere / 2
top_slice_radius = np.pi * (r_sphere + (thickness/2)) / 2
bot_slice_radius = np.pi * (r_sphere - (thickness/2)) / 2
slice_origin = template_bilayer.gen_random_slice_point(top_slice_radius)
# calculate slice indices
in_top_circular_slice = template_bilayer.circular_slice(slice_origin,top_slice_radius)
in_bot_circular_slice = template_bilayer.circular_slice(slice_origin,bot_slice_radius)
top_leaflet_ind = np.where(template_bilayer.metadata.leaflets == 1)[0]
bot_leaflet_ind = np.where(template_bilayer.metadata.leaflets == 0)[0]
# make slices
top_leaflet = template_bilayer.slice_pdb(np.intersect1d(in_top_circular_slice, top_leaflet_ind))
bot_leaflet = template_bilayer.slice_pdb(np.intersect1d(in_bot_circular_slice, bot_leaflet_ind))
top_copy    = template_bilayer.slice_pdb(np.intersect1d(in_top_circular_slice, top_leaflet_ind))
top_copy.append_pdb(bot_leaflet)
top_copy.write_pdb(outdir + 'sphere_top_and_bot_no_transform.pdb',position='center')

# scale slices to slice_radius
top_leaflet.coords = nrb.scale_coordinates_radial(top_leaflet.coords,
                     (slice_radius / top_slice_radius))
bot_leaflet.coords = nrb.scale_coordinates_radial(bot_leaflet.coords,
                     (slice_radius / bot_slice_radius))
        # merge and transform slices
top_leaflet.append_pdb(bot_leaflet)
top_leaflet.write_pdb(outdir + 'sphere_top_and_bot_same_size.pdb',position='center')
top_leaflet.coords = nrb.spherical_transform(top_leaflet.coords,r_sphere)
top_leaflet.write_pdb(outdir + 'sphere_cylindrical_transform.pdb',position='center_xy')
# ------------------------------------------------------------------------------
# torus
# ------------------------------------------------------------------------------
r_tube = 100
r_torus = 200

tube_circumference = 2 *  np.pi * r_tube
inner_tube_circumference = 2 * np.pi * (r_tube - (thickness/2))
outer_tube_circumference = 2 * np.pi * (r_tube + (thickness/2))
slice_min = r_torus - (tube_circumference / 4)
slice_max = r_torus + (tube_circumference / 4)
inner_slice_min = r_torus - (inner_tube_circumference/4)
outer_slice_min = r_torus - (outer_tube_circumference/4)
inner_slice_max = r_torus + (inner_tube_circumference/4)
outer_slice_max = r_torus + (outer_tube_circumference/4)

slice_origin = np.mean(template_bilayer.coords,axis=0)[0:2]
# calculate slice indices
in_top_circular_slice = template_bilayer.circular_slice(slice_origin,
                        outer_slice_max,exclude_radius=outer_slice_min)
in_bot_circular_slice = template_bilayer.circular_slice(slice_origin,
                        inner_slice_max, exclude_radius=inner_slice_min)
# difference from sphere, exclude center
top_leaflet_ind = np.where(template_bilayer.metadata.leaflets == 1)[0]
bot_leaflet_ind = np.where(template_bilayer.metadata.leaflets == 0)[0]
top_leaflet = template_bilayer.slice_pdb(np.intersect1d(in_top_circular_slice, top_leaflet_ind))
top_copy =    template_bilayer.slice_pdb(np.intersect1d(in_top_circular_slice, top_leaflet_ind))
bot_leaflet = template_bilayer.slice_pdb(np.intersect1d(in_bot_circular_slice, bot_leaflet_ind))
top_copy.append_pdb(bot_leaflet)
top_copy.write_pdb(outdir + 'torus_top_and_bot_no_transform.pdb',position='center')
        # scale slices to slice_radius
top_leaflet.coords = nrb.scale_coordinates_toroid(top_leaflet.coords,
                     [outer_slice_min,outer_slice_max],
                     [slice_min,slice_max])
bot_leaflet.coords = nrb.scale_coordinates_toroid(bot_leaflet.coords,
                     [inner_slice_min,inner_slice_max],
                     [slice_min,slice_max])

top_leaflet.append_pdb(bot_leaflet)
top_leaflet.write_pdb(outdir + 'torus_top_and_bot_same_size.pdb',position='center')
top_leaflet.coords = nrb.toroidal_transform(top_leaflet.coords,r_torus,r_tube)
top_leaflet.write_pdb(outdir + 'torus_cylindrical_transform.pdb',position='center_xy')

top_slice = top_leaflet.slice_pdb(np.where(top_leaflet.metadata.leaflets == 1)[0])
bot_slice = top_leaflet.slice_pdb(np.where(top_leaflet.metadata.leaflets == 0)[0])
