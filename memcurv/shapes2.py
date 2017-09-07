'''Redone shapes class

Will have the following:
    -argument parsing method
    -functions for generating shapes, with arguments and bilayer as paramaters
'''
import numpy as np
import rigid_body_transforms as rb
import nonrigid_coordinate_transformations as nrb
from copy import deepcopy,copy

class shapes:
    def gen_shape(args,bilayer):
        if args.shape == 'sphere':
            print('Making sphere with:\nradius={:f} Angstroms\nthickness={:f} Angstroms\n'.format(args.r_sphere,args.thickness))
            return shapes.sphere(bilayer,args.r_sphere,args.thickness)
        elif args.shape == 'cylinder':
            print('Making sphere with:\nradius={:f} Angstroms\nlength={:f} Angstroms\nthickness={:f} Angstroms\n'.format(args.r_cylinder,args.cylinder_length,args.thickness))
            return shapes.cylinder(bilayer,args.r_cylinder,args.cylinder_length,args.thickness)
        else:
            pass

    def semisphere(template_bilayer,r_sphere,thickness,contains_hole=False,
                   completeness=1):
        ''' returns molecules instance of semisphere'''
        # calculating slice radii
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
        top_leaflet = template_bilayer.slice_pdb(np.intersect1d(
                      in_top_circular_slice, top_leaflet_ind))
        bot_leaflet = template_bilayer.slice_pdb(np.intersect1d(
                      in_bot_circular_slice, bot_leaflet_ind))
        top_leaflet.write_pdb('top_init.pdb',position=False)
        bot_leaflet.write_pdb('bot_init.pdb',position=False)

        #stop
        # scale slices to slice_radius
        top_leaflet.coords = nrb.scale_coordinates_radial(top_leaflet.coords,
                             (slice_radius / top_slice_radius))
        bot_leaflet.coords = nrb.scale_coordinates_radial(bot_leaflet.coords,
                             (slice_radius / bot_slice_radius))
        top_leaflet.write_pdb('top_scaled.pdb',position=False)
        bot_leaflet.write_pdb('bot_scaled.pdb',position=False)
        # merge and transform slices
        top_leaflet.append_pdb(bot_leaflet)
        top_leaflet.write_pdb('merge.pdb',position=False)
        top_leaflet.coords = nrb.spherical_transform(top_leaflet.coords,r_sphere)
        top_leaflet.write_pdb('transform.pdb',position=False)
        return top_leaflet

    def sphere(template_bilayer,r_sphere,thickness,n_holes=0):
        # check for holes
        if n_holes == 0:
            pass
        elif n_holes == 1:
            pass
        elif n_holes  == 2:
            pass

        top_half = shapes.semisphere(template_bilayer,r_sphere,thickness,False)
        bot_half = copy(top_half)
        bot_half.coords = rb.rotate_coordinates(bot_half.coords,[180,0,0])
        top_half.append_pdb(bot_half,preserve_leaflets=True)
        return top_half

    def cylinder(template_bilayer,r_cylinder,l_cylinder,thickness,
                completeness=1):
        '''Makes a cylinder with given parameters.

        Completeness=0.5 for semicylinder,
        completeness=0.25 for a lateral junction
        '''
        # calculate slice lengths
        cylinder_slice_length = 2 * np.pi * r_cylinder * completeness
        slice_origin = np.min(template_bilayer.coords,axis=0)[0:2]+40
        outer_slice_length = 2 * np.pi * (r_cylinder + (thickness/2)) * completeness
        inner_slice_length = 2 * np.pi * (r_cylinder - (thickness/2)) * completeness
        # calculate slice indices
        xvals = [slice_origin[0],slice_origin[0] + l_cylinder]
        yvals_outer = [slice_origin[1],slice_origin[1] + outer_slice_length]
        yvals_inner = [slice_origin[1],slice_origin[1] + inner_slice_length]
        in_top_slice =  template_bilayer.rectangular_slice(xvals,yvals_outer)
        in_bot_slice =  template_bilayer.rectangular_slice(xvals,yvals_inner)
        top_leaflet_ind = np.where(template_bilayer.metadata.leaflets == 1)[0]
        bot_leaflet_ind = np.where(template_bilayer.metadata.leaflets == 0)[0]
        # make slices
        top_leaflet = template_bilayer.slice_pdb(np.intersect1d(
                      in_top_slice, top_leaflet_ind))
        bot_leaflet = template_bilayer.slice_pdb(np.intersect1d(
                      in_bot_slice, bot_leaflet_ind))
        # scale coordinates
        top_leaflet.coords = nrb.scale_coordinates_rectangular(
                             top_leaflet.coords,[1,cylinder_slice_length/outer_slice_length])
        bot_leaflet.coords = nrb.scale_coordinates_rectangular(
                             bot_leaflet.coords,[1,cylinder_slice_length/inner_slice_length])
        top_leaflet.append_pdb(bot_leaflet)
        top_leaflet.coords = nrb.cylindrical_transform(rb.center_coordinates_3D(top_leaflet.coords),r_cylinder)
        return top_leaflet

    def partial_torus(template_bilayer,r_torus,r_tube,thickness,partial='full'):
        '''Makes a partial torus with given parameters

        The flat circular slice of a half_torus with torus R ranges from
        (R-r') to (R+r'), where r' = circumference of tube / 4

        Can get a quarter torus as well using partial='inner' for the shorter
        junction, partial ='outer' for the larger

        partial tori are oriented so that the cylindrical edges point DOWN,
        and has a minimum at 0 z
        '''
        tube_circumference = 2 *  np.pi * r_tube
        inner_tube_circumference = 2 * np.pi * (r_tube - (thickness/2))
        outer_tube_circumference = 2 * np.pi * (r_tube + (thickness/2))

        if partial == 'full':
            slice_min = r_torus - (tube_circumference / 4)
            slice_max = r_torus + (tube_circumference / 4)
            inner_slice_min = r_torus - (inner_tube_circumference/4)
            outer_slice_min = r_torus - (outer_tube_circumference/4)
            inner_slice_max = r_torus + (inner_tube_circumference/4)
            outer_slice_max = r_torus + (outer_tube_circumference/4)
        elif partial == 'inside':  # inside half, set slice max to radius
            slice_min = r_torus - (tube_circumference / 4)
            slice_max = r_torus
            inner_slice_min = r_torus - (inner_tube_circumference/4)
            outer_slice_min = r_torus - (outer_tube_circumference/4)
            inner_slice_max = r_torus
            outer_slice_max = r_torus
        elif partial == 'outside': # outside half, slice min to radius
            slice_min = r_torus
            slice_max = r_torus + (tube_circumference / 4)
            inner_slice_min = r_torus
            outer_slice_min = r_torus
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
        top_leaflet = template_bilayer.slice_pdb(np.intersect1d(
                      in_top_circular_slice, top_leaflet_ind))
        bot_leaflet = template_bilayer.slice_pdb(np.intersect1d(
                      in_bot_circular_slice, bot_leaflet_ind))

        # scale slices to slice_radius
        top_leaflet.coords = nrb.scale_coordinates_toroid(top_leaflet.coords,
                             [outer_slice_min,outer_slice_max],
                             [slice_min,slice_max])
        bot_leaflet.coords = nrb.scale_coordinates_toroid(bot_leaflet.coords,
                             [inner_slice_min,inner_slice_max],
                             [slice_min,slice_max])


        top_leaflet.append_pdb(bot_leaflet)
        top_leaflet.coords = nrb.toroidal_transform(top_leaflet.coords,r_torus,r_tube)
        return top_leaflet

    def torus(template_bilayer,r_torus,r_tube,thickness,completeness=0.5):
        top_half = shapes.half_torus(template_bilayer,r_torus,r_tube,thickness)
        bot_half = copy(top_half)
        bot_half.coords = rb.rotate_coordinates(bot_half.coords,[180,0,0])
        top_half.append_pdb(bot_half)
        return top_half

    def semicylinder_plane(template_bilayer,r_cylinder,l_cylinder,r_junction,
                           thickness,area_matching=True):
        semicyl = shapes.cylinder(template_bilayer,r_cylinder,l_cylinder,
                                 thickness, completeness=0.5)
        junction = shapes.cylinder(template_bilayer,r_junction,l_cylinder,
                                 thickness, completeness=0.25)
        junction2 = copy(junction)

        # rotations and translations. 135 and 225 degrees gets junctions
        # rotated so that they max at 0 and taper to flat in the y direction
        # translation because they face the wrong direction and may not match
        # cylinder directions anyway
        junction.coords = rb.rotate_coordinates(junction.coords,[135,0,0]) - [0,2*r_cylinder,0]
        junction2.coords = rb.rotate_coordinates(junction2.coords,[225,0,0]) + [0,2*r_cylinder,0]

        if area_matching:
            # flat section will be size of xdim, ydim is arc length of cylinder
            y_flat = np.pi * r_cylinder
            flat_slice = template_bilayer.slice_pdb(
                         template_bilayer.rectangular_slice(
                         [20 ,l_cylinder + 20],[20, y_flat + 20]))
            flat_slice.coords = flat_slice.coords - np.mean(flat_slice.coords,axis=0)
            flat_slice.coords = flat_slice.coords + [0,r_junction+r_cylinder+(y_flat/2),-r_junction]
        semicyl.append_pdb(junction)
        semicyl.append_pdb(junction2)
        semicyl.append_pdb(flat_slice)
        return semicyl

    def mitochondrion(template_bilayer,r_cylinder,l_cylinder,r_junction,
                          thickness,box_xy,area_matching=True):
        #'''
        cyl = shapes.cylinder(template_bilayer,r_cylinder,l_cylinder,thickness,
                              completeness=1)
        cyl.coords = rb.rotate_coordinates(cyl.coords,[0,90,0])
        junction = shapes.half_torus(template_bilayer,r_cylinder+r_junction,
                                     r_cylinder,thickness,partial='inside')
        junction.coords = junction.coords +  [0,0,l_cylinder/2]
        junction_2 = copy(junction)
        junction_2.coords = rb.rotate_coordinates(junction_2.coords,[180,0,0])
        #'''
        flat_bilayer = template_bilayer.slice_pdb(template_bilayer.rectangular_slice(
                       [20,box_xy+20],[20,box_xy+20],r_cylinder+r_junction))

        flat_bilayer.coords = flat_bilayer.coords -np.mean(flat_bilayer.coords,axis=0) + [0,0,(l_cylinder/2) + r_junction]
        flat_bilayer_2 = copy(flat_bilayer)
        flat_bilayer_2.coords = rb.rotate_coordinates(flat_bilayer.coords,[180,0,0])
        cyl.append_pdb(junction)
        cyl.append_pdb(flat_bilayer)
        cyl.append_pdb(junction_2)
        cyl.append_pdb(flat_bilayer_2)
        return cyl
