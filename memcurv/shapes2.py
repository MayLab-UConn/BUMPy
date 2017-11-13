'''Redone shapes class

Will have the following:
    -argument parsing method
    -functions for generating shapes, with arguments and bilayer as paramaters
'''
import numpy as np
import rigid_body_transforms as rb
import nonrigid_coordinate_transformations as nrb
from copy import deepcopy,copy

class shape(object):
    '''
    Base class for all shapes to build on.
    '''
    @staticmethod
    def gen_shape():
        raise NotImplementedError

    @staticmethod
    def dimension_requirements():
        raise NotImplementedError

    @staticmethod
    def final_dimensions():
        raise NotImplementedError



class semisphere(shape):

    @staticmethod
    def dimension_requirements(r_sphere, buff=50):
        return np.array([(r_sphere + buff) * np.pi  , (r_sphere + buff) * np.pi ])

    @staticmethod
    def final_dimensions(r_sphere,buff=50):
        return np.array([2 * (r_sphere + buff), 2 *(r_sphere + buff), r_sphere +(2 * buff)])

    @ staticmethod
    def gen_shape(template_bilayer,zo,r_sphere,contains_hole=False,
               completeness=1):
        ''' returns molecules instance of semisphere'''
        # calculating slice radii
        slice_radius = np.pi * r_sphere / 2
        top_slice_radius = np.sqrt(2) * (r_sphere + zo)
        bot_slice_radius = np.sqrt(2) * (r_sphere - zo)
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
        # scale slices to slice_radius
        top_leaflet.coords = nrb.scale_coordinates_radial(top_leaflet.coords,
                             (slice_radius / top_slice_radius))
        bot_leaflet.coords = nrb.scale_coordinates_radial(bot_leaflet.coords,
                             (slice_radius / bot_slice_radius))
        # merge and transform slices
        top_leaflet.append_pdb(bot_leaflet)
        top_leaflet.coords = nrb.spherical_transform(top_leaflet.coords,r_sphere)
        return top_leaflet

class sphere(shape):

    @staticmethod
    def dimension_requirements(r_sphere, buff=50):
        return semisphere.dimension_requirements(r_sphere)

    @staticmethod
    def final_dimensions(r_sphere,buff=50):
        return np.array([2 * (r_sphere + buff)] * 3)

    @staticmethod
    def gen_shape(template_bilayer,zo,r_sphere,n_holes=0):
        top_half = semisphere.gen_shape(template_bilayer,zo,r_sphere,False)
        bot_half = copy(top_half)
        bot_half.coords = rb.rotate_coordinates(bot_half.coords,[180,0,0])
        top_half.append_pdb(bot_half,preserve_leaflets=True)
        return top_half

class cylinder(shape):

    @staticmethod
    def dimension_requirements(r_cylinder,l_cylinder,completeness=1, buff=50):
        return np.array([ l_cylinder + buff , (r_cylinder + buff) * 2 * np.pi * completeness ])

    @staticmethod
    def final_dimensions(r_cylinder,l_cylinder,buff=50):
        return np.array([l_cylinder, r_cylinder + buff, r_cylinder + buff])


    def gen_shape(template_bilayer,r_cylinder,l_cylinder,thickness,
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

        top_copy = copy(top_leaflet)
        top_copy.append_pdb(bot_leaflet)

        # scale coordinates
        top_leaflet.coords = nrb.scale_coordinates_rectangular(
                             top_leaflet.coords,[1,cylinder_slice_length/outer_slice_length])
        bot_leaflet.coords = nrb.scale_coordinates_rectangular(
                             bot_leaflet.coords,[1,cylinder_slice_length/inner_slice_length])

        top_leaflet.append_pdb(bot_leaflet)
        top_leaflet.coords = nrb.cylindrical_transform(rb.center_coordinates_3D(top_leaflet.coords),r_cylinder)
        return top_leaflet

class partial_torus(shape):
    @staticmethod
    def dimension_requirements(r_torus,r_tube, buff=50):
        return np.array([2 * (buff + r_torus + r_tube * np.pi) ] * 2)

    @staticmethod
    def final_dimensions(r_torus,r_tube,buff=50):
        return np.array([2 * (buff + r_torus + r_tube * np.pi) ] * 2 + [r_tube + buff])

    @staticmethod
    def gen_shape(template_bilayer,r_torus,r_tube,thickness,partial='full'):
        '''Makes a partial torus with given parameters

        The flat circular slice of a half_torus with torus R ranges from
        (R-r') to (R+r'), where r' = circumference of tube / 4

        Can get a quarter torus as well using partial='inner' for the shorter
        junction, partial ='outer' for the larger

        partial tori are oriented so that the cylindrical edges point DOWN,
        and has a minimum at 0 z

        I can't figure out good math on where to center / scale a quarter torus,
        in regards to the the effect of slicing on lipid ratios. Ie, can't just
        match length of slice to tube radius, as >>WHERE<< you slice changes
        total ratios, as opposed to cylinders (and spheres sorta). SO, will
        just do an additional slice of the half torus if quarter is selected
        '''
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
        #top_leaflet.write_pdb('torus_merge_notransform.pdb',position=False)
        # check quarter torus, use circular slice to cut off one side or other
        # the cutoff is r_torus. For inner, just take a circle that ends at r_torus
        # for outer, take circle larger than size of torus including everything,
        # then exclude up to r_torus
        top_leaflet.coords = nrb.toroidal_transform(top_leaflet.coords,r_torus,r_tube)
        #top_leaflet.write_pdb('torus_transform.pdb',position=False)
        if partial == 'inner':
            top_leaflet = top_leaflet.slice_pdb(top_leaflet.circular_slice(np.mean(top_leaflet.coords,axis=0),r_torus))
        elif partial == 'outer':
            top_leaflet = top_leaflet.slice_pdb(top_leaflet.circular_slice(np.mean(top_leaflet.coords,axis=0),r_torus+tube_circumference,exclude_radius=r_torus))

        return top_leaflet

class torus(shape):
    @staticmethod
    def dimension_requirements(r_torus,r_tube, buff=50):
        return partial_torus.dimension_requirements(r_torus,r_tube)

    @staticmethod
    def final_dimensions(r_torus,r_tube,buff=50):
        return np.array([2 * (buff + r_torus + r_tube * np.pi) ] * 2 + [2 *( r_tube + buff)])


    @staticmethod
    def gen_shape(template_bilayer,r_torus,r_tube,thickness,completeness=0.5):
        top_half = partial_torus(template_bilayer,r_torus,r_tube,thickness,partial='full')
        bot_half = copy(top_half)
        bot_half.coords = rb.rotate_coordinates(bot_half.coords,[180,0,0])
        top_half.append_pdb(bot_half)
        return top_half

class semicylinder_plane(shape):

    @staticmethod
    def dimension_requirements(r_cylinder,l_cylinder,r_junction,l_flat,buff=50):
        return

    @staticmethod
    def final_dimensions(r_cylinder,l_cylinder,r_junction,l_flat,buff=50):
        return


    @staticmethod
    def gen_shape(template_bilayer,zo,r_cylinder,l_cylinder,r_junction,l_flat):
        semicyl = shapes.cylinder(template_bilayer,r_cylinder,l_cylinder,
                                 thickness, completeness=0.5)
        junction = shapes.cylinder(template_bilayer,r_junction,l_cylinder,
                                 thickness, completeness=0.25)
        junction2 = copy(junction)

        # rotations and translations. 135 and 225 degrees gets junctions
        # rotated so that they max at 0 and taper to flat in the y direction
        # translation because they face the wrong direction and may not match
        # cylinder directions anyway
        junction.coords = rb.rotate_coordinates(junction.coords,[135,0,0]) - [0,r_junction + r_cylinder,0]
        junction2.coords = rb.rotate_coordinates(junction2.coords,[225,0,0]) + [0,r_junction + r_cylinder,0]

        ''' THIS IS INCOMPLETE AND WRONG'''
        flat_slice.coords = flat_slice.coords - np.mean(flat_slice.coords,axis=0)
        flat_slice.coords = flat_slice.coords + [0,r_junction+r_cylinder+(y_flat/2),-r_junction]
        semicyl.append_pdb(junction)
        semicyl.append_pdb(junction2)
        semicyl.append_pdb(flat_slice)
        return semicyl

class mitochondrion(shape):
    @staticmethod
    def dimension_requirements(r_cylinder,l_cylinder,r_junction,flat_dimension,buff=50):
        return

    @staticmethod
    def final_dimensions(r_cylinder,l_cylinder,r_junction,flat_dimension,buff=50):
        return


    @staticmethod
    def gen_shape(template_bilayer,zo,r_cylinder,l_cylinder,r_junction,flat_dimension):
        cyl = shapes.cylinder(template_bilayer,r_cylinder,l_cylinder,thickness,
                              completeness=1)
        cyl.coords = rb.rotate_coordinates(cyl.coords,[0,90,0])
        junction = shapes.partial_torus(template_bilayer,r_cylinder+r_junction,
                                    r_junction,thickness,partial='inner')
        junction.coords = junction.coords +  [0,0,l_cylinder/2]
        junction_2 = copy(junction)
        junction_2.coords = rb.rotate_coordinates(junction_2.coords,[180,0,0])
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

class elongated_vesicle(shape):
    @staticmethod
    def dimension_requirements(r_cylinder,l_cylinder,r_junction,l_flat,buff=50):
        return

    @staticmethod
    def final_dimensions(r_cylinder,l_cylinder,r_junction,l_flat,buff=50):
        return

    @staticmethod
    def gen_shape(template_bilayer,r,l_cylinder,thickness):
        ''' Two semispheres connected by a cylinder'''
        cyl = shapes.cylinder(template_bilayer,r,l_cylinder,thickness,completeness=1)
        semisphere1 = shapes.semisphere(template_bilayer,r,thickness,completeness=1)
        semisphere2 = copy(semisphere1)
        semisphere1.coords = rb.rotate_coordinates(semisphere1.coords,[0, 90,0])
        semisphere2.coords = rb.rotate_coordinates(semisphere2.coords,[0,270,0])
        semisphere1.coords[:,0] = semisphere1.coords[:,0] - l_cylinder / 2
        semisphere2.coords[:,0] = semisphere2.coords[:,0] + l_cylinder / 2
        cyl.append_pdb(semisphere1); cyl.append_pdb(semisphere2)
        return cyl
