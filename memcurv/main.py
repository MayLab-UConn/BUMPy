#!/usr/bin/env python

''' Main script for memcurv project'''

from argparse import ArgumentParser
import numpy as np
import pandas as pd
from time import time
from copy import copy

__version__ = '0.5'

# ------------------------------------------------------------------------------
# Rigid body stuff
# ------------------------------------------------------------------------------
class rb:
    @staticmethod
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

    @staticmethod
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
# ------------------------------------------------------------------------------
# NONRIGID BODY STUFF
# ------------------------------------------------------------------------------

class nrb:   # nonrigid body transformations
    @staticmethod
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

    @staticmethod
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
    @staticmethod
    def scale_coordinates_radial(coords,ratio):
        '''Coordinate scaling centered on 0'''
        meanvals = np.mean(coords,axis=0)
        coords = coords - [meanvals[0],meanvals[1],0]
        (theta,rho,z) = nrb.cart2pol(coords)
        rho = rho * ratio
        return nrb.pol2cart(theta,rho,z) + [meanvals[0],meanvals[1],0]

    @staticmethod
    def scale_coordinates_rectangular(coords,ratio):
        minvals = np.min(coords,axis=0)
        coords = coords - [minvals[0], minvals[1],0] # push to 0,0,0 for minima
        coords[:,0:2] = coords[:,0:2] * ratio
        return coords +  [minvals[0],minvals[1],0]

    @staticmethod
    def scale_coordinates_toroid(coords,current_range,new_range):
        '''Radial coordinate scaling from one range of spaces to another'''
        meanvals = np.mean(coords,axis=0)
        coords = coords - [meanvals[0],meanvals[1],0]
        (theta,rho,z) = nrb.cart2pol(coords)
        curr_range_size = current_range[1] - current_range[0]
        midpoint  = (current_range[0] + current_range[1]) / 2
        new_range_size  = new_range[1] - new_range[0]
        ratio = new_range_size / curr_range_size
        rho = (rho - midpoint)*ratio  + midpoint
        return nrb.pol2cart(theta,rho,z) + [meanvals[0],meanvals[1],0]


    # ------------------------------------------------------------------------------
    # Main curving transformations
    # ------------------------------------------------------------------------------
    @staticmethod
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

    @staticmethod
    def spherical_transform(xyz_coords,r,outer_leaflet='top'):
        ''' Transforms a circular segment of coordinates into a sphere with
            given radius r. Completeness of sphere will depend on radius of circular
            bilayer patch.

        '''
        xyz_coords = rb.center_coordinates_3D(xyz_coords)
        (theta,rho,z) = nrb.cart2pol(xyz_coords)

        radii = r + z
        arc_length_angle = rho / r
        rho_transform = radii * np.sin(arc_length_angle)
        z_transform   = radii * np.cos(arc_length_angle)
        return nrb.pol2cart(theta,rho_transform,z_transform)

    @staticmethod
    def toroidal_transform(xyz_coords,r_torus,r_tube):
        xyz_coords = rb.center_coordinates_3D(xyz_coords)
        (theta,rho,z) = nrb.cart2pol(xyz_coords)
        arc_length = rho - (r_torus - (np.pi * r_tube / 2))
        arc_length_angle = arc_length / r_tube
        radii = r_tube + z
        z_transform = radii * np.sin(arc_length_angle)
        rho_transform = r_torus + radii * np.sin(arc_length_angle - np.pi/2)

        return nrb.pol2cart(theta,rho_transform,z_transform)


# ------------------------------------------------------------------------------
# Molecules class
# ------------------------------------------------------------------------------

class Molecules:

    def __init__(self,infile=None,metadata=[],coords=[],boxdims=[]):
        '''Can initialize from file, or with manual inputs (like when slicing).
           Can also initialize with all objects blank, and wait for input
        '''
        if infile is not None:
            self.read_pdb(infile)
        else:
            self.coords   = coords
            self.metadata = metadata
            self.boxdims  = boxdims
    # -------------------------------------------------------------------------
    # Calculate and superficially change dataset properties
    # -------------------------------------------------------------------------

    def assign_leaflets(self):
        ''' Labels indices as top (1) or bottom (0) leaflet based on COM of
           entire residue
        '''
        coms = self.calc_residue_COMS()[:,2]
        start_indices = np.where(self.metadata.ressize > 0)[0]
        bilayer_com = np.mean(self.coords[:,2])
        leaflets = np.zeros(self.coords.shape[0],dtype=int)
        for i in range(coms.shape[0]):
            #leaflets[self.resid_list[i]] = int(coms[i,2] > bilayer_com)
            leaflets[start_indices[i]:start_indices[i]+self.metadata.ressize[start_indices[i]]] = int(coms[i] > bilayer_com)
        self.metadata.leaflets =leaflets


    def reorder_by_leaflet(self):
        ''' Switches up order of atoms so that top leaflet comes first,
            bot comes second
        '''
        new_index_order = np.append(np.where(self.metadata.leaflets == 1),
                                    np.where(self.metadata.leaflets == 0))
        self.coords = self.coords[new_index_order,:]
        self.metadata.index = new_index_order


    # -------------------------------------------------------------------------
    # adding and slicing pdb classes
    # -------------------------------------------------------------------------
    def append_pdb(self,new_pdb,preserve_leaflets=True):
        '''appends all information from new_pdb to end of current pdb '''
        # strings
        self.metadata = pd.concat((self.metadata,new_pdb.metadata))
        self.coords   = np.vstack((self.coords,  new_pdb.coords  ))
        # redo organizational arrays
        #self.reorganize_components(reset_leaflets= not preserve_leaflets)

    def slice_pdb(self,slice_indices,preserve_leaflets=True):
        ''' returns a new instance of current pdb class with sliced indices'''
        metadata_slice = self.metadata.ix[slice_indices]
        molecule_slice = Molecules(infile=None,
                                   metadata=metadata_slice,
                                   coords=self.coords[  slice_indices,:])

        #molecule_slice.reorganize_components(reset_leaflets= not preserve_leaflets)
        return molecule_slice

    def duplicate_laterally(self,nx,ny):

        # duplicate coordinates
        original_coords = np.copy(self.coords)
        nparts = original_coords.shape[0]
        new_coords = np.zeros((nparts*nx*ny,3))
        index = 0
        for i in range(nx):
            for j in range(ny):
                new_coords[index:index+nparts,:] = original_coords + [i*self.boxdims[0],j*self.boxdims[1],0]
                index+=nparts

        # duplicate metadata
        self.metadata = pd.DataFrame.from_dict({i:self.metadata[i].tolist() * nx * ny for i in self.metadata.columns})
        self.boxdims = [self.boxdims[0]*nx,self.boxdims[1]*ny,self.boxdims]
        self.coords = new_coords
    # -------------------------------------------------------------------------
    # calculating geometric slices
    # -------------------------------------------------------------------------
    def gen_slicepoint(self):
        ''' Within boundaries, randomize slicing region'''
        return np.min(self.coords[:,0:2],axis=0) + [20,20]

    def calc_residue_COMS(self):
        com_indices = np.where(self.metadata.ressize > 0)[0]
        sizes  = self.metadata.ressize[com_indices]
        res_coms = np.zeros((sizes.size,3)) # 0 will be empty
        count = 0
        for i in com_indices:
            res_coms[count,:] = np.mean(self.coords[i:i+sizes[i],:],axis=0)
            count +=1
        return res_coms

    def rectangular_slice(self,xvals,yvals,exclude_radius=0,partial_molecule='res_com'):
        '''Slices pdb to include only rectangular segment from x[0] to x[1] and
           y[0] to y[1]. Default is to exclude partial molecules, have option to
           include partial molecules or make whole and include.

           RETURNS INDICES, not actual slice
        '''
        res_starts = np.where(self.metadata.ressize > 0)[0]
        res_coms = self.calc_residue_COMS()
        indices_tokeep = []
        all_inrange = np.asarray(np.where( (res_coms[:,0] > xvals[0]) &
                                           (res_coms[:,0] < xvals[1]) &
                                           (res_coms[:,1] > yvals[0]) &
                                           (res_coms[:,1] < yvals[1]) ))[0]

        if exclude_radius > 0:
            centered_coords = res_coms - [np.mean(xvals),np.mean(yvals),0]
            theta,rho,z = nrb.cart2pol(centered_coords)
            rho_outside_exclusion = np.where(rho > exclude_radius)[0]
            all_inrange = np.intersect1d(all_inrange,rho_outside_exclusion)
            #stop
        if partial_molecule == 'exclude':
            print('excluding partial molecules')
            for i in range(len(self.resid_list)):
                resvals = np.array(self.resid_list[i])
                keep = True
                for j in resvals:
                    if not j in all_inrange: # every index of residue must be
                        keep = False         # in range
                        break
                if keep:
                    indices_tokeep =np.append(indices_tokeep,np.array(resvals))

        elif partial_molecule == 'res_com':
            print('Excluding based on residue COM cutoff')
            for i in all_inrange:
                indices_tokeep.extend(list(range(res_starts[i],res_starts[i] + self.metadata.ressize[res_starts[i]])))
        return np.asarray(indices_tokeep)

    def circular_slice(self,center,radius,exclude_radius=0,
                       partial_molecule='res_com'):
        indices_tokeep = []
        res_starts = np.where(self.metadata.ressize > 0)[0]
        res_coms = self.calc_residue_COMS()
        centered_coords = res_coms - [center[0],center[1],0]
        (theta,rho,z) = nrb.cart2pol(centered_coords)
        all_inrange = np.where((rho <= radius) & (rho>= exclude_radius))[0]
        if partial_molecule == 'exclude':
            for i in range(res_coms.size):
                resvals = np.array(self.resid_list[i])
                keep = True
                for j in resvals:
                    if not j in all_inrange: # every index of residue must be
                        keep = False         # in range
                        break
                if keep:
                    indices_tokeep =np.append(indices_tokeep,np.array(resvals))
        elif partial_molecule == 'res_com':
            print('Excluding based on residue COM cutoff')
            for i in all_inrange:
                indices_tokeep.extend(list(range(res_starts[i],res_starts[i] + self.metadata.ressize[res_starts[i]])))
        return np.asarray(indices_tokeep)

    # -------------------------------------------------------------------------
    # file i/o
    # -------------------------------------------------------------------------
    def read_pdb(self,pdbfile,reorganize=False):
        '''Read input pdb
        '''
        print('Reading in PDB file')
        t = time()
        with open(pdbfile,"r") as fid:
            # initialize temporary variables
            xcoord, ycoord, zcoord       = [],[],[]                   # 3D coordinates
            atomtype, atomname, resname  = [],[],[]                    # invariant labels
            curr_res , prev_res,ressize  = [],[],[]                    # residue indexing
            atomcount = 0                                             # counters

            for pdb_line in fid:
                if pdb_line.startswith("ATOM") or pdb_line.startswith("HETATM"):
                    # strings
                    atomtype.append(pdb_line[ 0: 6])
                    atomname.append(pdb_line[12:16])
                    resname.append( pdb_line[17:21])

                    # counting residues
                    if not prev_res: # first iteration
                        prev_res = pdb_line[22:26]
                    curr_res = pdb_line[22:26]
                    if curr_res == prev_res:
                        atomcount += 1
                    else:
                        ressize.append(atomcount)
                        ressize += [0] * (atomcount - 1)
                        atomcount = 1
                        prev_res = curr_res

                    # floats for coordinate array
                    xcoord.append(float(pdb_line[30:38]))
                    ycoord.append(float(pdb_line[38:46]))
                    zcoord.append(float(pdb_line[46:54]))
                elif pdb_line.startswith('CRYST1'):
                    self.boxdims = [float(i) for i in pdb_line.split()[1:4]]

            ressize.append(atomcount) # one more for final residue
            ressize += [0] * (atomcount - 1)
            # close file here

        zero_array = np.zeros(len(ressize),dtype=int)
        self.metadata = pd.DataFrame(data={
                                          'ressize':ressize,
                                         'leaflets':zero_array,
                                         'atomtype':atomtype,
                                         'atomname':atomname,
                                          'resname':resname})
        self.coords = np.array((xcoord,ycoord,zcoord)).T
        # assigning resid lengths and leaflets
        self.assign_leaflets()
        print('Finished reading in PDB file with {} atoms,time elapsed = {:.1f} seconds'.format(self.coords.shape[0],time()-t))

    def write_pdb(self,outfile,position='positive',reorder=True):
        print('Writing out PDB file')
        t = time()

        if reorder:
            self.reorder_by_leaflet()

        '''Outputs to pdb file, CRYST1 and ATOM lines only'''

        # default is to bring everything to positive regime, allows for up to
        #9999 angstroms in PDB file format
        if position == 'positive':
            out_coords = self.coords -np.min(self.coords,axis=0)
        elif position == 'positive_xy':
            out_coords = self.coords
            out_coords[:,0:2] = out_coords[:,0:2] - np.min(out_coords[:,0:2],axis=0)
        elif position == 'center':
            out_coords = self.coords - np.mean(self.coords,axis=0)
        elif position == 'center_xy':
            out_coords = self.coords
            out_coords[:,0:2] = out_coords[:,0:2] - np.mean(out_coords[:,0:2],axis=0)
        else:
            out_coords = self.coords

        nparts = self.coords.shape[0]
        # write out box dims


        with open(outfile,'w') as fout:
            fout.write('CRYST1{0:9.3f}{1:9.3f}{2:9.3f}{3:7.2f}{3:7.2f}{3:7.2f}\n'.format(*self.boxdims,90))

            # calculating residue and atom numbers
            resid = np.zeros(nparts,dtype=int)
            count = 1
            ressize = self.metadata.ressize.tolist()
            for i in range(nparts):
                size = ressize[i]
                if size > 0:
                    resid[i:i + size] = count
                    count += 1
            resid = np.mod(resid,10000)
            atomno = np.mod(np.arange(1,nparts+1),100000)

            fout.writelines(["{:6s}{:5d} {:4s} {:4s}{:5d}    {:8.3f}{:8.3f}{:8.3f}\n".format(
                                         i[0],i[1],i[2],i[3],i[4], i[5],i[6], i[7])  for i in zip(
                                         self.metadata.atomtype,
                                         atomno,
                                         self.metadata.atomname,
                                         self.metadata.resname,
                                         resid,
                                         out_coords[:,0],out_coords[:,1],out_coords[:,2] )])


        print('Finished writing PDB file with {} atoms,time elapsed = {:.1f} seconds'.format(nparts,time()-t))

    def write_topology(self,outfile):
        '''Writes out simple topology file (.top)'''
        with open(outfile,'w') as fout:
            fout.write("\n\n\n[ system ]\nmemcurv system\n\n[ molecules ]\n")


            reslist = self.metadata.resname[np.where(self.metadata.resid_length > 0)[0]]
            prev_res = reslist[0]
            counter = 0
            for res in reslist:
                if res == prev_res:
                    counter += 1
                else:
                    fout.write("{:4s} {:d}\n".format(prev_res,counter))
                    counter = 1
                    prev_res = res
            fout.write("{:4s} {:d}\n".format(prev_res,counter))

    def write_index(self,outfile,special_groups=None):
        '''Writes out index file (.ndx) with the following (hopefully useful)
           fields:
                    -system
                    -top_leaflet
                    -top_leaflet_component_1
                    -top_leaflet_component_2...
                    -bot_leaflet
                    -bot_leaflet_componenet_1...
        '''
        def write_index_unit(print_obj,name,indices):
            print_obj.write("[ {:s} ]\n".format(name))
            count = 1
            for i in indices:
                print_obj.write("{:d} ".format(i))
                count += 1
                if count > 15:  # 15 numbers per line
                    count = 1
                    print_obj.write("\n")
            print_obj.write("\n\n")

        with open(outfile,'w') as fout:
            top_atomno = np.where(self.metadata.leaflets == 1)[0] + 1
            bot_atomno = np.where(self.metadata.leaflets == 0)[0] + 1
            write_index_unit(fout,"system",self.metadata.atomno)
            write_index_unit(fout,"top leaflet",top_atomno)
            write_index_unit(fout,"bot leaflet",bot_atomno)
            #if special_groups is not None:
            #        for name in special_groups:
            #            name_indices =
            #for i in list(set(self.resname)):
                # this is where individual residues would go


# ------------------------------------------------------------------------------
# SHAPE REPOSITORY
# ------------------------------------------------------------------------------
class shapes:
    class shape:
        '''
        Base class for all shapes to build on. Must have all 3 methods
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
            slice_origin = np.mean(coords[:,0:2],axis=0)
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
            return shapes.semisphere.dimension_requirements(r_sphere)

        @staticmethod
        def final_dimensions(r_sphere,buff=50):
            return np.array([2 * (r_sphere + buff)] * 3)

        @staticmethod
        def gen_shape(template_bilayer,zo,r_sphere,n_holes=0):
            top_half = shapes.semisphere.gen_shape(template_bilayer,zo,r_sphere,False)
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


        def gen_shape(template_bilayer,zo,r_cylinder,l_cylinder,completeness=1):
            '''Makes a cylinder with given parameters.

            Completeness=0.5 for semicylinder,
            completeness=0.25 for a lateral junction
            '''
            # calculate slice lengths
            cylinder_slice_length = 2 * np.pi * r_cylinder * completeness
            slice_origin = np.min(template_bilayer.coords,axis=0)[0:2]+40
            outer_slice_length = 2 * np.pi * (r_cylinder + zo) * completeness
            inner_slice_length = 2 * np.pi * (r_cylinder - zo) * completeness
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
        def gen_shape(template_bilayer,zo,r_torus,r_tube,partial='full'):
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
            inner_tube_circumference = 2 * np.pi * (r_tube - zo)
            outer_tube_circumference = 2 * np.pi * (r_tube + zo)

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
            return shapes.partial_torus.dimension_requirements(r_torus,r_tube)

        @staticmethod
        def final_dimensions(r_torus,r_tube,buff=50):
            return np.array([2 * (buff + r_torus + r_tube * np.pi) ] * 2 + [2 *( r_tube + buff)])


        @staticmethod
        def gen_shape(template_bilayer,zo,r_torus,r_tube,completeness=0.5):
            top_half = shapes.partial_torus(template_bilayer,zo,r_torus,r_tube,partial='full')
            bot_half = copy(top_half)
            bot_half.coords = rb.rotate_coordinates(bot_half.coords,[180,0,0])
            top_half.append_pdb(bot_half)
            return top_half

    class semicylinder_plane(shape):

        @staticmethod
        def dimension_requirements(r_cylinder,l_cylinder,r_junction,l_flat,buff=50):
            cyldims = shapes.cylinder.dimension_requirements(r_cylinder,l_cylinder,completeness=0.5)
            flatdims = np.array([l_cylinder,l_flat])
            jdims  =  shapes.cylinder.dimension_requirements(r_junction,l_cylinder,completeness=0.5)
            return np.array([max([cyldims[0],flatdims[0],jdims[0]]), max([cyldims[1],flatdims[1],jdims[1]])])

        @staticmethod
        def final_dimensions(r_cylinder,l_cylinder,r_junction,l_flat,buff=50):
            return np.array([l_cylinder,2 * (r_cylinder + r_junction) + l_flat,r_cylinder + r_junction + 2* buff])


        @staticmethod
        def gen_shape(template_bilayer,zo,r_cylinder,l_cylinder,r_junction,l_flat):
            semicyl   = shapes.cylinder.gen_shape(template_bilayer, zo, r_cylinder, l_cylinder, completeness=0.5)
            junction  = shapes.cylinder.gen_shape(template_bilayer, zo, r_junction, l_cylinder, completeness=0.25)
            junction2 = copy(junction)

            # rotations and translations. 135 and 225 degrees gets junctions
            # rotated so that they max at 0 and taper to flat in the y direction
            # translation because they face the wrong direction and may not match
            # cylinder directions anyway
            junction.coords = rb.rotate_coordinates(junction.coords,[135,0,0]) - [0,r_junction + r_cylinder,0]
            junction2.coords = rb.rotate_coordinates(junction2.coords,[225,0,0]) + [0,r_junction + r_cylinder,0]

            ''' THIS IS INCOMPLETE AND WRONG'''
            slice_origin = template_bilayer.gen_slicepoint()
            flat_slice  = template_bilayer.slice_pdb(template_bilayer.rectangular_slice([slice_origin[0],slice_origin[0]+l_cylinder],[slice_origin[1],slice_origin[1]+l_flat]))
            flat_slice.coords = flat_slice.coords - np.mean(flat_slice.coords,axis=0) + [0,r_junction+r_cylinder+(l_flat/2),-r_junction]
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
            cyl = shapes.cylinder.gen_shape(template_bilayer,zo,r_cylinder,l_cylinder,
                                  completeness=1)
            cyl.coords = rb.rotate_coordinates(cyl.coords,[0,90,0])
            junction = shapes.partial_torus(template_bilayer,zo,r_cylinder+r_junction,
                                        r_junction,partial='inner')
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
        def gen_shape(template_bilayer,zo,r,l_cylinder):
            ''' Two semispheres connected by a cylinder'''
            cyl = shapes.cylinder(template_bilayer,zo,r,l_cylinder,completeness=1)
            semisphere1 = shapes.semisphere(template_bilayer,zo,r,completeness=1)
            semisphere2 = copy(semisphere1)
            semisphere1.coords = rb.rotate_coordinates(semisphere1.coords,[0, 90,0])
            semisphere2.coords = rb.rotate_coordinates(semisphere2.coords,[0,270,0])For strings, consider whether you really need them to appear in your source code in the first place. If you have plans to localize the script (i.e. translate the messages), you're almost certainly going to want them in a separate file anyway.
That said, wrapping things honestly isn't that bad. (And please don't use implicit concatenation; it's un-Pythonic - Explicit is better than implicit - and might conceivably go away; that's been proposed before.)
This is my preferred style - using the natural "concatenation" of the print function (I don't use 2.x any more, but it still works there):
            semisphere1.coords[:,0] = semisphere1.coords[:,0] - l_cylinder / 2
            semisphere2.coords[:,0] = semisphere2.coords[:,0] + l_cylinder / 2
            cyl.append_pdb(semisphere1); cyl.append_pdb(semisphere2)
            return cyl


def parse_command_lines():
    ''' Parses command line for parameters, returns parsed arguments '''
    prog_name =  'JASON'
    prog_description = 'Creating curved membrane systems with arbitrary geometry' + \
    ' and lipid composition, using a pivotal plane-based approach to appropriately'+ \
    ' match inter-leaflet area differences'

    geometry_description = 'Geometric arguments should be added as a series of ' + \
    'argument:value pairs separated by a colon.See the README for a list of required geometric arguments for a given shape.'


    parser = ArgumentParser(prog=prog_name,description=prog_description,
                            add_help=False,allow_abbrev=False)
    # groups
    required_inputs     = parser.add_argument_group('required inputs')
    geometric_arguments = parser.add_argument_group('geometric arguments',
                          geometry_description)
    optional_arguments  = parser.add_argument_group('optional arguments')
    output_arguments    = parser.add_argument_group('output arguments')
    # mandatory input
    required_inputs.add_argument('-s',help='Shape to make - see manual for a list of shapes',metavar='')
    required_inputs.add_argument('-f',help='Flat bilayer template to be used as a template',metavar='')
    required_inputs.add_argument('-z',type=float,help='Location of the pivotal plane (nm)',metavar='')

    # geometry
    geometric_arguments.add_argument('-g',nargs='*',required=True,help='Format is arg:value, ie r_cylinder:10, l_cylinder:20, ... for every geometric parameter in shape',metavar='')


    # optional arguments
    optional_arguments.add_argument('-h','--help', action='help', help='show this help message and exit')
    optional_arguments.add_argument('-outer',default='top',help='By default, top leaflet = outer leaflet. Set to "bot" to invert',metavar='')
    optional_arguments.add_argument('-uapl',metavar='',help='Slice top bilayer to achieve a specific area per lipid in final shape - not yet implemented')
    optional_arguments.add_argument('-lapl',metavar='',help='Slice bottom bilayer to achieve a specific area per lipid in final shape - not yet implemented')
    # output files
    output_arguments.add_argument('-o',help='Output structure - only PDBs for now',default='confout.pdb',metavar='')
    output_arguments.add_argument('-p',help='Simple .top topology file',metavar='') # optional
    output_arguments.add_argument('-n',help='Simple .ndx index file, separating leaflets',metavar='') # optional
    return parser.parse_args()

def display_parameters(cl_args):
    '''Displays selected parameters upon command line execution'''
    pass


# -----------------------------------------------------------------------------
# Command line start
# -----------------------------------------------------------------------------
def main():
    print('Parsing command line arguments\n')
    args = parse_command_lines()
    display_parameters(args)  # show user what they selected

    # parse arguments
    geometric_args = {garg.split(':')[0]:float(garg.split(':')[1]) for garg in args.g } # use a comprehension
    zo = args.z
    #shape_tobuild = getattr(shapes2.shapes,args.shape) This will be correct when compiled together
    shape_tobuild = getattr(shapes,args.s)

    # adjust size of template bilayer
    template_bilayer = Molecules(infile=args.f)

    if args.outer == 'bot':
        template_bilayer.coords = rb.rotate_coordinates(template_bilayer.coords,[180,0,0])
        template_bilayer.metadata.leaflets = np.invert(template_bilayer.metadata.leaflets)

    mult_factor = np.ceil( shape_tobuild.dimension_requirements(**geometric_args)/template_bilayer.boxdims[0:2]).astype(int)
    template_bilayer.duplicate_laterally(*mult_factor)

    # construct the shape
    shape = shape_tobuild.gen_shape(template_bilayer,zo,**geometric_args)
    shape.boxdims = shape_tobuild.final_dimensions(**geometric_args)
    # file output
    shape.write_pdb(args.o)
    if args.p:
        shape.write_topology(args.p)
    if args.n:
        shape.write_index(args.n)

if __name__ == '__main__':
    main()
