#!/usr/bin/env python

''' Main script for BUMPY project'''

import numpy as np
import inspect
import sys
from argparse import ArgumentParser
from time import time
from copy import deepcopy


__version__ = '0.5'


def cart2pol(cart_coords):
    ''' Converts cartesian to polar coordinates. Acts on first 2 columns (xy)
        returns tuple of (theta,rho,z)
    '''
    x, y, z = cart_coords[:, 0], cart_coords[:, 1], cart_coords[:, 2]
    theta    = np.arctan2(y, x)
    rho  = np.sqrt(x**2 + y**2)
    return(theta, rho, z)


def pol2cart(theta, rho, z):
    ''' Converts polar to cartesian coordinates, outputs as nparts * xyz array
    '''
    x = rho * np.cos(theta)
    y = rho * np.sin(theta)
    cart_coords = np.stack((x, y, z), axis=1)
    return cart_coords


# ------------------------------------------------------------------------------
# Molecules class
# ------------------------------------------------------------------------------

class Metadata:
    ''' Structure containing numpy arrays that contains atomnames, resnames, leaflet info, molecule info
        atomname and resname are dtype <U4, leaflets and ressize are dtype int.

        Supports appending, duplication, slicing, and reordering operations
    '''

    def __init__(self, atomname=[], resname=[], leaflets=[], ressize=[]):
        self.atomname = atomname
        self.resname  = resname
        self.leaflets = leaflets
        self.ressize  = ressize

    def append(self, other):
        self.atomname = np.concatenate((self.atomname, other.atomname))
        self.resname  = np.concatenate((self.resname,  other.resname))
        self.leaflets = np.concatenate((self.leaflets, other.leaflets))
        self.ressize  = np.concatenate((self.ressize,  other.ressize))

    def reorder(self, indices):
        self.atomname = self.atomname[indices]
        self.resname  = self.resname[indices]
        self.leaflets = self.leaflets[indices]
        self.ressize  = self.ressize[indices]

    def duplicate(self, n):
        self.atomname = np.tile(self.atomname, n)
        self.resname  = np.tile(self.resname,  n)
        self.leaflets = np.tile(self.leaflets, n)
        self.ressize  = np.tile(self.ressize,  n)

    def slice(self, indices):
        return Metadata(atomname=self.atomname[indices], resname=self.resname[indices],
                        leaflets=self.leaflets[indices], ressize=self.ressize[indices])


class Molecules:
    '''
        Giant bloated mess of a class containing the coordinates and type data for flat input bilayers as well as
        the shapes that they are transformed into. I probably could have organized this whole setup better, but as it
        stands this class contains all the i/o functionality, and the functions used to slice and merge shapes, which
        are then used by the shapes class.

        There are 2 major and 1 minor attributes in this class. Coords contain the 3D coordinates of a shape with
        dimensions of n-particles * 3. Metadata is an object containing atom names, residue names, a record of which
        leaflet the lipid belongs to, and a "ressize" array which I use to keep track of separate lipid molecules. For
        each atom in the system, if it is the first atom sequentially in the molecule, ressize for that atom is set to
        the number of atoms in the molecule, and every other atom in the molecule is set to 0. Boxdims are the last
        attribute, list of xyz box dimensions.
    '''

    def __init__(self, infile=None, metadata=[], coords=[], boxdims=[], ignore=[]):
        '''Can initialize from file, or with manual inputs (like when slicing).
           Can also initialize with all objects blank, and wait for input
        '''
        if infile is not None:
            self.read_input(infile, ignore=ignore)
        else:
            self.coords   = coords
            self.metadata = metadata
            self.boxdims  = boxdims

    # -----------------------------------------------------------------------------------------
    # Geometric transformations
    # -----------------------------------------------------------------------------------------
    def center_on_zero(self, ztype='bilayer_interface'):
        meanvals = self.coords.mean(axis=0)
        if ztype == 'bilayer_interface':    # other option is 'mean', which is just the simple average done above
            meanvals[2] = self.get_bilayer_center()
        self.coords -= meanvals

    def translate(self, offset):
        ''' offset is [x,y,z] translation'''
        self.coords += offset

    def rotate(self, rotation_angles, com=False, unit='degrees'):
        ''' Rotational transform of coordinate system. Default rotation center is
            about the origin (com=False). If com is set to True, rotation will occur
            about center of mass of coordinates (with no weighting of particles).
            Default input is in degrees, other option is 'radians'
        '''
        if unit == 'degrees':
            rotation_angles = np.radians(rotation_angles)
        if com:
            init_com = self.coords.mean(axis=0)
            self.coords -= init_com

        r_x = np.array([[1,                       0,                              0],
                        [0, np.cos(rotation_angles[0]), -np.sin(rotation_angles[0])],
                        [0, np.sin(rotation_angles[0]), np.cos(rotation_angles[0])]])

        r_y = np.array([[np.cos(rotation_angles[1]) , 0, np.sin(rotation_angles[1])],
                        [0                       ,    1,                        0  ],
                        [-np.sin(rotation_angles[1]), 0, np.cos(rotation_angles[1])]])

        r_z = np.array([[np.cos(rotation_angles[2]), -np.sin(rotation_angles[2]), 0],
                        [np.sin(rotation_angles[2]),  np.cos(rotation_angles[2]), 0],
                        [0,                       0,                              1]])

        rotmat = r_x.dot(r_y).dot(r_z)       # net rotation matrix
        self.coords = np.dot(self.coords, rotmat)   # do rotation
        if com:
            self.coords += init_com

    def scale_coordinates_radial(self, ratio):
        '''Coordinate scaling centered on 0'''
        meanvals = np.mean(self.coords, axis=0)
        self.coords -= [meanvals[0], meanvals[1], 0]
        (theta, rho, z) = cart2pol(self.coords)
        rho = rho * ratio
        self.coords = pol2cart(theta, rho, z)

    def scale_coordinates_rectangular(self, ratio):
        minvals = np.min(self.coords, axis=0)
        self.coords -= [minvals[0], minvals[1], 0]  # push to 0,0,0 for minima
        self.coords[:, 0:2] = self.coords[:, 0:2] * ratio

    def scale_coordinates_toroidal(self, current_range, new_range):
        ''' Radial coordinate scaling from one range of spaces to another. This function is a mess but I don't feel
            like reworking it for now
        '''
        meanvals = np.mean(self.coords, axis=0)
        self.coords -= [meanvals[0], meanvals[1], 0]
        (theta, rho, z) = cart2pol(self.coords)                      # center and turn to polar coordinates
        curr_range_size = current_range[1] - current_range[0]
        midpoint  = (current_range[0] + current_range[1]) / 2
        new_range_size  = new_range[1] - new_range[0]
        ratio = new_range_size / curr_range_size
        rho = (rho - midpoint) * ratio  + midpoint
        self.coords = pol2cart(theta, rho, z) + [meanvals[0], meanvals[1], 0]

    def cylindrical_transform(self, r, outer_leaflet='top'):
        ''' Transforms a rectangular segment of coordinates into a cylinder with
            given radius r. Completeness of cylinder will depend on length of y
            dimension, x dimension is long axis of cylinder.
        '''
        self.center_on_zero()
        radii = r + self.coords[:, 2]                                           # r is the desired radius, and is used
        arc_length_angle = self.coords[:, 1]  /  r                              # to convert cartesian coordinate to an
        y_transform = radii * np.sin(arc_length_angle)                          # angle, whereas the "radii" variable is
        z_transform = radii * np.cos(arc_length_angle)                          # used for the actual transformation
        self.coords = np.stack((self.coords[:, 0], y_transform, z_transform), axis=1)

    def spherical_transform(self, r, outer_leaflet='top'):
        ''' Transforms a circular segment of coordinates into a sphere with
            given radius r. Completeness of sphere will depend on radius of circular
            bilayer patch.
        '''
        self.center_on_zero()
        (theta, rho, z) = cart2pol(self.coords)
        radii = r + z
        arc_length_angle = rho / r
        rho_transform = radii * np.sin(arc_length_angle)
        z_transform   = radii * np.cos(arc_length_angle)
        self.coords = pol2cart(theta, rho_transform, z_transform)

    def toroidal_transform(self, r_torus, r_tube):
        self.center_on_zero()
        (theta, rho, z) = cart2pol(self.coords)
        arc_length = rho - (r_torus - (np.pi * r_tube / 2))
        arc_length_angle = arc_length / r_tube
        radii = r_tube + z
        z_transform = radii * np.sin(arc_length_angle)
        rho_transform = r_torus + radii * np.sin(arc_length_angle - np.pi / 2)
        self.coords = pol2cart(theta, rho_transform, z_transform)

    # -------------------------------------------------------------------------
    # Calculate and superficially change dataset properties
    # -------------------------------------------------------------------------

    def assign_leaflets(self):
        ''' Labels indices as top (1) or bottom (0) leaflet based on COM of
           entire residue
        '''
        coms = self.calc_residue_COMS()[:, 2]
        start_indices = np.where(self.metadata.ressize > 0)[0]
        bilayer_com = np.mean(self.coords[:, 2])
        leaflets = np.zeros(self.coords.shape[0], dtype=int)
        for i in range(coms.shape[0]):
            ind = start_indices[i]
            leaflets[ind:ind + self.metadata.ressize[ind]] = int(coms[i] > bilayer_com)
        self.metadata.leaflets = leaflets

    def get_bilayer_center(self, method='per_residue', nparts=None):
        ''' The center of a bilayer is easy to calculate in a symmetric system but not in an asymmetric system, as
            the thickenss of individual monolayers may differ, leaving no good point of comparison between the two.
            The solution here is to calculate the closest points to the bilayer center for each leaflet, and take the
            average of those
        '''
        if method == 'per_residue':
            top_min = []
            bot_max = []
            for i in np.where(self.metadata.ressize > 0)[0]:
                if self.metadata.leaflets[i]  == 1:
                    top_min.append(self.coords[i:i + self.metadata.ressize[i], 2].min())
                else:
                    bot_max.append(self.coords[i:i + self.metadata.ressize[i], 2].max())

        elif method == 'first_nparts':
            if not nparts:
                nparts = [np.sum((self.metadata.leaflets == 1) & (self.metadata.ressize > 0)),  # default 1 per residue
                          np.sum((self.metadata.leaflets == 0) & (self.metadata.ressize > 0)) ]

            top_min = np.sort(self.coords[self.metadata.leaflets == 1, 2])[:nparts[0]]
            bot_max = np.sort(self.coords[self.metadata.leaflets == 0, 2])[-nparts[1]:]
        return (np.array(top_min).mean() + np.array(bot_max).mean()) / 2

    def reorder_by_leaflet(self):
        ''' Switches up order of atoms so that top leaflet comes first,
            bot comes second
        '''
        new_index_order = np.append(np.where(self.metadata.leaflets == 1), np.where(self.metadata.leaflets == 0))
        self.coords = self.coords[new_index_order, :]
        self.metadata.reorder(new_index_order)

    def reorder_within_leaflet(self):
        outer_resnames = np.unique(self.metadata.resname[self.metadata.leaflets == 1])
        inner_resnames = np.unique(self.metadata.resname[self.metadata.leaflets == 0])
        new_order = []
        for resname in outer_resnames:
            new_order += list(np.where((self.metadata.resname == resname) & (self.metadata.leaflets == 1))[0])
        for resname in inner_resnames:
            new_order += list(np.where((self.metadata.resname == resname) & (self.metadata.leaflets == 0))[0])
        self.coords = self.coords[new_order, :]
        self.metadata.reorder(new_order)

    def reorder_with_dummies_in_back(self, dummy_name):
        not_dummy_indices = np.where(self.metadata.resname != dummy_name)[0]
        dummy_indices     = np.where(self.metadata.resname == dummy_name)[0]
        new_order = np.append(not_dummy_indices, dummy_indices)
        self.coords = self.coords[new_order, :]
        self.metadata.reorder(new_order)

    # -------------------------------------------------------------------------
    # adding and slicing pdb classes
    # -------------------------------------------------------------------------
    def append(self, new_pdb, preserve_leaflets=True):
        '''appends all information from new_pdb to end of current pdb. Doesn't change boxdims at this point '''
        # strings
        self.metadata.append(new_pdb.metadata)
        self.coords = np.vstack((self.coords, new_pdb.coords))

    def slice_pdb(self, slice_indices):
        ''' returns a new instance of current pdb class with sliced indices'''
        metadata_slice = self.metadata.slice(slice_indices.astype(bool))
        molecule_slice = Molecules(metadata=metadata_slice, coords=self.coords[slice_indices.astype(bool), :],
                                   boxdims=[0, 0, 0])
        return molecule_slice

    def duplicate_laterally(self, nx, ny):
        ''' Duplicates a flat bilayer nx times in the x dimension, ny times in the y dimension, similarly to
            gmx genconf. Works really fast, so it's much better to load in a small bilayer and multiply it a lot,
            rather than load in a giant bilayer.
        '''
        original_coords = np.copy(self.coords)
        nparts = original_coords.shape[0]
        new_coords = np.zeros((nparts * nx * ny, 3))
        index = 0
        for i in range(nx):
            for j in range(ny):
                new_coords[index:index + nparts, :] = original_coords + [i * self.boxdims[0], j * self.boxdims[1], 0]
                index += nparts

        # duplicate metadata
        self.metadata.duplicate(nx * ny)
        self.boxdims = np.array([self.boxdims[0] * nx, self.boxdims[1] * ny, self.boxdims[2]])
        self.coords = new_coords

    # -------------------------------------------------------------------------
    # calculating geometric slices
    # -------------------------------------------------------------------------
    def gen_slicepoint(self, buff=20):
        ''' Selects corner of a rectangular bilayer, with a slight buffer'''
        return np.min(self.coords[:, 0:2], axis=0) + [buff, buff]

    def calc_residue_COMS(self):
        com_indices = np.where(self.metadata.ressize > 0)[0]
        sizes  = self.metadata.ressize[com_indices]
        res_coms = np.zeros((sizes.size, 3))  # 0 will be empty
        for i, (ind, size) in enumerate(zip(com_indices, sizes)):
            res_coms[i, :] = np.mean(self.coords[ind:ind + size, :], axis=0)
        return res_coms

    def rectangular_slice(self, xvals, yvals, exclude_radius=0, cutoff_method='com', nlips=None):
        ''' The exclude method is broken right now. Doesn't actually slice bilayer, just returns indices to slice
        '''
        atom_bool_to_keep = np.zeros(self.coords.shape[0], dtype=bool)
        if cutoff_method == 'com':
            res_starts = np.where(self.metadata.ressize > 0)[0]
            res_coms = self.calc_residue_COMS()
            res_bool_to_keep = ((res_coms[:, 0] > xvals[0]) & (res_coms[:, 0] < xvals[1]) &
                                (res_coms[:, 1] > yvals[0]) & (res_coms[:, 1] < yvals[1]))
            if exclude_radius > 0:
                centered_coords = res_coms - res_coms[res_bool_to_keep].mean(axis=0)
                _, rho, _ = cart2pol(centered_coords)
                res_bool_to_keep = (res_bool_to_keep) & (rho > exclude_radius)
            if cutoff_method == 'com':
                for resind, atomind in enumerate(res_starts):
                    atom_bool_to_keep[res_starts[resind]:res_starts[resind] +
                                      self.metadata.ressize[atomind]] = res_bool_to_keep[resind]
        return atom_bool_to_keep

    def circular_slice(self, center, radius, exclude_radius=0, cutoff_method='com'):
        ''' Doesn't actually slice bilayer, just returns indices to slice'''
        atom_bool_to_keep = np.zeros(self.coords.shape[0], dtype=bool)
        res_starts = np.where(self.metadata.ressize > 0)[0]
        res_coms = self.calc_residue_COMS()
        centered_coords = res_coms - [center[0], center[1], 0]
        (_, rho, _) = cart2pol(centered_coords)

        res_bool_to_keep = (rho <= radius) & (rho >= exclude_radius)

        if cutoff_method == 'com':
            for resind, atomind in enumerate(res_starts):
                atom_bool_to_keep[res_starts[resind]:res_starts[resind] +
                                  self.metadata.ressize[atomind]] = res_bool_to_keep[resind]
        return atom_bool_to_keep

    # -------------------------------------------------------------------------
    # file i/o
    # -------------------------------------------------------------------------
    def read_input(self, filein, fmt='pdb', reorganize=False, ignore=[]):
        '''Read input pdb or gro. .gro files are free format, which is tricky for reading if fields are missing or
           don't have whitespace (ie as occurs in manual format between resnumber and resname).
        '''
        with open(filein, "r") as fid:

            if filein[-4:] == '.gro':   # default is pdb
                # remember to convert to angstroms
                stringin = fid.readlines()
                stringin.pop(0)      # remove title
                stringin.pop(0)      # remove atom count, can't use with exclusions

                temp = stringin.pop()  # get tail end for box dims
                while temp.isspace():
                    temp = stringin.pop()    # take care of any potential trailing whitespace
                self.boxdims = [10 * float(i) for i in temp.split()[0:3]]   # final line is box coordinates

                if ignore:
                    string_processed = [line for line in stringin if not line[11:15].strip() in ignore]
                else:
                    string_processed = stringin
                n_atoms = len(string_processed)

                self.coords = np.empty((n_atoms, 3))

                atomname = np.empty(n_atoms, dtype="<U4")
                resname  = np.empty(n_atoms, dtype="<U4")
                curr_res , prev_res, ressize  = [], [], []                   # residue indexing

                atomcount = 0
                for i, line in enumerate(string_processed):
                    atomname[i] = line[11:15]
                    resname[i] =  line[5:9]
                    self.coords[i, :] = [10 * float(j) for j in line[20:43].split()]
                    if not prev_res:  # first iteration
                        prev_res = line[0:5]
                    curr_res = line[0:5]
                    if curr_res == prev_res:
                        atomcount += 1
                    else:
                        ressize.append(atomcount)
                        ressize += [0] * (atomcount - 1)
                        atomcount = 1
                        prev_res = curr_res
                ressize.append(atomcount)  # one more for final residue
                ressize += [0] * (atomcount - 1)
                self.metadata = Metadata(atomname=atomname,
                                         resname=resname,
                                         leaflets=np.zeros(len(ressize), dtype=int),
                                         ressize=np.array(ressize, dtype=int))

            else:
                xcoord, ycoord, zcoord        = [], [], []                   # 3D coordinates
                atomname, resname             = [], []                       # invariant labels
                curr_res , prev_res, ressize  = [], [], []                   # residue indexing
                atomcount = 0                                                # counters
                for pdb_line in fid:
                    if pdb_line.startswith("ATOM") and not pdb_line[17:21].strip() in ignore:
                        # strings
                        atomname.append(pdb_line[12:16])
                        resname.append( pdb_line[17:21])

                        # counting residues
                        if not prev_res:  # first iteration
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

                ressize.append(atomcount)  # one more for final residue
                ressize += [0] * (atomcount - 1)
                self.metadata = Metadata(atomname=np.array(atomname, dtype="<U4"),
                                         resname=np.array(resname, dtype="<U4"),
                                         leaflets=np.zeros(len(ressize), dtype=int),
                                         ressize=np.array(ressize, dtype=int))
                self.coords = np.array((xcoord, ycoord, zcoord)).T

        # assigning resid lengths and leaflets
        self.assign_leaflets()

    def write_coordinates(self, outfile, position='positive', reorder=True, buff=8192, header=None, dummy_name=None,
                          chunksize=100000):

        if reorder:
            self.reorder_by_leaflet()
            self.reorder_within_leaflet()
            if dummy_name:
                self.reorder_with_dummies_in_back(dummy_name)

        # default is to bring everything to positive regime, allows for up to
        # 9999 angstroms in PDB file format
        if position == 'positive':
            self.coords -= np.min(self.coords, axis=0)
        elif position == 'positive_xy':
            self.coords[:, 0:2] = self.coords[:, 0:2] - np.min(self.coords[:, 0:2], axis=0)
        elif position == 'center':
            self.coords -= np.mean(self.coords, axis=0)
        elif position == 'center_xy':
            self.coords[:, 0:2] = self.coords[:, 0:2] - np.mean(self.coords[:, 0:2], axis=0)

        # generating resid list
        nparts = self.coords.shape[0]
        resid = np.zeros(nparts, dtype=int)
        count = 1
        ressize = self.metadata.ressize
        for i in range(nparts):
            size = ressize[i]
            if size > 0:
                resid[i:i + size] = count
                count += 1

        with open(outfile, 'w', buff) as fout:
            chunks = list(range(0, nparts, chunksize))
            chunks.append(nparts + 1)   # so slicing is from 2ndtolast:size

            if outfile[-4:] == '.gro':  # defaults to pdb otherwise
                self.coords /= 10    # internal is angstroms, need to get back to nm
                if header:
                    fout.write('Generated using command: {:s}\n'.format(header))
                else:
                    fout.write('created using BUMPY\n')
                fout.write(' {:d}\n'.format(nparts))
                resid = np.mod(resid, 100000)    # 99,999 max for gro
                atomno = np.mod(np.arange(1, nparts + 1), 100000)
                for startind, stopind in zip(chunks[:-1], chunks[1:]):
                    fout.writelines(["{:5d}{:>5s}{:5s}{:5d}{:8.3f}{:8.3f}{:8.3f}\n".format(
                                    i[0], i[1], i[2], i[3], i[4], i[5], i[6]) for i in zip(
                                    resid[startind:stopind],    # noqa
                                    self.metadata.resname[startind:stopind],
                                    self.metadata.atomname[startind:stopind],
                                    atomno[startind:stopind],
                                    self.coords[startind:stopind, 0],
                                    self.coords[startind:stopind, 1],
                                    self.coords[startind:stopind, 2])])
                fout.write(' {:9.5f} {:9.5f} {:9.5f}\n'.format(self.boxdims[0] / 10, self.boxdims[1] / 10,
                                                               self.boxdims[2] / 10))
            else:
                if header:
                    fout.write('REMARK  - Generated using command: {:s}\n'.format(header))
                fout.write('CRYST1{0:9.3f}{1:9.3f}{2:9.3f}{3:7.2f}{3:7.2f}{3:7.2f}\n'.format(*self.boxdims, 90))

                # calculating residue and atom numbers
                resid = np.mod(resid, 10000)   # 9,999 max for pdb
                atomno = np.mod(np.arange(1, nparts + 1), 100000)

                for startind, stopind in zip(chunks[:-1], chunks[1:]):
                    fout.writelines(["ATOM  {:5d} {:4s} {:4s}{:5d}    {:8.3f}{:8.3f}{:8.3f}\n".format(
                                    i[0], i[1], i[2], i[3], i[4], i[5], i[6]) for i in zip(
                                    atomno[startind:stopind],  # noqa
                                    self.metadata.atomname[startind:stopind],
                                    self.metadata.resname[startind:stopind],
                                    resid[startind:stopind],
                                    self.coords[startind:stopind, 0],
                                    self.coords[startind:stopind, 1],
                                    self.coords[startind:stopind, 2])])

    def write_topology(self, outfile):
        '''Writes out simple topology file (.top)'''
        with open(outfile, 'w') as fout:
            fout.write("\n\n\n[ system ]\nBUMPy system\n\n[ molecules ]\n")
            reslist = self.metadata.resname[np.where(self.metadata.ressize > 0)[0]]
            prev_res = reslist[0]
            counter = 0
            for res in reslist:
                if res == prev_res:
                    counter += 1
                else:
                    fout.write("{:4s} {:d}\n".format(prev_res, counter))
                    counter = 1
                    prev_res = res
            fout.write("{:4s} {:d}\n".format(prev_res, counter))

    def write_index(self, outfile, dummy_name=''):
        ''' Writes out index file (.ndx) with the following (hopefully useful) fields:
                    -system
                    -top_leaflet
                    -top_leaflet_component_1
                    -top_leaflet_component_2...
                    -bot_leaflet
                    -bot_leaflet_componenet_1...
            If dummy particles are present, will NOT include in "top_leaflet" or "bot_leaflet", but will add the
            following sections
                    -not_dummy
                    -dummy
                    -top_dummy
                    -bot_dummy

        '''
        def write_index_unit(print_obj, name, indices):
            print_obj.write("[ {:s} ]\n".format(name))
            count = 1
            for i in indices:
                print_obj.write("{:d} ".format(i))
                count += 1
                if count > 15:  # 15 numbers per line
                    count = 1
                    print_obj.write("\n")
            print_obj.write("\n\n")

        with open(outfile, 'w') as fout:
            top_atomno = np.where((self.metadata.leaflets == 1) & (self.metadata.resname != dummy_name))[0] + 1
            bot_atomno = np.where((self.metadata.leaflets == 0) & (self.metadata.resname != dummy_name))[0] + 1
            write_index_unit(fout, "system", np.arange(1, self.coords.shape[0] + 1))
            write_index_unit(fout, "top_leaflet", top_atomno)
            write_index_unit(fout, "bot_leaflet", bot_atomno)

            # individual lipids by leaflet
            for rname in np.unique(self.metadata.resname):
                if rname != dummy_name:
                    top = np.where((self.metadata.leaflets == 1) & (self.metadata.resname == rname))[0] + 1
                    bot = np.where((self.metadata.leaflets == 0) & (self.metadata.resname == rname))[0] + 1
                    if top.size > 0:
                        write_index_unit(fout, "top_" + rname, top)
                    if bot.size > 0:
                        write_index_unit(fout, "bot_" + rname, bot)

            # dummy info
            if dummy_name:
                write_index_unit(fout, "not_" + dummy_name, np.where(self.metadata.resname != dummy_name)[0] + 1)
                write_index_unit(fout, dummy_name,          np.where(self.metadata.resname == dummy_name)[0] + 1)
                top = np.where((self.metadata.resname == dummy_name) & (self.metadata.leaflets == 1))[0] + 1
                bot = np.where((self.metadata.resname == dummy_name) & (self.metadata.leaflets == 0))[0] + 1
                write_index_unit(fout, "top_" + dummy_name, top)
                write_index_unit(fout, "bot_" + dummy_name, bot)


# ------------------------------------------------------------------------------
# SHAPE REPOSITORY
# ------------------------------------------------------------------------------
class shapes:
    ''' This is the repository for all of the shapes that can be built. The 3 basic shapes are semispheres, cylinders
        and partial tori (junctions), as well as flat bilayers which I don't have a class for, we just take slices out
        of the template. Every shape more complex than that should be a combination of the basic shapes combined
        with some translations and rotations.

        Each shape needs to have 3 static methods:
            1. dimension_requirements - returns an xy dimension which is the minimum size the flat template can be.
            2. final_dimensions       - returns an xyz array of box dimensions for the final shape. This is important
                                        mostly for shapes that have at least 1 periodic connection that need exact
                                        boundaries. It's not so important for closed shapes, and anyone making them will
                                        likely want to change the boundaries of those to establish a certain buffer
                                        prior to solvation
            3. gen_shape              - returns the desired shape in the form of a Molecules instance.

            For the first 2 functions, the geometric arguments are passed in as inputs. For the last one, the correct
            and necessary order is template_bilayer, zo, then geometric arguments
    '''

    class shape:
        '''
        Base class for all shapes to build on. Must have all 3 methods
        '''
        @staticmethod
        def gen_shape(template_bilayer, zo, **geometric_args):
            raise NotImplementedError

        @staticmethod
        def dimension_requirements(**geometric_args):
            raise NotImplementedError

        @staticmethod
        def final_dimensions(**geometric_args):
            raise NotImplementedError

    class flat_bilayer(shape):
        ''' Three uses for a flat bilayer class. First, to make things simpler when building different shapes that
            encompass a flat bilayer. Second, if one wants to make a flat system of a certain size with a dummy grid.
            A third use would simply be expanding an existing bilayer to specified dimensions. One can do something like
            that with gmx editconf, but only in multiples of the original box size
        '''
        @staticmethod
        def dimension_requirements(x_dimension, y_dimension, buff=50):
            return np.array([x_dimension + buff  , y_dimension + buff])

        @staticmethod
        def final_dimensions(x_dimension, y_dimension, buff=50):
            return np.array([x_dimension, y_dimension, 2 * buff])

        @ staticmethod
        def gen_shape(template_bilayer, zo, x_dimension, y_dimension, r_hole=0, cutoff_method='com',
                      print_intermediates=False):
            slice_origin = template_bilayer.gen_slicepoint()
            flat_slice  = template_bilayer.slice_pdb(template_bilayer.rectangular_slice(
                                                     [slice_origin[0], slice_origin[0] + x_dimension],
                                                     [slice_origin[1], slice_origin[1] + y_dimension],
                                                     r_hole))
            flat_slice.center_on_zero()

            if print_intermediates:
                flat_slice.write_coordinates(print_intermediates, position=False)
            return flat_slice

    class semisphere(shape):
        @staticmethod
        def dimension_requirements(r_sphere, buff=50):
            return np.array([(r_sphere + buff) * np.pi  , (r_sphere + buff) * np.pi ])

        @staticmethod
        def final_dimensions(r_sphere, buff=50):
            return np.array([2 * (r_sphere + buff), 2 * (r_sphere + buff), r_sphere + (2 * buff)])

        @ staticmethod
        def gen_shape(template_bilayer, zo, r_sphere, r_hole=0, cutoff_method='com', print_intermediates=False):
            ''' returns molecules instance of semisphere'''
            # calculating slice radii
            slice_radius = np.pi * r_sphere / 2
            top_slice_radius = np.sqrt(2) * (r_sphere + zo[0])
            bot_slice_radius = np.sqrt(2) * (r_sphere - zo[1])
            slice_origin = np.mean(template_bilayer.coords[:, 0:2], axis=0)
            # calculate slice indices
            bool_in_top_slice = template_bilayer.circular_slice(slice_origin, top_slice_radius, exclude_radius=r_hole)
            bool_in_bot_slice = template_bilayer.circular_slice(slice_origin, bot_slice_radius, exclude_radius=r_hole)

            # make slices
            top_leaflet = template_bilayer.slice_pdb(bool_in_top_slice &  template_bilayer.metadata.leaflets)
            bot_leaflet = template_bilayer.slice_pdb(bool_in_bot_slice &  np.invert(template_bilayer.metadata.leaflets))

            # scale slices to slice_radius
            top_leaflet.scale_coordinates_radial(slice_radius / top_slice_radius)
            bot_leaflet.scale_coordinates_radial(slice_radius / bot_slice_radius)
            # merge and transform slices
            top_leaflet.append(bot_leaflet)

            if print_intermediates:
                top_leaflet.write_coordinates(print_intermediates, position=None)

            top_leaflet.spherical_transform(r_sphere)
            return top_leaflet

    class cylinder(shape):

        @staticmethod
        def dimension_requirements(r_cylinder, l_cylinder, completeness=1, buff=50):
            return np.array([ l_cylinder + buff , (r_cylinder + buff) * 2 * np.pi * completeness ])

        @staticmethod
        def final_dimensions(r_cylinder, l_cylinder, buff=50):
            return np.array([l_cylinder, 2 * (r_cylinder + buff), 2 * (r_cylinder + buff)])

        def gen_shape(template_bilayer, zo, r_cylinder, l_cylinder, completeness=1,
                      cutoff_method='com', print_intermediates=False):
            '''Makes a cylinder with given parameters.

            Completeness=0.5 for semicylinder,
            completeness=0.25 for a lateral junction
            '''
            # calculate slice lengths
            cylinder_slice_length = 2 * np.pi * r_cylinder * completeness
            slice_origin = np.min(template_bilayer.coords, axis=0)[0:2] + 40
            outer_slice_length = 2 * np.pi * (r_cylinder + zo[0]) * completeness
            inner_slice_length = 2 * np.pi * (r_cylinder - zo[1]) * completeness
            # calculate slice indices
            xvals = [slice_origin[0], slice_origin[0] + l_cylinder]
            yvals_outer = [slice_origin[1], slice_origin[1] + outer_slice_length]
            yvals_inner = [slice_origin[1], slice_origin[1] + inner_slice_length]
            bool_in_top_slice =  template_bilayer.rectangular_slice(xvals, yvals_outer)
            bool_in_bot_slice =  template_bilayer.rectangular_slice(xvals, yvals_inner)

            # make slices
            top_leaflet = template_bilayer.slice_pdb(bool_in_top_slice &     template_bilayer.metadata.leaflets)
            bot_leaflet = template_bilayer.slice_pdb(bool_in_bot_slice &  np.invert(template_bilayer.metadata.leaflets))
            # scale coordinates
            top_leaflet.scale_coordinates_rectangular([1, cylinder_slice_length / outer_slice_length])
            bot_leaflet.scale_coordinates_rectangular([1, cylinder_slice_length / inner_slice_length])
            top_leaflet.append(bot_leaflet)
            if print_intermediates:
                top_leaflet.write_coordinates(print_intermediates, position=None)

            top_leaflet.cylindrical_transform(r_cylinder)
            return top_leaflet

    class partial_torus(shape):
        @staticmethod
        def dimension_requirements(r_torus, r_tube, buff=50):
            return np.array([2 * (buff + r_torus + r_tube * np.pi) ] * 2)

        @staticmethod
        def final_dimensions(r_torus, r_tube, buff=50):
            return np.array([2 * (buff + r_torus + r_tube * np.pi) ] * 2 + [r_tube + buff])

        @staticmethod
        def gen_shape(template_bilayer, zo, r_torus, r_tube, partial='full',
                      cutoff_method='com', print_intermediates=False):
            '''Makes a partial torus with given parameters

            The flat circular slice of a half_torus with torus R ranges from(R-r') to (R+r'),
            where r' = circumference of tube / 4

            Can get a quarter torus as well using partial='inner' for the shorter junction,
            partial ='outer' for the larger

            partial tori are oriented so that the cylindrical edges point DOWN, and has a minimum at 0 z

            I can't figure out good math on where to center / scale a quarter torus, in regards to the the effect of
            slicing on lipid ratios. Ie, can't just match length of slice to tube radius, as >>WHERE<< you slice changes
            total ratios, as opposed to cylinders (and spheres sorta). SO, will just do an additional slice of the half
            torus if quarter is selected
            '''

            tube_circumference = 2 *  np.pi * r_tube
            inner_tube_circumference = 2 * np.pi * (r_tube - zo[1])
            outer_tube_circumference = 2 * np.pi * (r_tube + zo[0])

            slice_min = r_torus - (tube_circumference / 4)
            slice_max = r_torus + (tube_circumference / 4)
            inner_slice_min = r_torus - (inner_tube_circumference / 4)
            outer_slice_min = r_torus - (outer_tube_circumference / 4)
            inner_slice_max = r_torus + (inner_tube_circumference / 4)
            outer_slice_max = r_torus + (outer_tube_circumference / 4)

            slice_origin = np.mean(template_bilayer.coords, axis=0)[0:2]
            # calculate slice indices
            bool_in_top_slice = template_bilayer.circular_slice(slice_origin, outer_slice_max,
                                                                exclude_radius=outer_slice_min)
            bool_in_bot_slice = template_bilayer.circular_slice(slice_origin, inner_slice_max,
                                                                exclude_radius=inner_slice_min)
            # difference from sphere, exclude center
            top_leaflet = template_bilayer.slice_pdb(bool_in_top_slice &     template_bilayer.metadata.leaflets)
            bot_leaflet = template_bilayer.slice_pdb(bool_in_bot_slice &  np.invert(template_bilayer.metadata.leaflets))

            # scale slices to slice_radius
            top_leaflet.scale_coordinates_toroidal([outer_slice_min, outer_slice_max], [slice_min, slice_max])
            bot_leaflet.scale_coordinates_toroidal([inner_slice_min, inner_slice_max], [slice_min, slice_max])

            top_leaflet.append(bot_leaflet)
            # top_leaflet.write_coordinates('torus_merge_notransform.pdb',position=False)
            # check quarter torus, use circular slice to cut off one side or other
            # the cutoff is r_torus. For inner, just take a circle that ends at r_torus
            # for outer, take circle larger than size of torus including everything,
            # then exclude up to r_torus

            if print_intermediates:
                if partial == 'inner':
                    top_leaflet.write_coordinates(print_intermediates)

            # top_leaflet.write_coordinates('torus_transform.pdb',position=False)
            if partial == 'inner':
                top_leaflet = top_leaflet.slice_pdb(top_leaflet.circular_slice(
                                                    np.mean(top_leaflet.coords, axis=0), r_torus))
            elif partial == 'outer':
                top_leaflet = top_leaflet.slice_pdb(top_leaflet.circular_slice(np.mean(top_leaflet.coords, axis=0),
                                                    r_torus + tube_circumference, exclude_radius=r_torus))

            if print_intermediates:
                top_leaflet.write_coordinates(print_intermediates, position=None)
            top_leaflet.toroidal_transform(r_torus, r_tube)

            return top_leaflet

    class sphere(shape):

        @staticmethod
        def dimension_requirements(r_sphere, buff=50):
            return shapes.semisphere.dimension_requirements(r_sphere)

        @staticmethod
        def final_dimensions(r_sphere, buff=50):
            return np.array([2 * (r_sphere + buff)] * 3)

        @staticmethod
        def gen_shape(template_bilayer, zo, r_sphere, n_holes=0):
            top_half = shapes.semisphere.gen_shape(template_bilayer, zo, r_sphere, False)
            bot_half = deepcopy(top_half)
            bot_half.rotate([180, 0, 0])
            top_half.append(bot_half, preserve_leaflets=True)
            return top_half

    class torus(shape):
        @staticmethod
        def dimension_requirements(r_torus, r_tube, buff=50):
            return shapes.partial_torus.dimension_requirements(r_torus, r_tube)

        @staticmethod
        def final_dimensions(r_torus, r_tube, buff=50):
            return np.array([2 * (buff + r_torus + r_tube * np.pi) ] * 2 + [2 * (r_tube + buff)])

        @staticmethod
        def gen_shape(template_bilayer, zo, r_torus, r_tube, completeness=0.5):

            if r_tube >= r_torus:
                raise UserWarning("r_torus should be less than r_tube for a ring torus")

            top_half = shapes.partial_torus.gen_shape(template_bilayer, zo, r_torus, r_tube, partial='full')
            bot_half = deepcopy(top_half)
            bot_half.rotate([180, 0, 0])
            top_half.append(bot_half)
            return top_half

    class semicylinder_plane(shape):

        @staticmethod
        def dimension_requirements(r_cylinder, l_cylinder, r_junction, l_flat, buff=50):
            cyldims = shapes.cylinder.dimension_requirements(r_cylinder, l_cylinder, completeness=0.5)
            flatdims = np.array([l_cylinder, l_flat])
            jdims  =  shapes.cylinder.dimension_requirements(r_junction, l_cylinder, completeness=0.25)
            return np.array([max([cyldims[0], flatdims[0], jdims[0]]), max([cyldims[1], flatdims[1], jdims[1]])])

        @staticmethod
        def final_dimensions(r_cylinder, l_cylinder, r_junction, l_flat, buff=50):
            return np.array([l_cylinder, 2 * (r_cylinder + r_junction) + l_flat, r_cylinder + r_junction + 2 * buff])

        @staticmethod
        def gen_shape(template_bilayer, zo, r_cylinder, l_cylinder, r_junction, l_flat):
            semicyl   = shapes.cylinder.gen_shape(template_bilayer, zo, r_cylinder, l_cylinder, completeness=0.5)

            template_2 = deepcopy(template_bilayer)
            template_2.rotate([180, 0, 0])     # junctions show inner leaflet to top, so reverse coordinates
            template_2.metadata.leaflets = 1 - template_2.metadata.leaflets   # fix leaflet description

            junction  = shapes.cylinder.gen_shape(template_2, zo, r_junction, l_cylinder, completeness=0.25)
            junction.metadata.leaflets = 1 - junction.metadata.leaflets   # now reverse leaflets again to match "top"
            junction2 = deepcopy(junction)

            # rotations and translations. 135 and 225 degrees gets junctions rotated so that they max at 0 and taper to
            # flat in the y direction translation because they face the wrong direction and may not match cylinder
            # directions anyway
            junction.rotate([135, 0, 0])
            junction.translate([0, -(r_junction + r_cylinder), 0])
            junction2.rotate([225, 0, 0])
            junction2.translate([0, r_junction + r_cylinder, 0])

            flat_slice = shapes.flat_bilayer.gen_shape(template_bilayer, zo, l_cylinder, l_flat)
            flat_slice.translate( [0, r_junction + r_cylinder + (l_flat / 2), - r_junction])

            semicyl.append(junction)
            semicyl.append(junction2)
            semicyl.append(flat_slice)
            return semicyl

    class semisphere_plane(shape):

        @staticmethod
        def dimension_requirements(r_sphere, r_junction, l_flat, buff=50):
            sphdims = shapes.semisphere.dimension_requirements(r_sphere)
            flatdims = np.array([l_flat, l_flat])
            jdims  =  shapes.partial_torus.dimension_requirements(r_sphere + r_junction, r_junction)
            return np.array([max([sphdims[0], flatdims[0], jdims[0]]), max([sphdims[1], flatdims[1], jdims[1]])])

        @staticmethod
        def final_dimensions(r_sphere, r_junction, l_flat, buff=50):
            return np.array([l_flat, l_flat, r_sphere + r_junction + 2 * buff])

        @staticmethod
        def gen_shape(template_bilayer, zo, r_sphere, r_junction, l_flat):
            semisph   = shapes.semisphere.gen_shape(template_bilayer, zo, r_sphere)

            template_2 = deepcopy(template_bilayer)
            template_2.rotate([180, 0, 0])     # junctions show inner leaflet to top, so reverse coordinates
            template_2.metadata.leaflets = 1 - template_2.metadata.leaflets   # fix leaflet description

            junction  = shapes.partial_torus.gen_shape(template_2, zo,
                                                       r_sphere + r_junction, r_junction, partial='inner')
            junction.metadata.leaflets = 1 - junction.metadata.leaflets   # now reverse leaflets again to match "top"

            junction.rotate([180, 0, 0])

            flat_slice = shapes.flat_bilayer.gen_shape(template_bilayer, zo, l_flat, l_flat, r_sphere + r_junction)
            flat_slice.translate([0, 0, -r_junction])

            semisph.append(junction)
            semisph.append(flat_slice)
            return semisph

    class double_bilayer_cylinder(shape):
        @staticmethod
        def dimension_requirements(r_cylinder, l_cylinder, r_junction, l_flat, buff=50):
            cyldims = shapes.cylinder.dimension_requirements(r_cylinder, l_cylinder)
            flatdims = np.array([l_flat] * 2)
            jdims  =  shapes.cylinder.dimension_requirements(r_junction, l_cylinder, completeness=0.5)
            return np.array([max([cyldims[0], flatdims[0], jdims[0]]), max([cyldims[1], flatdims[1], jdims[1]])])

        @staticmethod
        def final_dimensions(r_cylinder, l_cylinder, r_junction, l_flat, buff=50):
            return np.array([l_flat, l_flat, l_cylinder + 2 * (r_junction + buff)])

        def gen_shape(template_bilayer, zo, r_cylinder, l_cylinder, r_junction, l_flat):
            ''' "Top" leaflet will be IMS facing leaflet, so: TOP of flat region (no change), INSIDE of cylinder (flip),
                and OUTSIDE of torus (no change)
            '''
            if (r_cylinder + r_junction) > (l_flat / 2):
                raise UserWarning("Flat region too small for cylinder/junction radii")

            junction = shapes.partial_torus.gen_shape(template_bilayer, zo, r_cylinder + r_junction, r_junction,
                                                      partial='inner')
            junction.translate([0, 0, l_cylinder / 2])
            junction_2 = deepcopy(junction)
            junction_2.rotate( [180, 0, 0])
            flat_bilayer = shapes.flat_bilayer.gen_shape(template_bilayer, zo, l_flat, l_flat, r_cylinder + r_junction)

            flat_bilayer.translate([0, 0, (l_cylinder / 2) + r_junction])
            flat_bilayer_2 = deepcopy(flat_bilayer)
            flat_bilayer_2.rotate([180, 0, 0])

            template_2 = deepcopy(template_bilayer)
            template_2.rotate([180, 0, 0])    # flip
            template_2.metadata.leaflets = 1 - template_2.metadata.leaflets
            cyl = shapes.cylinder.gen_shape(template_2, zo, r_cylinder, l_cylinder, completeness=1)
            cyl.metadata.leaflets = 1 - cyl.metadata.leaflets
            cyl.rotate([0, 90, 0])

            cyl.append(junction)
            cyl.append(flat_bilayer)
            cyl.append(junction_2)
            cyl.append(flat_bilayer_2)
            return cyl

    class capped_cylinder(shape):
        @staticmethod
        def dimension_requirements(r_cylinder, l_cylinder, buff=50):
            cyldims = shapes.cylinder.dimension_requirements(r_cylinder, l_cylinder)
            sphdims = shapes.semisphere.dimension_requirements(r_cylinder)
            return np.array([max([cyldims[0], sphdims[0]]), max([cyldims[1], sphdims[1]])])

        @staticmethod
        def final_dimensions(r_cylinder, l_cylinder, buff=50):
            return np.array([l_cylinder + 2 * (r_cylinder + buff), 2 * (r_cylinder + buff), 2 * r_cylinder + buff ])

        @staticmethod
        def gen_shape(template_bilayer, zo, r_cylinder, l_cylinder):
            ''' Two semispheres connected by a cylinder'''
            cyl = shapes.cylinder.gen_shape(template_bilayer, zo, r_cylinder, l_cylinder, completeness=1)
            semisphere1 = shapes.semisphere.gen_shape(template_bilayer, zo, r_cylinder)
            semisphere2 = deepcopy(semisphere1)
            semisphere1.rotate([0, 90, 0])
            semisphere2.rotate([0, 270, 0])
            semisphere1.coords[:, 0] = semisphere1.coords[:, 0] - l_cylinder / 2
            semisphere2.coords[:, 0] = semisphere2.coords[:, 0] + l_cylinder / 2
            cyl.append(semisphere1)
            cyl.append(semisphere2)
            return cyl

    class sphere_cylinder(shape):
        @staticmethod
        def dimension_requirements(r_sphere, r_cylinder, l_cylinder, r_junction):
            cyldims = shapes.cylinder.dimension_requirements(r_cylinder, l_cylinder, completeness=0.5)
            sphdims = shapes.semisphere.dimension_requirements(r_sphere)
            jdims  =  shapes.cylinder.dimension_requirements(r_junction, l_cylinder, completeness=0.5)
            return np.array([max([cyldims[0], sphdims[0], jdims[0]]), max([cyldims[1], sphdims[1], jdims[1]])])

        @staticmethod
        def final_dimensions(r_sphere, r_cylinder, l_cylinder, r_junction, buff=50):
            yzdim = 2 * (r_sphere + buff)
            xdim  = l_cylinder + (2 * np.sqrt(r_sphere ** 2 - (r_cylinder + r_junction) ** 2)) + (2 * r_junction)
            return np.array([xdim, yzdim, yzdim])

        @staticmethod
        def gen_shape(template_bilayer, zo, r_sphere, r_cylinder, l_cylinder, r_junction):
            cyl = shapes.cylinder.gen_shape(template_bilayer, zo, r_cylinder, l_cylinder, completeness=1)
            sph1 = shapes.semisphere.gen_shape(template_bilayer, zo, r_sphere, r_hole=r_cylinder + r_junction)
            sph2 = deepcopy(sph1)
            sph1.rotate([0,  90, 0])
            sph2.rotate([0, 270, 0])
            junc1 = shapes.partial_torus.gen_shape(template_bilayer, zo, r_cylinder + r_junction,
                                                   r_junction, partial='inner')
            junc2 = deepcopy(junc1)
            junc1.rotate([0,  90, 0])
            junc2.rotate([0, 270, 0])

            # the lateral distance from the center of the sphere to the hole is sqrt(Rtotal^2 - Rhole^2) by trig
            d_sphere = np.sqrt(r_sphere ** 2 - (r_cylinder + r_junction) ** 2)

            junc1.coords += [d_sphere + r_junction, 0, 0]
            junc2.coords -= [d_sphere + r_junction, 0, 0]
            cyl.coords += [(l_cylinder / 2) + d_sphere + r_junction , 0, 0 ]
            cyl.append(sph1)
            cyl.append(sph2)
            cyl.append(junc1)
            cyl.append(junc2)
            return cyl


def parse_command_lines():
    ''' Parses command line for parameters, returns parsed arguments '''
    prog_name =  'BUMPY'
    prog_description = 'Creating curved membrane systems with arbitrary geometry and lipid composition, using a ' + \
                       'pivotal plane-based approach to appropriately match inter-leaflet area differences'

    geometry_description = 'Geometric arguments should be added as a series of argument:value pairs separated by a ' + \
                           'colon. Run this program with the --list flag for a list of supported shapes and their '  + \
                           'respective geometric arguments.'

    parser = ArgumentParser(prog=prog_name, description=prog_description, add_help=False, allow_abbrev=False,
                            usage='')
    # groups
    required_inputs     = parser.add_argument_group('required inputs')
    geometric_arguments = parser.add_argument_group('geometric arguments', geometry_description)
    optional_arguments  = parser.add_argument_group('optional arguments')
    dummy_arguments     = parser.add_argument_group('dummy particle options')
    output_arguments    = parser.add_argument_group('output arguments')

    # mandatory input
    required_inputs.add_argument('-s', help='Shape to make - see manual for a list of shapes', metavar='')
    required_inputs.add_argument('-f', help='Flat bilayer template to be used as a template',  metavar='')
    required_inputs.add_argument('-z', metavar='', help='Location of the pivotal plane (nm). Just one value, or' +
                                                        'outer_zo:inner_zo')

    # geometry
    geometric_arguments.add_argument('-g', nargs='*', help='Format is arg:value, ie r_cylinder:10 ' +
                                     ' l_cylinder:20  ... for every geometric parameter in shape', metavar='')

    # optional arguments
    optional_arguments.add_argument('-h', '--help', action='help', help='show this help message and exit')
    optional_arguments.add_argument('-l', '--list', default=False, action='store_true',
                                    help='List current repository of shapes and their geometric arguments')
    optional_arguments.add_argument('--outer', default='top', help='By default, top leaflet = outer leaflet. ' +
                                    'Set to "bot" to invert', metavar='')
    optional_arguments.add_argument('--apl', metavar='', help='Slice top bilayer to achieve a specific area per ' +
                                    'lipid in final shape - not yet implemented', default=None)
    optional_arguments.add_argument('--ignore_resnames', metavar='', help='colon separated list of resnames to ignore' +
                                    'when reading in a structure file, for example to exclude water', default=[],
                                    nargs="*")

    dummy_arguments.add_argument('--gen_dummy_particles', action='store_true',
                                 help='Add a grid of dummy particles surrounding bilayer' )
    dummy_arguments.add_argument('--dummy_name', metavar='', type=str, default='DUMY',
                                 help='Dummy particle atomname/resname. Defaults to DUMY')
    dummy_arguments.add_argument('--dummy_grid_thickness', metavar='', type=float,
                                 help='Create dummy array with thickness specified')
    dummy_arguments.add_argument('--dummy_grid_spacing', metavar='', type=float, default=5,
                                 help='dummy grid spacing distance')

    # output files
    output_arguments.add_argument('-o', help='Output structure - only PDBs for now', default='confout.pdb', metavar='')
    output_arguments.add_argument('-p', help='Simple .top topology file', metavar='')                    # optional
    output_arguments.add_argument('-n', help='Simple .ndx index file, separating leaflets', metavar='')  # optional

    return parser.parse_args()


def gen_dummy_grid(lateral_distance=5, thickness=50, atomname='DUMY', resname='DUMY'):
    ''' Creates a  10 * 10 * 2 grid of dummy particles, with two sheet separated by thickness and
        inter-particle distances separated by lateral distance
    '''
    coords = np.array(np.meshgrid(np.arange(10) * lateral_distance, np.arange(10) * lateral_distance,
                                  [ -thickness / 2, thickness / 2])).T.reshape(-1, 3)   # thanks @pv, stackoverflow
    atomname = np.array([atomname] * 200, dtype="<U4")
    resname  = np.array([resname]  * 200, dtype="<U4")
    ressize = np.ones(200, dtype=int)
    leaflets = coords[:, 2] > 0
    meta = Metadata(atomname=atomname, resname=resname, leaflets=leaflets, ressize=ressize)
    return Molecules(metadata=meta, coords=coords,
                     boxdims=np.array([10 * lateral_distance, 10 * lateral_distance, thickness]))


def display_parameters(cl_args):
    '''Displays selected parameters upon command line execution'''
    pass


# -----------------------------------------------------------------------------
# Command line start
# -----------------------------------------------------------------------------
def main():
    cli = " ".join(sys.argv)
    args = parse_command_lines()
    display_parameters(args)  # show user what they selected

    if args.list:
        print('Shapes in repository, with listed arguments:\n')
        for shape in inspect.getmembers(shapes):
            if (not shape[0].startswith('_')) and  (shape[0] != 'shape'):
                print('{:s}'.format(shape[0]))
                sig = inspect.signature(shape[1].gen_shape)
                for param in sig.parameters.values():
                    if param.default == param.empty and param.name != 'zo' and param.name != 'template_bilayer':
                        print('    {:s}'.format(param.name))
        exit()

    # parse arguments
    geometric_args = {garg.split(':')[0] : float(garg.split(':')[1]) for garg in args.g }

    zo = [float(i) for i in args.z.split(':')]
    if len(zo) == 1:    # if one value given, apply to both leaflets
        zo *= 2

    shape_tobuild = getattr(shapes, args.s)

    # read in pdb
    print('Reading in PDB file  ... ', end='', flush=True)
    t = time()
    template_bilayer = Molecules(infile=args.f, ignore=args.ignore_resnames)
    print('Finished reading PDB with {:d} atoms - time elapsed = {:.1f} seconds'.format(
          template_bilayer.coords.shape[0], time() - t))

    # flip bilayer if bottom is selected, scale if selected, multiply laterally to appropriate size
    if args.outer == 'bot':
        template_bilayer.rotate([180, 0, 0], com=True)
        template_bilayer.metadata.leaflets = 1 - template_bilayer.metadata.leaflets
    if args.apl:
        currarea = template_bilayer.boxdims[0] * template_bilayer.boxdims[1]
        newarea = args.apl * template_bilayer.coords.shape[0] / 2  # 2 leaflets
        ratio = np.sqrt(newarea / currarea)
        template_bilayer.scale_coordinates_rectangular(ratio)

    mult_factor = (np.ceil(shape_tobuild.dimension_requirements(**geometric_args) /
                           template_bilayer.boxdims[0:2]).astype(int))
    template_bilayer.duplicate_laterally(*mult_factor)

    # make shape
    print('Generating shape     ... ', end='', flush=True)
    t = time()
    # construct the shape
    shape = shape_tobuild.gen_shape(template_bilayer, zo, **geometric_args)
    shape.boxdims = shape_tobuild.final_dimensions(**geometric_args)
    print('Finished - time elapsed = {:.1f} seconds'.format(time() - t))

    if args.gen_dummy_particles:
        print('Creating dummy particles')
        dummy_name = args.dummy_name[0:4]   # shorten to first 4 letters
        dummy_template = gen_dummy_grid(thickness=args.dummy_grid_thickness, lateral_distance=args.dummy_grid_spacing,
                                        atomname=dummy_name, resname=dummy_name)
        mult_factor = (np.ceil( shape_tobuild.dimension_requirements(**geometric_args) /
                                dummy_template.boxdims[0:2]).astype(int))
        dummy_template.duplicate_laterally(*mult_factor)
        dummy_shape = shape_tobuild.gen_shape(dummy_template, 2 * [args.dummy_grid_thickness / 2], **geometric_args)
        shape.append(dummy_shape)
    else:
        dummy_name = ''

    # file output
    print('Writing out PDB file ... ', end='', flush=True)
    t = time()
    shape.write_coordinates(args.o, header=cli, dummy_name=dummy_name)
    if args.p:
        shape.write_topology(args.p)
    if args.n:
        shape.write_index(args.n, dummy_name=dummy_name)

    print('Finished writing PDB file with {} atoms - time elapsed = {:.1f} seconds'.format(
          shape.coords.shape[0], time() - t))


if __name__ == '__main__':
    main()
