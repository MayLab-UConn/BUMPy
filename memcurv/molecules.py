''' Molecules class has coordinates, pdb i/o stuff, manipulation'''

import numpy as np
import nonrigid_coordinate_transformations as nrb
import rigid_body_transforms as rb
import time

class Molecules:

    def __init__(self,infile=None,string_info=[],atomno=[],resid=[],coords=[],
                 leaflets=[],resid_list=[]):
        '''Can initialize from file, or with manual inputs (like when slicing).
           Can also initialize with all objects blank, and wait for input
        '''
        # read from file
        if infile is not None:
            self.read_pdb(infile)
        # assign manually
        else:
            self.string_info = string_info  # dictionary, nparts size
            self.atomno    = atomno     #  np int,
            self.resid     = resid      #  np int,
            self.coords    = coords     # np float, 3D array of nparts * xyz
            # organizational variables
            self.leaflets  = leaflets   # 1 for top, 0 for bot
            self.resid_list= resid_list # nres size dictionary


    # -------------------------------------------------------------------------
    # Calculate and superficially change dataset properties
    # -------------------------------------------------------------------------

    def get_current_dims(self):
        ''' Returns magnitude of coordinate range'''
        return np.ptp(self.coords,axis=0)

    def get_current_minmax(self):
        ''' Returns minima and maxima in form of 2*3 matrix. Min first row,
            max second row
        '''
        return np.vstack((np.min(self.coords,axis=0),np.max(self.coords,axis=0)))
    def calc_thickness(self):
        pass
    def assign_leaflets(self):
        ''' Labels indices as top (1) or bottom (0) leaflet based on COM of
           entire residue
        '''
        self.leaflets = np.empty_like(self.atomno)
        bilayer_com = np.mean(self.coords[:,2])

        for i in self.resid_list.keys():
            res_indices = self.resid_list[i]
            res_com = np.mean(self.coords[res_indices,2])
            if res_com > bilayer_com:
                self.leaflets[res_indices] = 1
            else:
                self.leaflets[res_indices] = 0

    def renumber_atoms(self):
        ''' Renumber atoms from 1 to natoms'''
        self.atomno = np.arange(1,self.atomno.size + 1)

    def renumber_resids(self):
        ''' Renumber residues from 1 to n_residues '''
        counter = 1                 # what new resid will be
        i = 0                       # indexing variable
        curr_resid = self.resid[i]  # check before overriding
        while i < (self.resid.size):
            self.resid[i]= counter
            if (i+1) < self.resid.size:
                if self.resid[i+1] != curr_resid:
                    counter += 1
                    curr_resid = self.resid[i+1]
            i +=1

    def reorder_by_leaflet(self):
        new_index_order = np.append(np.where(self.leaflets == 1),
                                    np.where(self.leaflets == 0))
        self.coords = self.coords[new_index_order,:]
        self.atomno = self.atomno[new_index_order]
        self.resid =  self.resid[new_index_order]
        self.leaflets = self.leaflets[new_index_order]
        self.string_info['atomtype'] = [self.string_info['atomtype'][i] for i in new_index_order]
        self.string_info['atomname'] = [self.string_info['atomname'][i] for i in new_index_order]
        self.string_info['resname'] =  [self.string_info['resname'][i]  for i in new_index_order]
        self.string_info['chain'] =    [self.string_info['chain'][i]    for i in new_index_order]
        self.string_info['junk'] =     [self.string_info['junk'][i]     for i in new_index_order]

    def assign_resid_list(self,wipe=True):
        if wipe:
            self.resid_list = dict()                   # start from scratch
            for i in np.unique(self.resid):
                self.resid_list[i] = []
        for i in range(0,len(self.atomno)):
            self.resid_list[self.resid[i]].append(i)

    def reorganize_components(self):
        '''Renumbers atoms and resids according to leaflet, also recalculates
           organizational arrays
        '''
        self.renumber_atoms()
        self.renumber_resids()
        self.assign_resid_list()
        self.assign_leaflets()
        self.reorder_by_leaflet()
        self.renumber_atoms()
        self.renumber_resids()
        self.assign_resid_list()
    # -------------------------------------------------------------------------
    # adding and slicing pdb classes
    # -------------------------------------------------------------------------
    def append_pdb(self,new_pdb):
        '''appends all information from new_pdb to end of current pdb '''
        # strings
        self.string_info['atomtype'].extend(new_pdb.string_info['atomtype'])
        self.string_info['atomname'].extend(new_pdb.string_info['atomname'])
        self.string_info['resname'].extend(new_pdb.string_info['resname'])
        self.string_info['chain'].extend(new_pdb.string_info['chain'])
        self.string_info['junk'].extend(new_pdb.string_info['junk'])

        # numpy arrays
        self.atomno   = np.append(self.atomno,new_pdb.atomno)
        self.resid    = np.append(self.resid, new_pdb.resid)
        self.coords   = np.vstack((self.coords,new_pdb.coords))
        self.leaflets = np.append(self.leaflets,new_pdb.leaflets)
        # redo organizational arrays
        self.reorganize_components()
    def slice_pdb(self,slice_indices):
        ''' returns a new instance of current pdb class with sliced indices'''
        string_slice = {'atomtype':[self.string_info['atomtype'][i] for i in slice_indices],
                        'atomname':[self.string_info['atomname'][i] for i in slice_indices],
                        'resname' :[self.string_info['resname'][i]  for i in slice_indices],
                        'chain'   :[self.string_info['chain'][i]    for i in slice_indices],
                        'junk'    :[self.string_info['junk'][i]     for i in slice_indices]}
        molecule_slice = Molecules(infile=None,
                                   string_info=string_slice,
                                   atomno=self.atomno[  slice_indices],
                                   resid=self.resid[   slice_indices],
                                   leaflets=self.leaflets[slice_indices],
                                   coords=self.coords[  slice_indices,:])
        molecule_slice.reorganize_components()
        return molecule_slice


    # -------------------------------------------------------------------------
    # calculating geometric slices
    # -------------------------------------------------------------------------
    def gen_random_slice_point(self,slice_r):
        ''' Within boundaries, randomize slicing region'''
        return np.mean(self.coords,axis=0)[0:2]

    def rectangular_slice(self,xvals,yvals,partial_molecule='exclude'):
        '''Slices pdb to include only rectangular segment from x[0] to x[1] and
           y[0] to y[1]. Default is to exclude partial molecules, have option to
           include partial molecules or make whole and include.

           RETURNS INDICES, not actual slice
        '''
        indices_tokeep = np.array([],dtype=int)
        all_inrange = np.asarray(np.where( (self.coords[:,0] > xvals[0]) &
                                           (self.coords[:,0] < xvals[1]) &
                                           (self.coords[:,1] > yvals[0]) &
                                           (self.coords[:,1] < yvals[1]) ))
        if partial_molecule == 'exclude':
            print('excluding partial molecules')
            for i in self.resid_list.keys():
                resvals = np.array(self.resid_list[i])
                keep = True
                for j in resvals:
                    if not j in all_inrange: # every index of residue must be
                        keep = False         # in range
                        break
                if keep:
                    indices_tokeep =np.append(indices_tokeep,np.array(resvals))
        return indices_tokeep

    def circular_slice(self,center,radius,exclude_radius=0,
                       partial_molecule='exclude'):
        indices_tokeep = np.array([],dtype=int)
        centered_coords = self.coords - [center[0],center[1],0]
        (theta,rho,z) = nrb.cart2pol(centered_coords)
        all_inrange = np.asarray(np.where((rho <= radius) & (rho>= exclude_radius)))
        if partial_molecule == 'exclude':
            for i in self.resid_list.keys():
                resvals = np.array(self.resid_list[i])
                keep = True
                for j in resvals:
                    if not j in all_inrange: # every index of residue must be
                        keep = False         # in range
                        break
                if keep:
                    indices_tokeep =np.append(indices_tokeep,np.array(resvals))

        return indices_tokeep

    # -------------------------------------------------------------------------
    # file i/o
    # -------------------------------------------------------------------------

    def read_pdb(self,pdbfile,reorganize=True):
        '''Read input pdb, default is to renumber and reorderatoms and resids
           based on leaflets (top first)
        '''
        print('Reading in PDB file')
        t = time.time()
        fid = open(pdbfile,"r")
        # initialize temporary variables
        self.string_info ={'atomtype':[],'atomname':[],'resname':[],'chain':[],'junk':[]}
        atomnum = []
        resnum  = []
        xcoord  = []
        ycoord  = []
        zcoord  = []
        toadd_resid = 0
        toadd_atomno= 0
        prev_res = []
        for pdb_line in fid:
            if pdb_line.startswith("ATOM") or pdb_line.startswith("HETATM"):
                # strings
                self.string_info['atomtype'].append(pdb_line[ 0: 6])
                self.string_info['atomname'].append(pdb_line[12:16])
                self.string_info['resname'].append( pdb_line[16:21])
                self.string_info['chain'].append(   pdb_line[21:22])
                self.string_info['junk'].append(    pdb_line[54:  ])
                # integers
                atomnum.append(int(pdb_line[ 6:11]) + toadd_atomno)
                resnum.append( int(pdb_line[22:26]) + toadd_resid)
                # account for pdb truncation
                if pdb_line[6:11]   == '99999':
                    toadd_atomno += 100000
                if pdb_line[22:26]  == '9999' and pdb_line[22:26] != prev_res:
                    toadd_resid  += 10000
                # floats for coordinate array
                xcoord.append(float(pdb_line[30:38]))
                ycoord.append(float(pdb_line[38:46]))
                zcoord.append(float(pdb_line[46:54]))
                prev_res = pdb_line[22:26]
        # now turn numbers into numpy
        self.atomno = np.array(atomnum)
        self.resid =  np.array(resnum)
        x = np.array(xcoord) # this is a stupid way to concatenate things,
        y = np.array(ycoord) # but I suck with lists
        z = np.array(zcoord)
        self.coords = np.stack((x,y,z),axis=1)
        fid.close()   # test this?
        print('Finished reading in PDB file with {} atoms,time elapsed = {:.1f} seconds'.format(self.atomno.size,time.time()-t))
        if reorganize:
            t = time.time(); print('Reformatting PDB input')
            self.reorganize_components()
            print('Finished formatting PDB input, time required = {:.1f} seconds'.format(time.time()-t))
    def write_pdb(self,outfile):
        print('Writing out PDB file')
        t = time.time()
        '''Outputs to pdb file, CRYST1 and ATOM lines only'''
        # will need modulus for atomno(max=99,999) AND resno (max = 9,999)
        fout = open(outfile,'w')
        nparts = len(self.atomno)
        dims = self.get_current_dims()
        # write out box dims
        line = 'CRYST1{0:9.3f}{1:9.3f}{2:9.3f}{3:7.2f}{3:7.2f}{3:7.2f}'.format(
                dims[0],dims[1],dims[2],90)
        fout.write(line + '\n')

        for i in range(nparts):
            line =("{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s} "
                  + "{:8.3f}{:8.3f}{:8.3f}{:s}")
            fout.write(line.format(self.string_info['atomtype'][i],
                                   np.mod(self.atomno[i],100000), # cap at 10^5
                                   self.string_info['atomname'][i],
                                   " ", self.string_info['resname'][i],
                                   self.string_info['chain'][i],
                                   np.mod(self.resid[i],10000),   # cap at 10^4
                                   " ", self.coords[i,0],
                                   self.coords[i,1],
                                   self.coords[i,2],
                                   self.string_info['junk'][i]))
        fout.close()
        print('Finished writing PDB file with {} atoms,time elapsed = {:.1f} seconds'.format(self.atomno.size,time.time()-t))

    def write_topology(self,outfile):
        '''Writes out simple topology file (.top)'''
        pass
    def write_index(self,outfile):
        '''Writes out index file (.ndx) with the following (hopefully useful)
           fields:
                    -top_leaflet
                    -top_leaflet_component_1
                    -top_leaflet_component_2...
                    -bot_leaflet
                    -bot_leaflet_componenet_1...
        '''
        pass
