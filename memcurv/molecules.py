''' Molecules class has coordinates, pdb i/o stuff, manipulation'''

import numpy as np
import pandas as pd
import nonrigid_coordinate_transformations as nrb
import rigid_body_transforms as rb
import time

class Molecules:

    def __init__(self,infile=None,metadata=[],coords=[]):
        '''Can initialize from file, or with manual inputs (like when slicing).
           Can also initialize with all objects blank, and wait for input
        '''
        # read from file
        if infile is not None:
            self.read_pdb(infile)
        # assign manually
        else:
            self.coords      = coords     # np float, 3D array of nparts * xyz
            # organizational variables
            self.metadata    = metadata

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

    def assign_leaflets(self):
        ''' Labels indices as top (1) or bottom (0) leaflet based on COM of
           entire residue
        '''

        ####coms = self.calc_residue_COMS()
        res_starts = np.where(self.metadata.resid_length > 0)[0]
        coms = self.coords[res_starts,2]
        bilayer_com = np.mean(self.coords[:,2])
        leaflets = np.zeros(self.coords.shape[0],dtype=int)
        for i in range(len(self.resid_list)):
            #leaflets[self.resid_list[i]] = int(coms[i,2] > bilayer_com)
            leaflets[self.resid_list[i]] = int(coms[i] > bilayer_com)
        self.metadata.leaflets =leaflets

    def renumber_resids(self):
        ''' Renumber residues from 1 to n_residues '''
        reslist_ind = np.where(self.metadata.resid_length > 0)[0]
        lengths     = self.metadata.resid_length[reslist_ind]
        counts = list(range(1,reslist_ind.size+1))
        new_resids = np.zeros(self.metadata.resid_length.size,dtype=int)
        for ind,l_ind,count in zip(reslist_ind,lengths,counts):
            new_resids[ind:ind+l_ind] = count
        self.metadata.resid = new_resids

    def reorder_by_leaflet(self):
        ''' Switches up order of atoms so that top leaflet comes first,
            bot comes second
        '''
        new_index_order = np.append(np.where(self.metadata.leaflets == 1),
                                    np.where(self.metadata.leaflets == 0))
        self.coords = self.coords[new_index_order,:]
        self.metadata.reindex(new_index_order)

    def assign_resids(self,wipe=True):
        '''Initializes the "resid_length" section of metadata, where the first
           atom of each residue contains the length of that residue.
        '''
        resid_list,resid_start,resid_length = np.unique(self.metadata['resid'],
                                              return_index=True,
                                              return_counts=True)
        self.metadata.loc[resid_start,'resid_length'] = resid_length

    def gen_resid_list(self):
        '''Generates the resid_list list, which for each residue contains the
           list of atoms in that residue. This can be inferred from the "resid_length"
           section of metadata instead, but would have to be recalculated every
           time a residue-based cutoff is made.
        '''
        reslist_ind = np.where(self.metadata.resid_length > 0)[0]
        lengths = self.metadata.resid_length[reslist_ind]
        self.resid_list = [list(range(ind,ind +l_ind)) for ind,l_ind in zip(reslist_ind,lengths)]


    def reorganize_components(self,assign_resids=True,reset_leaflets=False,
                              reorder_by_leaflets=True,renumber_resids=True):
        '''Renumbers atoms and resids according to leaflet, also recalculates
           organizational arrays
        '''
        self.metadata.reset_index(drop=True,inplace=True)

        if assign_resids:
            self.assign_resids()
        self.gen_resid_list()
        if reset_leaflets:
            self.assign_leaflets()   # need resid list for this
        if reorder_by_leaflets:
            self.reorder_by_leaflet()
        if renumber_resids:
            self.renumber_resids()
    # -------------------------------------------------------------------------
    # adding and slicing pdb classes
    # -------------------------------------------------------------------------
    def append_pdb(self,new_pdb,preserve_leaflets=False):
        '''appends all information from new_pdb to end of current pdb '''
        # strings
        self.metadata.append(new_pdb.metadata)
        self.coords   = np.vstack((self.coords,new_pdb.coords))
        # redo organizational arrays
        self.reorganize_components(preserve_leaflets)
    def slice_pdb(self,slice_indices):
        ''' returns a new instance of current pdb class with sliced indices'''
        metadata_slice = self.metadata.ix[slice_indices]
        molecule_slice = Molecules(infile=None,
                                   metadata=metadata_slice,
                                   coords=self.coords[  slice_indices,:])
        molecule_slice.reorganize_components()
        return molecule_slice


    # -------------------------------------------------------------------------
    # calculating geometric slices
    # -------------------------------------------------------------------------
    def gen_random_slice_point(self,slice_r):
        ''' Within boundaries, randomize slicing region'''
        return np.mean(self.coords,axis=0)[0:2]

    def calc_residue_COMS(self):
        n_res = np.unique(self.metadata.resid).size
        res_coms = np.zeros((n_res,3)) # 0 will be empty
        for i in range(n_res):
            res_coms[i,:] = np.mean(self.coords[self.resid_list[i],:],axis=0)
        return res_coms

    def rectangular_slice(self,xvals,yvals,partial_molecule='res_com'):
        '''Slices pdb to include only rectangular segment from x[0] to x[1] and
           y[0] to y[1]. Default is to exclude partial molecules, have option to
           include partial molecules or make whole and include.

           RETURNS INDICES, not actual slice
        '''
        res_coms = self.calc_residue_COMS()
        indices_tokeep = []
        all_inrange = np.asarray(np.where( (res_coms[:,0] > xvals[0]) &
                                           (res_coms[:,0] < xvals[1]) &
                                           (res_coms[:,1] > yvals[0]) &
                                           (res_coms[:,1] < yvals[1]) ))[0]
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
                    indices_tokeep.extend(self.resid_list[i])
        return np.asarray(indices_tokeep)

    def circular_slice(self,center,radius,exclude_radius=0,
                       partial_molecule='res_com'):
        indices_tokeep = []
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
                indices_tokeep.extend(self.resid_list[i])
        return np.asarray(indices_tokeep)

    # -------------------------------------------------------------------------
    # file i/o
    # -------------------------------------------------------------------------
    def read_pdb(self,pdbfile,reorganize=False):
        '''Read input pdb, default is to renumber and reorder atoms and resids
           based on leaflets (top first)
        '''
        print('Reading in PDB file')
        t = time.time()
        fid = open(pdbfile,"r")
        # initialize temporary variables
        resnum = []; xcoord = []; ycoord = []; zcoord = []; prev_res = []
        atomtype=[]; atomname = []; resname = []; chain = []; junk = []
        residcount = 0; atom_in_res = 0; prev_res = []; natoms_per_res = [];
        for pdb_line in fid:
            if pdb_line.startswith("ATOM") or pdb_line.startswith("HETATM"):
                # strings
                atomtype.append(pdb_line[ 0: 6])
                atomname.append(pdb_line[12:16])
                resname.append( pdb_line[17:21])
                chain.append(   pdb_line[21:22])
                junk.append(    pdb_line[54:  ])

               # for resnum, have to account for pdb reset at 10K
               # also, if start of residue, say how long the residue is
                curr_res = int(pdb_line[22:26])
                if curr_res != prev_res:

                    atom_in_res = 0
                    residcount += 1
                    prev_res = curr_res
                resnum.append(residcount)

                # floats for coordinate array
                xcoord.append(float(pdb_line[30:38]))
                ycoord.append(float(pdb_line[38:46]))
                zcoord.append(float(pdb_line[46:54]))

        fid.close()
        # now turn numbers into numpy
        resid =  np.array(resnum)
        atomno = np.arange(resid.size) + 1
        zero_array = np.zeros(resid.size,dtype=int)

        self.metadata = pd.DataFrame(data={'atomno':atomno,
                                            'resid':resid,
                                     'resid_length':zero_array,
                                         'leaflets':zero_array,
                                         'atomtype':atomtype,
                                         'atomname':atomname,
                                          'resname':resname,
                                            'chain':chain,
                                             'junk':junk})

        x = np.array(xcoord) # this is a stupid way to concatenate things,
        y = np.array(ycoord) # but I suck with lists
        z = np.array(zcoord)
        self.coords = np.stack((x,y,z),axis=1)

        # assigning resid lengths and leaflets
        self.assign_resids()
        self.reorganize_components(reset_leaflets=True)

        print('Finished reading in PDB file with {} atoms,time elapsed = {:.1f} seconds'.format(atomno.size,time.time()-t))
        if reorganize:
            t = time.time(); print('Reformatting PDB input')
            self.reorganize_components()
            print('Finished formatting PDB input, time required = {:.1f} seconds'.format(time.time()-t))

    def write_pdb(self,outfile,position='positive'):
        print('Writing out PDB file')
        t = time.time()
        '''Outputs to pdb file, CRYST1 and ATOM lines only'''
        # will need modulus for atomno(max=99,999) AND resno (max = 9,999)
        fout = open(outfile,'w')

        # bring everything to positive regime, allows for up to 9999 angstroms
        # in PDB file format
        if position == 'positive':
            out_coords = np.copy(self.coords -np.min(self.coords,axis=0))
        elif position == 'positive_xy':
            out_coords = np.copy(self.coords)
            out_coords[:,0:2] = out_coords[:,0:2] - np.min(out_coords[:,0:2],axis=0)
        elif position == 'center':
            out_coords = np.copy(self.coords - np.mean(self.coords,axis=0))
        elif position == 'center_xy':
            out_coords = np.copy(self.coords)
            out_coords[:,0:2] = out_coords[:,0:2] - np.mean(out_coords[:,0:2],axis=0)
        else:
            out_coords = np.copy(self.coords)

        nparts = self.coords.shape[0]
        dims = self.get_current_dims()
        # write out box dims
        line = 'CRYST1{0:9.3f}{1:9.3f}{2:9.3f}{3:7.2f}{3:7.2f}{3:7.2f}'.format(
                dims[0],dims[1],dims[2],90)
        fout.write(line + '\n')

        # dataframes are no good at individual access for millions of times,
        # so convert to a dictionary of lists, speedup of each access is
        # like in the 50 ns range vs 15 us for the dataframe
        '''
        dict_out = self.metadata.to_dict('list')
        atomtype =

        fout.writelines("{:6s}{:5d} {:4s} {:4s}{:1s}{:4d}    {:8.3f}{:8.3f}{:8.3f}{:s}\n".format(
                                     dict_out['atomtype'][i],
                                     np.mod(i+1,100000), # cap at 10^5
                                     dict_out['atomname'][i],
                                     dict_out['resname'][i],
                                     dict_out['chain'][i],
                                     np.mod(dict_out['resid'][i],10000),   # cap at 10^4
                                     out_coords[i,0],
                                     out_coords[i,1],
                                     out_coords[i,2],
                                     dict_out['junk'][i]) for i in range(nparts))
        '''
        atomtype = list(self.metadata.atomtype)
        atomname = list(self.metadata.atomname)
        resname  = list(self.metadata.resname)
        resid    = list(self.metadata.resid)
        chain    = list(self.metadata.chain)
        junk     = list(self.metadata.junk)
        fout.writelines("{:6s}{:5d} {:4s} {:4s}{:1s}{:4d}    {:8.3f}{:8.3f}{:8.3f}{:s}".format(
                                     atomtype[i],
                                     np.mod(i+1,100000), # cap at 10^5
                                     atomname[i],
                                     resname[i],
                                     chain[i],
                                     np.mod(resid[i],10000),   # cap at 10^4
                                     out_coords[i,0],
                                     out_coords[i,1],
                                     out_coords[i,2],
                                     junk[i]) for i in range(nparts))
        fout.close()
        print('Finished writing PDB file with {} atoms,time elapsed = {:.1f} seconds'.format(nparts,time.time()-t))

    def write_topology(self,outfile):
        '''Writes out simple topology file (.top)'''
        fout = open(outfile,'w')
        fout.write("\n\n\n[ system ]\nmemcurv system\n\n[ molecules ]\n")
        '''count = 1
        checker = 0
        curr_res = self.string_info['resname'][0]
        curr_resid= self.resid[0]
        for i in np.arange(self.resid.size):
            if self.resid[i] != curr_resid:
                res = self.string_info['resname'][i]

                if res == curr_res:
                    count += 1
                    curr_resid = self.resid[i]
                else:
                    fout.write("{:4s} {:d}\n".format(curr_res,count))
                    count = 1
                    curr_res = res
                    curr_resid = self.resid[i]

        fout.write("{:4s} {:d}\n".format(curr_res,count))
        '''

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
        fout.close()

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

        fout = open(outfile,'w')
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
        fout.close()
