''' Molecules class has coordinates, pdb i/o stuff, manipulation'''

import numpy as np
import pandas as pd
import nonrigid_coordinate_transformations as nrb
import rigid_body_transforms as rb
import time

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
        self.metadata.reindex(new_index_order)


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
    def gen_random_slice_point(self,slice_r):
        ''' Within boundaries, randomize slicing region'''
        return np.mean(self.coords,axis=0)[0:2]

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
        '''Read input pdb
        '''
        print('Reading in PDB file')
        t = time.time()
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
        print('Finished reading in PDB file with {} atoms,time elapsed = {:.1f} seconds'.format(self.coords.shape[0],time.time()-t))

    def write_pdb(self,outfile,position='positive',reorder=True):
        print('Writing out PDB file')
        t = time.time()

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
        dims = self.get_current_dims()
        # write out box dims


        with open(outfile,'w') as fout:
            line = 'CRYST1{0:9.3f}{1:9.3f}{2:9.3f}{3:7.2f}{3:7.2f}{3:7.2f}'.format(
                    dims[0],dims[1],dims[2],90)
            fout.write(line + '\n')

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


        print('Finished writing PDB file with {} atoms,time elapsed = {:.1f} seconds'.format(nparts,time.time()-t))

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
