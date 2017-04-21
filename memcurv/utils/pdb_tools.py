''' bilayer_pdb class has coordinates, pdb i/o stuff, manipulation'''

import numpy as np

class bilayer_pdb:

    def __init__(self,atomtype=[],atomno=[],atomname=[],resname=[],chain=[],
                resid=[],junk=[],coords=[],leaflets=[],resid_list=[]):
        # things needed to write pdb
        self.atomtype  = atomtype   #   str, 0-5
        self.atomno    = atomno     #   int, 6-10
        self.atomname  = atomname   #   str, 12-15
        self.resname   = resname    #   str, 16-20
        self.chain     = chain      #   str  22
        self.resid     = resid      #   int, 23-26
        self.junk      = junk       #   str, rest of pdb    pass
        self.coords    = coords     # numpy, 3D array of nparts * xyz
        # organizational variables
        self.leaflets  = leaflets   # 1 for top, 0 for bot
        self.resid_list= resid_list # nres size list
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
        i = 0
        self.leaflets = np.zeros_like(self.resid)
        bilayer_mean = np.mean(self.coords[:,2])
        while i < (len(self.resid)-1):
            j = i + 1
            while self.resid[j] == self.resid[i]:
                j += 1
            rescoords = self.coords[i:j,2]
            if np.mean(rescoords) >  bilayer_mean:
                self.leaflets[i:j+1] = 1
            i = j

    def renumber_atoms(self):
        ''' Renumber atoms from 1 to natoms'''
        self.atomno = np.arange(1,self.atomno.size + 1)
        self.assign_resid_list()  # update resid list
        
    def renumber_resids(self):
        ''' Renumber residues from 1 to n_residues '''
        for i in range(self.resid_list):
            resid[resid_list[i,:]]= i + 1
            self.assign_resid_list

    def reorder_byleaflet(self):
        pass

    def assign_resid_list(self):
        reslist = np.unique(self.resid)
        self.resid_list = list()
        for i in range(1,len(reslist)+1):         # make a bunch of empty lists
            self.resid_list.append(list())
        for i in range(0,len(self.atomno)):            # 0 will be unused
            self.resid_list[self.resid[i]].append(i)

    # -------------------------------------------------------------------------
    # adding and slicing pdb classes
    # -------------------------------------------------------------------------
    def append_pdb(self,new_pdb):
        '''appends all information from new_pdb to end of current pdb '''
        # strings
        self.atomtype.append(new_pdb.atomtype)
        self.atomname.append(new_pdb.atomname)
        self.resname.append( new_pdb.resname)
        self.chain.append(   new_pdb.chain)
        self.junk.append(    new_pdb.junk)
        # numpy arrays
        self.atomno   = np.append(self.atomno,new_pdb.atomno)
        self.resid    = np.append(self.resid, new_pdb.resid)
        self.coords   = np.vstack((self.coords,new_pdb.coords))
        self.leaflets = np.append(self.leaflets,new_pdb.leaflets)

    def slice_pdb(self,slice_indices):
        ''' returns a new instance of current pdb class with sliced indices'''
        return bilayer_pdb(atomtype=[self.atomtype[x] for x in slice_indices],
                           atomno=  self.atomno[  slice_indices],
                           atomname=[self.atomname[x] for x in slice_indices],
                           resname= [self.resname[x]  for x in slice_indices],
                           chain=   [self.chain[x]    for x in slice_indices],
                           resid=   self.resid[   slice_indices],
                           junk=    [self.junk[x]     for x in slice_indices],
                           leaflets=self.leaflets[slice_indices],
                           coords=  self.coords[  slice_indices,:])

    def reorder_byleaflet(self):
        pass
    # -------------------------------------------------------------------------
    # calculating geometric slices
    # -------------------------------------------------------------------------
    def rectangular_slice(self,xrange,yrange,partial_molecule='exclude'):
        '''Slices pdb to include only rectangular segment from x[0] to x[1] and
           y[0] to y[1]. Default is to exclude partial molecules, have option to
           include partial molecules or make whole and include
        '''
        xinrange=self.coords[:,0] > xrange[0] & self.coords[:,0] < xrange[1]
        yinrange=self.coords[:,1] > yrange[0] & self.coords[:,1] < yrange[1]

    def circular_slice(self,radius,partial_molecule='exclude'):
        pass

    # -------------------------------------------------------------------------
    # file i/o
    # -------------------------------------------------------------------------

    def read_pdb(self,pdbfile):
        fid = open(pdbfile,"r")
        # initialize temporary variables
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
                self.atomtype.append(pdb_line[ 0: 6])
                self.atomname.append(pdb_line[12:16])
                self.resname.append( pdb_line[16:21])
                self.junk.append(    pdb_line[54:  ])
                self.chain.append(   pdb_line[21:22])
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

    def write_pdb(self,outfile):
        # will need modulus for atomno(max=99,999) AND resno (max = 9,999)
        fout = open(outfile,'w')
        nparts = len(self.atomname)
        dims = self.get_current_dims()
        # write out box dims
        line = 'CRYST1{0:9.3f}{1:9.3f}{2:9.3f}{3:7.2f}{3:7.2f}{3:7.2f}'.format(
                dims[0],dims[1],dims[2],90)
        fout.write(line + '\n')

        for i in range(nparts):
            line =("{:6s}{:5d} {:^4s}{:1s}{:3s} {:1s}{:4d}{:1s} "
            + "{:8.3f}{:8.3f}{:8.3f}{:s}")
            fout.write(line.format(self.atomtype[i],
                                   np.mod(self.atomno[i],100000), # cap at 10^5
                                   self.atomname[i],
                                   " ", self.resname[i],
                                   self.chain[i],
                                   np.mod(self.resid[i],10000),   # cap at 10^4
                                   " ", self.coords[i,0],
                                   self.coords[i,1],
                                   self.coords[i,2],
                                   self.junk[i]))
        fout.close()
