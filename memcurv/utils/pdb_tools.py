''' File i/o stuff '''





import numpy as np


class bilayer_pdb:

    def __init__(self):
        self.atomtype = []   #   str, 0-5
        self.atomno   = []   #   int, 6-10
        self.atomname = []   #   str, 12-15
        self.resname  = []   #   str, 16-20
        self.resid    = []   #   int, 22-25
        self.junk     = []   #   str, rest of pdb    pass
        self.coords   = []   # numpy, 3D array of nparts * xyz

    def read_pdb(self,pdbfile):
        fid = open(pdbfile,"r")
        # initialize temporary variables
        atomnum = []
        resnum  = []
        xcoord  = []
        ycoord  = []
        zcoord  = []
        for pdb_line in fid:
            if pdb_line.startswith("ATOM") or pdb_line.startswith("HETATM"):
                # strings
                self.atomtype.append(pdb_line[ 0: 6])
                self.atomname.append(pdb_line[12:16])
                self.resname.append( pdb_line[16:21])
                self.junk.append(    pdb_line[54:])
                # integers
                atomnum.append(int(pdb_line[ 6:11]))
                resnum.append( int(pdb_line[22:26]))
                # floats for coordinate array
                xcoord.append(float(pdb_line[30:38]))
                ycoord.append(float(pdb_line[38:46]))
                zcoord.append(float(pdb_line[46:54]))
        # now turn numbers into numpy
        self.atomno = np.array(atomnum)
        self.resid =  np.array(resnum)
        x = np.array(xcoord) # this is a stupid way to concatenate things,
        y = np.array(ycoord) # but I suck with lists
        z = np.array(zcoord)
        self.coords = np.stack((x,y,z),axis=1)


    def write_pdb():
        # will need modulus for atomno(max=99,999) AND resno (max = 9,999)
        pass
