import numpy as np

def gen_random_coords(nparts,xrange,yrange,zrange):
    ''' generates a random set of cartesian coordinates, with number of
        particles=nparts, and values ranging from 0 to range
    '''
    coords = np.random.rand(nparts,3)
    return coords * np.array([xrange,yrange,zrange])
