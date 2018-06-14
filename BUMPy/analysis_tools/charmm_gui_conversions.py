import numpy as np

''' Series of scripts for converting between zo values and arbitrary charmm-gui
    size and ratio estimations
'''



# p values
# DPPE - 53.5    - io = 0.62315    - zo =  0.9416
# DPPC - 59.5    - io = 0.59147    - zo =  1.0443
# DOPE - 57.2    - io = 0.60333    - zo =  1.0052
# DOPC - 62.1    - io = 0.57850    - zo =  1.0879


def io_from_charmm(R,p,n=1.25):
    return (R ** n) / (R ** n + p ** n)

def calc_zo(radius,io_ratio,geometry,correction_factor_on=False):
    ''' Calculates pivotal plane from given radius and inner to outer lipid
        ratio.

        parameters:
            radius               - radius of sphere or cylinder
            io_ratio             - inner to outer ratio, < 1
            geometry             - 'sphere' or 'cylinder'
        returns:
            z - pivotal plane, in same units as radius
    '''
    oi_ratio = 1 / io_ratio
    if geometry == 'sphere':
        return radius * (np.sqrt(oi_ratio) - 1) / (np.sqrt(oi_ratio) + 1)
    elif geometry == 'cylinder':
        return radius * (oi_ratio -1) / (oi_ratio + 1)
