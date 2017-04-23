# for symmetry, faces stored in a 3X2 array (xyz * (updown/frontback/leftright))
# options are 'flat','circle','half_circle','closed'
# -----------------------------------------------------------------------------
# shapes with lateral symmetry
# -----------------------------------------------------------------------------
class plane:
    def __init__(xdim,ydim):
        self.get_area(xdim,ydim)
        self.faces = ['']
    def set_area(x,y):
        self.area =  x * y

class semi_cylinder:
    def __init__(radius,ydim):
        self.get_area(radius,ydim)
        self.lat_sym = [[1 0 0],[-1 0 0]] # x open, both sides
        self.rad_sym = [[0 1 0],[0 -1 0]] # y open, both sides
    def set_area(r,y):
        self.area = math.pi * r * y
# -----------------------------------------------------------------------------
# Closed shapes
# -----------------------------------------------------------------------------

class sphere:
    def __init__(self,radius):
        self.get_area(radius)
        self.latsym = []       # no connections available
        self.radsym = []
    def get_area(radius):
        self.area =  4 * math.pi * (radius ** 2)
class cylinder:

# -----------------------------------------------------------------------------
# Shapes with radial symmetry
# -----------------------------------------------------------------------------

class semi_sphere:
    def __init__(self,radius):
        self.get_area(radius)
        self.latsym = []   # no flat connections
        self.
    def get_area(radius):
        self.area =  2 * math.pi * (radius ** 2)
