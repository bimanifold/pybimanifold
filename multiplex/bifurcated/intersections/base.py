from dolfin import *
from mshr import *
from multiplex.meshelement import MeshElement
from multiplex.boundary import Boundary

class BaseIntersection(MeshElement):
    """Intersection consists of upper, inner and lower wall, two exits and one inlet"""

    def __init__(self, dim):
        super().__init__()

        # Declare varibales
        self.origin = dim[0]
        self.li = dim[1][0]
        self.lo = dim[1][1]
        self.lm = dim[1][2]
        self.wi = dim[1][3]
        self.wo = dim[1][4]
        self.wm = dim[1][5]

        #Wall Boundaries
        self.upper  = Boundary()
        self.inner  = Boundary()
        self.lower  = Boundary()

        #Opening Boundaries
        self.inlet      = Boundary()
        self.lower_exit = Boundary()
        self.upper_exit = Boundary()

        #Adding points to Boundaries
        x = self.origin[0]
        y = self.origin[1]
        li = self.li
        lo = self.lo
        lm = self.lm
        wi = self.wi
        wo = self.wo
        wm = self.wm
        x_outer = dim[2]
        y_outer = dim[3]
        x_inner = dim[4]
        y_inner = dim[5]

        self.lower.add(Point(x,y-wi/2))
        for xval,yval in zip(x_outer,y_outer):
            self.lower.add(Point(x+xval,y-yval))
        self.lower.add(Point(x+li+lo, y-lm/2))

        self.lower_exit.add(Point(x+li+lo, y-lm/2))
        self.lower_exit.add(Point(x+li+lo, y-lm/2+wo))

        self.inner.add(Point(x+li+lo, y-lm/2+wo))
        for xval,yval in zip(x_inner,y_inner):
            self.inner.add(Point(x+xval,y-yval))
        for xval,yval in zip(x_inner[::-1],y_inner[::-1]):
            self.inner.add(Point(x+xval,y+yval))
        self.inner.add(Point(x+li+lo, y+lm/2-wo))

        self.upper_exit.add(Point(x+li+lo, y+lm/2-wo))
        self.upper_exit.add(Point(x+li+lo, y+lm/2))

        self.upper.add(Point(x+li+lo, y+lm/2))
        for xval,yval in zip(x_outer[::-1],y_outer[::-1]):
            self.upper.add(Point(x+xval,y+yval))
        self.upper.add(Point(x,y+wi/2))

        self.inlet.add(Point(x,y+wi/2))
        self.inlet.add(Point(x,y-wi/2))

    def inflow(self, point):
        """Boundary of inlet"""
        return self.inlet.on_boundary(point)

    def outflow(self, point):
        """"Boundary of exits"""
        if self.lower_exit.on_boundary(point): return True
        if self.upper_exit.on_boundary(point): return True
        return False

    def walls(self, point):
        """"Boundary of walls"""
        if self.upper.on_boundary(point): return True
        if self.inner.on_boundary(point): return True
        if self.lower.on_boundary(point): return True
        return False

    # Depreciated (kept for debugging purposes during development)
    # def plot(self,color=(0,0,0)):
    #     self.upper.plot(color)
    #     self.inner.plot(color)
    #     self.lower.plot(color)

    def vertices(self):
        vertices = [self.lower.vertices, self.inner.vertices, self.upper.vertices]
        vertices = [vertex for sublist in vertices for vertex in sublist]
        return vertices

    # Depreciated (kept for debugging purposes during development)
    # def listall(self):
    #     """list all vertices"""
    #     for point in self.vertices():
    #         print(point.x(), point.y())

    def outlets(self):
        """get the y cordinates for each outlet"""

        x = self.origin[0]
        y = self.origin[1]
        li = self.li
        lo = self.lo
        lm = self.lm
        wi = self.wi
        wo = self.wo
        wm = self.wm

        return x+li+lo, [[y-lm/2, y-lm/2+wo],
                         [y+lm/2-wo, y+lm/2]];
