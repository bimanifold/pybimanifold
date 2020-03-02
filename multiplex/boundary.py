from dolfin import *
from mshr import *
import numpy as np
import matplotlib.pyplot as plt

class Boundary:
    """Handels the the interstion walls, inlets and exits"""

    def __init__(self):
        self.vertices = []

    def add(self, point):
        """Add Dolfin point to vertices"""
        self.vertices.append(point)

    def on_boundary(self, point):
        """Check if point is on Linear Fit between vertices"""
        for index in range(1, len(self.vertices)):
            x,y = point[0], point[1]
            x0,y0 = self.vertices[index-1].x(), self.vertices[index-1].y()
            x1,y1 = self.vertices[index].x(),   self.vertices[index].y()
            if near(x0,x1, 1e-15):
                """Check if line is vertical;
                Avoid divison by zero"""
                if near(x0, x, 1e-15):
                    "Check if point is on Vertical Line"
                    if y <= np.max([y0, y1]):
                        if np.min([y0,y1]) <= y:
                            return True
            else:
                "Check if point is on Linear Fit between 2 vertices"
                m = (y1-y0)/(x1-x0)
                if near(y, m*(x-x1)+y1, 1e-15):
                    return True
        return False

    # def vertices(self):
    #     return self.vertices

    # def listall(self):
    #     """list alll vertices"""
    #     for point in self.vertices:
    #         print(point.x(), point.y())

    # def plot(self,color=(0,0,0)):

    #     """list boundary plot"""
    #     xval = []
    #     yval = []
    #     for point in self.vertices:
    #         xval.append(point.x())
    #         yval.append(point.y())
    #     plt.plot(xval, yval, color=color)
    #     plt.axis('equal')
