from dolfin import *
from mshr import *
from multiplex.intersection import Intersection
import numpy as np

class Ellipse():
    def __init__(self, x, y, a, b):
        """Ellipse with center (x,y) and axis a,b"""
        self.x = x
        self.y = y
        self.a = a
        self.b = b

    def get_arc_verticies_1(self,n):
        """return arc verticies with postive x and negative y"""
        theta = np.linspace(3*pi/2,2*pi,n+2)[1:-1]
        return list(self.a*np.cos(theta)+self.x), list(self.b*np.sin(theta)+self.y)

    def get_arc_verticies_2(self,n):
        """return arc verticies with negative x and positive y"""
        theta = np.linspace(pi,pi/2,n+2)[1:-1]
        return list(self.a*np.cos(theta)+self.x), list(self.b*np.sin(theta)+self.y)


class CurvedIntersection(Intersection):
    """Curved Intersection"""

    def __init__(self, origin, li, lo, lm, wi, wo, wm):

        a1,a2,a3,a4 = 0.4*wi,0.4*wi,0.4*wi,0.4*wi
        b1,b2,b3,b4 = 0.4*wo,0.4*wo,0.4*wo,0.4*wo

        resolution = 10

        x_outer = [li-a1]
        y_outer = [wi/2]
        P1 = Ellipse(li-a1,wi/2+a2,a1,a2)
        x,y = P1.get_arc_verticies_1(resolution)
        x_outer.extend(x)
        y_outer.extend(y)
        x_outer.extend([li,li])
        y_outer.extend([wi/2+a2,lm/2-b2])
        P2 = Ellipse(li+b1,lm/2-b2,b1,b2)
        x,y = P2.get_arc_verticies_2(resolution)
        x_outer.extend(x)
        y_outer.extend(y)
        x_outer.extend([li+b1])
        y_outer.extend([lm/2])

        x_inner = [li+wm+b3]
        y_inner = [lm/2-wo]
        P3 = Ellipse(li+wm+b3,lm/2-wo-b4,b3,b4)
        x,y = P3.get_arc_verticies_2(resolution)
        x_inner.extend(x[::-1])
        y_inner.extend(y[::-1])
        x_inner.extend([li+wm,li+wm])
        y_inner.extend([lm/2-wo-b4,a4])
        P4 = Ellipse(li+wm-a3,a4,a3,a4)
        x,y = P4.get_arc_verticies_1(resolution)
        x_inner.extend(x[::-1])
        y_inner.extend(y[::-1])


        super().__init__([origin, [li,lo, lm, wi, wo, wm],x_outer,y_outer,x_inner,y_inner])
        self.domain = Polygon(self.vertices())
