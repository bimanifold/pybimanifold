from dolfin import *
from mshr import *
from multiplex.bifurcated.intersections.base import BaseIntersection

class DigitizedIntersection(BaseIntersection):
    """Digitized Intersection"""

    def __init__(self, origin, li, lo, lm, wi, wo, wm):

        a1,a2,a3,a4 = 0.2*wi,0.2*wi,0.2*wi,0.2*wi
        b1,b2,b3,b4 = 0.2*wo,0.2*wo,0.2*wo,0.2*wo

        x_outer = [li-a1,li-a1,li,li,li+b1,li+b1]
        y_outer = [wi/2,wi/2+a2,wi/2+a2,lm/2-b2,lm/2-b2,lm/2]
        x_inner = [li+wm+b3,li+wm+b3,li+wm,li+wm,li+wm-a3]
        y_inner = [lm/2-wo,lm/2-wo-b4,lm/2-wo-b4,a4,a4]

        super().__init__([origin,[li,lo,lm,wi,wo,wm],x_outer,y_outer,x_inner,y_inner])
        self.domain = Polygon(self.vertices())
