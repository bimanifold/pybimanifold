from dolfin import *
from mshr import *
from multiplex.bifurcated.intersections.base import BaseIntersection

class RectangularIntersection(BaseIntersection):
    """Rectangular Intersection"""

    def __init__(self, origin, li, lo, lm, wi, wo, wm):

        x_outer = [li, li]
        y_outer = [wi/2, lm/2]
        x_inner = [li+wm]
        y_inner = [lm/2-wo]

        super().__init__([origin, [li,lo, lm, wi, wo, wm],x_outer,y_outer,x_inner,y_inner])
        self.domain = Polygon(self.vertices())
