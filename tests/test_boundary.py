import multiplex
from dolfin import *

from multiplex.boundary import Boundary

class TestBoundaryClass:

   def test_one(self):
      boundary = Boundary()
      boundary.add(Point(4,5))
      boundary.add(Point(3,4))
      boundary.listall()
      assert(len(boundary.vertices)==2)

   def test_two(self):
      boundary = Boundary()
      boundary.add(Point(5,5))
      boundary.add(Point(3,5))
      boundary.listall()
      assert(boundary.on_boundary(Point(4,5)))