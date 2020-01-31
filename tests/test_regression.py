from multiplex.manifold import Manifold
import numpy as np
import os
import shutil
import pkg_resources

class TestClass:
   def test_one(self):
      if os.path.exists('test_data') and os.path.isdir('multiplex'):
          shutil.rmtree('test_data')

      mymanifold = Manifold('./test_data/test_01','check', verbose=True)
      mymanifold.load(pkg_resources.resource_filename('multiplex','default.yaml'))
      mymanifold.solve()
      assert np.sum(mymanifold.get_Q()) - 1.0215314871585195e-05 < 1e-6

   def test_two(self):
      mymanifold = Manifold('./test_data/test_01','check', verbose=True)
      assert np.sum(mymanifold.get_Q()) - 1.0215314871585195e-05 < 1e-6
