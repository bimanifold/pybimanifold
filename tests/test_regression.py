#from multiplex import BifurcatedManifold
import multiplex
import numpy as np
import os
import shutil
import pkg_resources
from os.path import isdir, isfile

class TestClass:
   def test_one(self):
      if os.path.exists('test_data') and os.path.isdir('multiplex'):
          shutil.rmtree('test_data')

      mymanifold = multiplex.BifurcatedManifold('./test_data/','test01', verbose=True)
      mymanifold.load(pkg_resources.resource_filename('multiplex','default.yaml'))
      mymanifold.solve()
      assert np.sum(mymanifold.get_Q()) - 1.0215314871585195e-05 < 1e-6

   # Test changing the mass density
   def test_two(self):
      mymanifold = multiplex.BifurcatedManifold('./test_data/','test02', verbose=True)
      mymanifold.load(pkg_resources.resource_filename('multiplex','default.yaml'))
      mymanifold.change('mass_density',500)
      mymanifold.solve()
      assert np.sum(mymanifold.get_Q()) > 0 #- 1.0215314871585195e-05 < 1e-6

   def test_three(self):
      mymanifold = multiplex.BifurcatedManifold('./test_data/','test03', verbose=True)
      mymanifold.load(pkg_resources.resource_filename('multiplex','default.yaml'))
      mymanifold.change('mesh_type','digitized')
      mymanifold.solve()
      assert np.sum(mymanifold.get_Q()) > 0 #- 1.0215314871585195e-05 < 1e-6

   def test_four(self):
      mymanifold = multiplex.BifurcatedManifold('./test_data/','test04', verbose=True)
      mymanifold.load(pkg_resources.resource_filename('multiplex','default.yaml'))
      mymanifold.change('mesh_type','triangular')
      mymanifold.solve()
      assert np.sum(mymanifold.get_Q()) > 0 #- 1.0215314871585195e-05 < 1e-6

   def test_five(self):
      mymanifold = multiplex.BifurcatedManifold('./test_data/','test05', verbose=True)
      mymanifold.load(pkg_resources.resource_filename('multiplex','default.yaml'))
      mymanifold.change('mesh_type','curved')
      mymanifold.solve()
      assert np.sum(mymanifold.get_Q()) > 0 #- 1.0215314871585195e-05 < 1e-6
      mymanifold.hdf2pvd()
      assert isdir("./test_data/test05/")
      mymanifold.plot()
      assert isfile("./test_data/foo.pdf")

   # def test_six(self):
   #    mymanifold.hdf2pvd()
   #    assert isdir("./test_data/test05/")
   #    del mymanifold



   # def test_six(self):
   #    mymanifold.hdf2pvd()
   #    assert isdir("./test_data/test05")