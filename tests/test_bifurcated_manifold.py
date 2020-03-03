#from multiplex import BifurcatedManifold
import multiplex
import numpy as np
import os
import shutil
import pkg_resources
from os.path import isdir, isfile
import pytest

class TestBifurcatedManifoldClass:

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
      assert np.sum(mymanifold.get_Q()) > 0 
      
   def test_three(self):
      mymanifold = multiplex.BifurcatedManifold('./test_data/','test03', verbose=True)
      mymanifold.load(pkg_resources.resource_filename('multiplex','default.yaml'))
      mymanifold.change('mesh_type','digitized')
      mymanifold.solve()
      assert np.sum(mymanifold.get_Q()) > 0 

   def test_four_1(self):
      mymanifold = multiplex.BifurcatedManifold('./test_data/','test04', verbose=True)
      mymanifold.load(pkg_resources.resource_filename('multiplex','default.yaml'))
      mymanifold.change('mesh_type','triangular')
      mymanifold.solve()
      assert np.sum(mymanifold.get_Q()) > 0

   def test_four_2(self):
      mymanifold = multiplex.BifurcatedManifold('./test_data/','test04', verbose=True)
      assert len(mymanifold.get_time_intervals()) > 0 

   def test_four_3(self):
      mymanifold = multiplex.BifurcatedManifold('./test_data/','test04', verbose=True)
      u = mymanifold.getVelocity(10)
      assert np.mean(u(0,0)) > 0 

   def test_four_4(self):
      mymanifold = multiplex.BifurcatedManifold('./test_data/','test04', verbose=True)
      u = mymanifold.getVelocity()
      assert np.mean(u(0,0)) > 0 

   def test_five_1(self):
      mymanifold = multiplex.BifurcatedManifold('./test_data/','test05', verbose=True)
      mymanifold.load(pkg_resources.resource_filename('multiplex','default.yaml'))
      mymanifold.change('mesh_type','curved')
      mymanifold.solve()
      assert np.sum(mymanifold.get_Q()) > 0 

   def test_five_2(self):
      mymanifold = multiplex.BifurcatedManifold('./test_data/','test05', verbose=True)
      mymanifold.hdf2pvd()
      assert isdir("./test_data/test05/")

   def test_five_3(self):
      mymanifold = multiplex.BifurcatedManifold('./test_data/','test05', verbose=True)
      assert mymanifold.get_Qin() > 0

   def test_five_4(self):
      mymanifold = multiplex.BifurcatedManifold('./test_data/','test05', verbose=True)
      mymanifold.solve(50)
      assert np.sum(mymanifold.get_Q()) > 0

   def test_five_4(self):
      mymanifold = multiplex.BifurcatedManifold('./test_data/','test05', verbose=True)
      mymanifold.solve(50)
      assert np.sum(mymanifold.get_Q()) > 0

   ################# Exception tests #################

   def test_five_5(self):
      mymanifold = multiplex.BifurcatedManifold('./test_data/','test05', verbose=True)
      mymanifold.plot()
      assert(isfile("./test_data/foo.pdf"))
      with pytest.raises(Exception):
         mymanifold.change('inlet_width0', 5)

   def test_five_6(self):
      mymanifold = multiplex.BifurcatedManifold('./test_data/','test05', verbose=True)
      with pytest.raises(Exception):
         mymanifold.load(pkg_resources.resource_filename('multiplex','default.yaml'))

   def test_six(self):
      mymanifold = multiplex.BifurcatedManifold('./test_data/','test06',verbose=True)
      with pytest.raises(Exception):
         mymanifold.change('inlet_width', 5)
      with pytest.raises(Exception):
         mymanifold.plot()
      with pytest.raises(Exception):
         mymanifold.solve()
      with pytest.raises(Exception):
         mymanifold.hdf2pvd()
      with pytest.raises(Exception):
         mymanifold.get_Q()
      with pytest.raises(Exception):
         mymanifold.get_Qin()
      mymanifold.load(pkg_resources.resource_filename('multiplex','default.yaml'))
      with pytest.raises(Exception):
         mymanifold.change('mesh_type', 'rectangular_false')

   