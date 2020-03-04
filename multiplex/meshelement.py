from dolfin import *
from mshr import *
import numpy as np
import time
import yaml
from os.path import isfile, isdir, split, join
from os import mkdir, makedirs,getcwd, remove
from re import findall
from datetime import datetime
import matplotlib.pyplot as plt
from abc import ABCMeta, abstractmethod

class MeshElement:

    def __init__(self,cwd=getcwd(),name=datetime.now().strftime("%d_%m_%Y_%H_%M_%S"),verbose=False):
        """
        Mesh Element

        """
        self.path              = cwd
        self.name              = name
        self.verbose           = verbose
        self.initialized       = False  # bool if the manifold is able to run
        self.mesh_isgenerated  = False
        self.mesh              = None
        self.meta_data         = None

        # Create working directory (+ path) if it not exists
        # Load data from YAML file if it exists
        if isfile(join(self.path,self.name)+'.yaml'):

            # Read from file
            if verbose==True:
                print("working directory: '"+cwd+"'")
                print(join(self.path,self.name)+'.yaml exits!')

            stream = open(join(self.path,self.name)+'.yaml', 'r')
            self.meta_data = yaml.load(stream, Loader=yaml.FullLoader)

            # here, we don't change the .h5, .xml or pvd/ folder becaue data is still valid
            self.__validate_meta_data_and_write_to_file__()

        else:
            if isdir(cwd):
                if verbose==True:
                    print("working directory: '"+cwd+"'")
            else:
                makedirs(cwd)
                if verbose==True:
                    print("working directory: '"+cwd+"' created!")



    def load(self,yamlFile):
        """
        Load meta data of manifold from YAML file and creates new yaml file with name MANIFOLD_NAME.yaml.

        Parameters:
        yamlFile (str):     yaml file path / name 
        
        Returns:
        None

        Example:
        mymanifold.load("sample.yaml")

        Here the "sample.yaml" is inside the current working directory.
        """

        if isfile(join(self.path,self.name)+'.xml'):
            raise Exception(".xml file already exits for this instance in current woring directory.")

        if isfile(join(self.path,self.name)+'.h5'):
            raise Exception(".h5 file already exits for this instance in current woring directory.")

        if isdir(join(self.path,self.name)):
            raise Exception("PVD folder already exits for this instance in current woring directory.")

        stream = open(yamlFile, 'r')
        self.meta_data = yaml.load(stream, Loader=yaml.FullLoader)

        self.__validate_meta_data_and_write_to_file__()



    def change(self,item,new_value):
        """
        Convinent function to change a single parameter after laoding from YAML file.

        This function automatically updates the corresponding NAME.yaml file, NAME.xml file
        and deltes the NAME.h5 file. If you only want to update the number of snap shot please refer to solve(snaps). 
        solve(snaps) keeps the mesh file NAME.xml and data file NAME.h5 unchanged.

        Parameters:
        item (str):             name of meta_data item to change in yaml file
        new_value (str,int):    new value for item

        Returns:
        None

        Examples:
        mymanifold.change("inlet_width", 5)
        mymanifold.change("mesh_type", "curved")

        This changes the inlet width to 5 and the mesh type to the cuved manifold.

        """

        if self.initialized==False:
            raise Exception("Please initialize all parameters with load() before changing a parameter.")

        if item not in self.meta_data.keys():
            raise Exception("Item to change not found in object meta data. Please check the spelling in the YAML file")

        self.meta_data[item] = new_value

        # delete .h5 and .xml file and ensure that no pvd folder exists
        if isdir(join(self.path,self.name)):
            raise Exception("PVD folder already exits for this instance in current woring directory.")

        if isfile(join(self.path,self.name)+'.xml'):
            remove(join(self.path,self.name)+'.xml')

        if isfile(join(self.path,self.name)+'.h5'):
            remove(join(self.path,self.name)+'.h5')

        self.__validate_meta_data_and_write_to_file__()



    def __validate_meta_data_and_write_to_file__(self):

        if 'iterations' not in self.meta_data.keys():
            raise Exception("YAML file needs to contain file 'iterations'")
        if 'snaps' not in self.meta_data.keys():
            raise Exception("YAML file needs to contain file 'snaps'")
        if 'mesh_resolution' not in self.meta_data.keys():
            raise Exception("YAML file needs to contain file 'mesh_resolution'")
        if 'time_step' not in self.meta_data.keys():
            raise Exception("YAML file needs to contain file 'time_step'")
        if 'viscosity' not in self.meta_data.keys():
            raise Exception("YAML file needs to contain file 'viscosity'")
        if 'mass_density' not in self.meta_data.keys():
            raise Exception("YAML file needs to contain file 'mass_density'")
        if 'reynolds_number' not in self.meta_data.keys():
            raise Exception("YAML file needs to contain file 'reynolds_number'")

        self.__initialize_domain__()
        self.initialized = True

        # Write data to file
        with open(join(self.path,self.name)+'.yaml', 'w') as file:
            documents = yaml.dump(self.meta_data, file)

        # Print data to screen
        if self.verbose==True:
            snaps        = self.meta_data['snaps']
            iterations   = self.meta_data['iterations']
            resolution   = self.meta_data['mesh_resolution']
            dt           = self.meta_data['time_step']
            mu           = self.meta_data['viscosity']
            rho          = self.meta_data['mass_density']
            Re           = self.meta_data['reynolds_number']

            print("sanps:\t\t\t\t\t", snaps)
            print("iterations:\t\t\t\t",  iterations)
            print("mesh_resolution:\t\t\t",resolution)
            print("time_step (seconds):\t\t\t", dt)
            print("viscosity (Pa):\t\t\t\t", mu)
            print("mass density (kg/cubic meter):\t\t", rho)
            print("Reynolds number:\t\t\t", Re)

        self.__print_domain_meta_data_to_screen__()

        if isfile(join(self.path,self.name)+'.xml'):

            self.mesh = Mesh(join(self.path,self.name)+'.xml')
            self.mesh_isgenerated = True;
            if self.verbose:
                print("")
                print("Mesh loaded from file: ", join(self.path,self.name)+'.xml')

        else:

            self.mesh = generate_mesh(self.domain, resolution)
            self.mesh_isgenerated = True
            File(join(self.path,self.name)+'.xml') << self.mesh
            if self.verbose:
                print("")
                print("Mesh newly generated!")



    def plot(self,filename="foo.pdf"):
        """
        Plot the mesh as a collection of dots.
        Inflow dots appear in green.
        Outlfow dots appear in blue.
        Wall dots appear in red. 

        Parameters:
        filename (str):         filename of PDF file, will be saved in current working directory

        Returns:
        None

        Example:
        mymanifold.plot("foo.pdf")

        """
        import matplotlib.pyplot as plt

        if self.mesh_isgenerated == False:
            raise Exception("MESH IS NOT GENERATED")

        plt.figure()

        mesh_coordinates = self.mesh.coordinates()
        grey = np.array([0.8,0.8,0.8])
        plt.plot(mesh_coordinates[:,0], mesh_coordinates[:,1],',', color=grey);

        for v in vertices(self.mesh):
            x = v.point().x()
            y = v.point().y()
            if self.inflow([x,y]):  plt.plot(x,y,',', color='green')
            if self.outflow([x,y]): plt.plot(x,y,',', color='blue')
            if self.walls([x,y]):   plt.plot(x,y,',', color='red')

        plt.axis('equal');
        plt.savefig(join(self.path,filename), bbox_inches='tight')


    def solve(self,snaps=None):
        """

        Solving time dependent Navier Stokes equation for a given number of snap shots.
        Creates a .h5 file with the solution in the current working directory. 

        If the number of snaps is not the same as in YAML file, the higher number is chosen and the YAML file is updated.

        If a prior .h5 file already exists, the simulation picks up at the final snap shot number.

        Parameters:
        snaps (int):    number of snaps for the simulation to run

        Returns:
        None

        Example:
        >>> mymanifold.BifurcatedManifold()
        >>> mymanifold.load("sample.yaml")
        >>> mymanifold.change('iterations', 30)
        >>> mymanifold.solve()                          # run snaps 0  to 30 
        >>> mymanifold.solve(40)                        # run snaps 31 to 40 

        """

        #Referece: https://github.com/hplgit/fenics-tutorial/blob/master/pub/python/vol1/ft08_navier_stokes_cylinder.py

        if self.initialized == False:
            raise RuntimeError("Manifold is not initialized!")

        meta_data = self.meta_data
        width  = float(meta_data['inlet_width'])
        d      = float(meta_data['inlet_length'])
        ddd    = float(meta_data['device_to_device_distance'])
        layers = int(meta_data['number_of_layers'])
        scale  = float(meta_data['scale'])
        gamma  = float(meta_data['gamma'])
        origin = meta_data['origin']

        mesh_type   = meta_data['mesh_type']
        res         = int(meta_data['mesh_resolution'])
        snapsF      = int(meta_data['snaps'])
        num_iter    = int(meta_data['iterations'])
        dt          = float(meta_data['time_step'])
        mu          = float(meta_data['viscosity'])
        rho         = float(meta_data['mass_density'])
        Re          = float(meta_data['reynolds_number'])

        if(snaps != None):
            if(snapsF < snaps):
                self.meta_data['snaps'] = snaps
                # here we don't delete .h5, .xml and pvd/ folder bc data is still valid
                self.__validate_meta_data_and_write_to_file__()

        num_snaps = self.meta_data['snaps']

        v_mean = Re*mu/rho/width
        inflow_profile=(str(v_mean*scale*3/2/(width*scale/2)**2)+'*(x[1]+'+str(width*scale/2)+')*('+str(width*scale/2)+'-x[1])', '0')


        # Rescale mu and rho 
        mu = mu/scale
        rho = rho/scale**3


        def inflow(x):
            return self.inflow(x)

        def walls(x):
            return self.walls(x)

        def outflow(x):
            return self.outflow(x)


        # Define function spaces
        V = VectorFunctionSpace(self.mesh, 'P', 2)
        Q = FunctionSpace(self.mesh, 'P', 1)


        # Define boundary conditions
        bcu_inflow = DirichletBC(V, Expression(inflow_profile, degree=2),inflow)
        bcu_walls = DirichletBC(V, Constant((0, 0)), walls)
        bcp_outflow = DirichletBC(Q, Constant(0), outflow)
        bcu = [bcu_inflow, bcu_walls]
        bcp = [bcp_outflow]

        # Define trial and test functions
        u = TrialFunction(V)
        v = TestFunction(V)
        p = TrialFunction(Q)
        q = TestFunction(Q)

        # Define functions for solutions at previous and current time steps
        u_n = Function(V)
        u_  = Function(V)
        p_n = Function(Q)
        p_  = Function(Q)

        # Define expressions used in variational forms
        U  = 0.5*(u_n + u)
        n  = FacetNormal(self.mesh)
        f  = Constant((0, 0))
        k  = Constant(dt)
        mu = Constant(mu)
        rho = Constant(rho)

        # Define symmetric gradient
        def epsilon(u):
            return sym(nabla_grad(u))

        # Define stress tensor
        def sigma(u, p):
            return 2*mu*epsilon(u) - p*Identity(len(u))

        # Define variational problem for step 1
        F1 = rho*dot((u - u_n) / k, v)*dx \
            + rho*dot(dot(u_n, nabla_grad(u_n)), v)*dx \
            + inner(sigma(U, p_n), epsilon(v))*dx \
            + dot(p_n*n, v)*ds - dot(mu*nabla_grad(U)*n, v)*ds \
            - dot(f, v)*dx
        a1 = lhs(F1)
        L1 = rhs(F1)

        # Define variational problem for step 2
        a2 = dot(nabla_grad(p), nabla_grad(q))*dx
        L2 = dot(nabla_grad(p_n), nabla_grad(q))*dx - (1/k)*div(u_)*q*dx

        # Define variational problem for step 3
        a3 = dot(u, v)*dx
        L3 = dot(u_, v)*dx - k*dot(nabla_grad(p_ - p_n), v)*dx

        # Assemble matrices
        A1 = assemble(a1)
        A2 = assemble(a2)
        A3 = assemble(a3)

        # Apply boundary conditions to matrices
        [bc.apply(A1) for bc in bcu]
        [bc.apply(A2) for bc in bcp]

        if isfile(join(self.path,self.name)+'.h5'):
            hdf = HDF5File(self.mesh.mpi_comm(), join(self.path,self.name)+'.h5', "r")
            attr = hdf.attributes("p")
            dataset = "u/vector_%d"%(attr['count']-2)
            hdf.read(u_n,dataset)
            dataset = "p/vector_%d"%(attr['count']-2)
            hdf.read(p_n,dataset)
            dataset = "u/vector_%d"%(attr['count']-1)
            hdf.read(u_,dataset)
            dataset = "p/vector_%d"%(attr['count']-1)
            start_index = attr['count']-1
            hdf.read(p_,dataset)
            hdf.close()
            Hdf=HDF5File(self.mesh.mpi_comm(), join(self.path,self.name)+'.h5', "a")
        else:
            Hdf=HDF5File(self.mesh.mpi_comm(), join(self.path,self.name)+'.h5', "w");
            Hdf.write(self.mesh, "mesh")
            Hdf.write(u_, "u", 0)
            Hdf.write(p_, "p", 0)
            start_index = 0

        textfile = join(self.path,self.name)+'.txt'

        # Time-stepping
        t = start_index*num_iter*dt
        start = time.time()
        for n in range(start_index, num_snaps):

            for l in range(num_iter):
                # Update current time
                t += dt

                # Step 1: Tentative velocity step
                b1 = assemble(L1)
                [bc.apply(b1) for bc in bcu]
                solve(A1, u_.vector(), b1, 'bicgstab', 'hypre_amg')

                # Step 2: Pressure correction step
                b2 = assemble(L2)
                [bc.apply(b2) for bc in bcp]
                solve(A2, p_.vector(), b2, 'bicgstab', 'hypre_amg')

                # Step 3: Velocity correction step
                b3 = assemble(L3)
                solve(A3, u_.vector(), b3, 'cg', 'sor')

                # Update previous solution
                u_n.assign(u_)
                p_n.assign(p_)

            # Display error if required
            if self.verbose:
                print('n:', n, 't:', np.round(t,3), 'u max:', u_.vector().vec().max(),'time: ',time.time()-start)

            file1 = open(textfile,"a+")
            file1.write('n: '+str(n)+', t: '+str(np.round(t,3))+', u max: '+str(u_.vector().vec().max())+', time: '+str(time.time()-start)+'\n')
            file1.close()

            Hdf.write(u_, "u", t)
            Hdf.write(p_, "p", t)
            Hdf.flush()

        Hdf.close()

    def hdf2pvd(self):
        """
        Converts HDF file in the current working directory to PVD files.

        Parameters:
        None

        Returns:
        None

        Examle:
        >>> mymanifold.hdf2pvd()

        """
        if self.initialized == False or isfile(join(self.path,self.name)+".h5")==False:
            raise RuntimeError("Manifold is not initialized or .h5 file not found!")

        destination = join(self.path,self.name)
        source      = self.path

        if not isdir(destination):
            mkdir(destination)

        # Define function spaces
        V = VectorFunctionSpace(self.mesh, 'P', 2)
        Q = FunctionSpace(self.mesh, 'P', 1)

        u_  = Function(V)
        p_  = Function(Q)

        file_u = File(destination+'/velocity.pvd')
        file_p = File(destination+'/pressure.pvd')

        hdf = HDF5File(self.mesh.mpi_comm(), join(source,self.name)+'.h5', "r")
        attr = hdf.attributes("p")

        nsteps = attr['count']
        for i in range(nsteps):

            dataset = "p/vector_%d"%i
            hdf.read(p_,dataset)
            dataset = "u/vector_%d"%i
            hdf.read(u_,dataset)

            attr = hdf.attributes(dataset)
            t = attr['timestamp']
            file_p << (p_, t)
            file_u << (u_, t)

        if self.verbose:
            print("")
            print("Success HDF to PVD")
            print("")

        hdf.close()

    def getVelocity(self, timestep=None):
        """
        Retrieve an interpolation function of the velocity profile

        Parameter: 
        timstep (int): Snap shot number for the velocity profile

        Returns:
        u (dolfin.function.function): Velocity profile

        Example:

        >>>  u = channel.getVelocity()
        >>>  u(0,1)
        np.array([1.354545, 0.0012])

        This is the velocity vector [v_x, v_y] at point (x=0,y=1) at the final time. It is also possible to retrieve the velocity profile at a different timestep. 

        >>>  u = channel.getVelocity(15)
        >>>  u(0,1)
        np.array([1.354545, 0.0032])

        This retrieves the velocity profile at iteration timestep=15. To find the the correponding time use get_time_intervals() function.
        """

        V = VectorFunctionSpace(self.mesh, 'P', 2)
        Q = FunctionSpace(self.mesh, 'P', 1)

        u_  = Function(V)
        p_  = Function(Q)

        hdf = HDF5File(self.mesh.mpi_comm(), join(self.path+self.name)+'.h5', "r")
        attr = hdf.attributes("p")

        if (timestep==None):
            timestep = attr['count']-1

        dataset = "u/vector_%d"%timestep
        hdf.read(u_,dataset)
        dataset = "p/vector_%d"%timestep
        hdf.read(p_,dataset)

        hdf.close()

        return interpolate(u_, V)

    def get_time_intervals(self):
        """
        Retrieve the time intervals of the simulation

        Parameter: 
        None

        Returns:
        times (int list): List with the time steps

        Example:

        >>> times = channel.get_time_intervals()
        >>>  times
        [1.2, 2.4, 3.6, ....]

        The number of snap shots of the simulation can be obtained by len(times).
        """

        V = VectorFunctionSpace(self.mesh, 'P', 2)
        u_  = Function(V)
        times = []

        hdf = HDF5File(self.mesh.mpi_comm(), join(self.path+self.name)+'.h5', "r")
        attr = hdf.attributes("p")

        nsteps = attr['count']
        for i in range(nsteps):

            dataset = "u/vector_%d"%i
            hdf.read(u_,dataset)

            attr = hdf.attributes(dataset)
            t = attr['timestamp']
            times.append(t)

        hdf.close()

        return times
