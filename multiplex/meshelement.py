from dolfin import *
from mshr import *
import numpy as np
import time
from os.path import isfile, join
import matplotlib.pyplot as plt

class MeshElement:

    def __init__(self):
        self.mesh_isgenerated = False
        self.mesh = None

    def generate_mesh(self, resolution):
        self.mesh = generate_mesh(self.domain, resolution)
        self.mesh_isgenerated = True

    def mesh(self):
        if self.mesh_isgenerated == False:
            raise Exception("MESH IS NOT GENERATED")
        return self.mesh

    def domain(self):
        return self.domain

    def plot(self,filename="foo.pdf"):

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

    def solver(self,path, hdfFile, num_snaps, num_iter, mu, rho, dt, inflow_profile, overwrite, verbose=True, consecutive_run_allowed=True):

        """it is based on https://github.com/hplgit/fenics-tutorial/blob/master/pub/python/vol1/ft08_navier_stokes_cylinder.py"""

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

        # Create XDMF files for visualization output
        # file_u = File(path+'/velocity.pvd')
        # file_p = File(path+'/pressure.pvd')

        # Save mesh to file (for use in reaction_system.py)
        #File(path+'/channel.xml.gz') << self.mesh
        first_run = True;
        if isfile(join(path,hdfFile)+'.h5'):
            first_run = False;

        if first_run==False and overwrite==False:
            if consecutive_run_allowed==True:
                hdf = HDF5File(self.mesh.mpi_comm(), join(path,hdfFile)+'.h5', "r")
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
            else:
                raise Exception("HDF file already exists")

        if first_run == True:
            Hdf=HDF5File(self.mesh.mpi_comm(), join(path,hdfFile)+'.h5', "w");
            Hdf.write(self.mesh, "mesh")
            Hdf.write(u_, "u", 0)
            Hdf.write(p_, "p", 0)
            start_index = 0
        else:
            Hdf=HDF5File(self.mesh.mpi_comm(), join(path,hdfFile)+'.h5', "a")


        # Hdf.write(mu, "mu")
        # Hdf.write(rho, "rho")
        #comm = self.mesh.mpi_comm()
        #mpiRank = MPI.rank(comm)
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

            # Display error from master slave rank == 0
            if verbose:
                print('n:', n, 't:', np.round(t,3), 'u max:', u_.vector().vec().max(),'time: ',time.time()-start)

            file1 = open(textfile,"a+")
            file1.write('n: '+str(n)+', t: '+str(np.round(t,3))+', u max: '+str(u_.vector().vec().max())+', time: '+str(time.time()-start)+'\n')
            file1.close()


            # PVD file
            """https://fenicsproject.org/qa/6675/hdf-file-read-write-for-time-series/"""
            # file_u << (u_, t)
            # file_p << (p_, t)
            Hdf.write(u_, "u", t)
            Hdf.write(p_, "p", t)
            Hdf.flush()

        Hdf.close()

    def hdf_to_pvd(self, destination_path, source_path, hdfFile, verbose):

        # Define function spaces
        V = VectorFunctionSpace(self.mesh, 'P', 2)
        Q = FunctionSpace(self.mesh, 'P', 1)

        u_  = Function(V)
        p_  = Function(Q)

        file_u = File(destination_path+'/velocity.pvd')
        file_p = File(destination_path+'/pressure.pvd')

        #xdmffile_u = XDMFFile('navier_stokes_cylinder/velocity.xdmf')
        #xdmffile_p = XDMFFile('navier_stokes_cylinder/pressure.xdmf')

        hdf = HDF5File(self.mesh.mpi_comm(), join(source_path,hdfFile)+'.h5', "r")
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
            #xdmffile_u.write(u_, t)
            #xdmffile_p.write(p_, t)

        if verbose:
            print("")
            print("Success HDF to PVD")
            print("")

        hdf.close()

    def get_u(self, step=None):

        V = VectorFunctionSpace(self.mesh, 'P', 2)
        Q = FunctionSpace(self.mesh, 'P', 1)

        u_  = Function(V)
        p_  = Function(Q)

        hdf = HDF5File(self.mesh.mpi_comm(), join(self.path+self.name)+'.h5', "r")
        attr = hdf.attributes("p")

        if (step==None):
            step = attr['count']-1

        dataset = "u/vector_%d"%step
        hdf.read(u_,dataset)
        dataset = "p/vector_%d"%step
        hdf.read(p_,dataset)

        hdf.close()

        return interpolate(u_, V)

    def get_meta_data(self):

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

        return times, nsteps
