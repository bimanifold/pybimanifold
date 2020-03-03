from dolfin import *
from mshr import *
import numpy as np
import yaml
from os.path import isfile, isdir, split,join
from os import mkdir, makedirs,getcwd
from re import findall
from datetime import datetime
from multiplex.bifurcated.intersections.rectangular import RectangularIntersection as RI
from multiplex.bifurcated.intersections.triangular import TriangularIntersection as TI
from multiplex.bifurcated.intersections.digitized import DigitizedIntersection as DI
from multiplex.bifurcated.intersections.curved import CurvedIntersection as CI
from multiplex.meshelement import MeshElement

class BifurcatedManifold(MeshElement):
    """
    Bifurcated Manifold class which sets up a bifurcated manifold.
    You can initalized with:

        BifurcatedManifold()
        BifurcatedManifold(path="working_directory")
        BifurcatedManifold(name="unique_name")
        BifurcatedManifold(path="working_directory",name="unique_name")
        BifurcatedManifold(path="working_directory",name="unique_name",verbose=False)

    """

    def __init__(self,cwd=getcwd(),name=datetime.now().strftime("%d_%m_%Y_%H_%M_%S"),verbose=False):
        super().__init__()
        self.path = cwd
        self.name = name
        self.verbose = verbose
        self.initialized       = False  # bool if the manifold is able to run
        self.started_running   = False  # bool if the manifold has been started to run

        # Create path to working directory if it not exists
        # Load data from YAML file if it exists
        if isfile(join(self.path,self.name)+'.yaml'):
            # Read from file
            if verbose==True:
                print(join(self.path,self.name)+'.yaml exits!')
            stream = open(join(self.path,self.name)+'.yaml', 'r')
            self.meta_data = yaml.load(stream, Loader=yaml.FullLoader)
            self.__initialize_domain__()
            self.__update_meta_data_to_file__()
            self.initialized = True
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
        Load meta data of manifold from YAML file.
        """
        if isfile(join(self.path,self.name)+'.xml'):
            raise Exception("XML file already exits with this name in the woring directory.")
        stream = open(yamlFile, 'r')
        self.meta_data = yaml.load(stream, Loader=yaml.FullLoader)
        self.__initialize_domain__()
        self.__update_meta_data_to_file__()
        self.initialized = True


    def change(self,item,new_value):
        """
        Change a single item and rewrite the data to file
        """
        if self.initialized==False:
            raise Exception("Please initialize all parameters with load() before changing a parameter.")

        if item not in self.meta_data.keys():
            raise Exception("Item to change not found in object meta data. Please check the spelling in the YAML file")

        self.meta_data[item] = new_value
        self.__initialize_domain__()
        self.__update_meta_data_to_file__()


    def __initialize_domain__(self):

        meta_data    = self.meta_data
        width        = meta_data['inlet_width']
        d            = meta_data['inlet_length']
        ddd          = meta_data['device_to_device_distance']
        layers       = meta_data['number_of_layers']
        scale        = meta_data['scale']
        gamma        = meta_data['gamma']
        origin       = meta_data['origin']
        mtype        = meta_data['mesh_type']
        res          = meta_data['mesh_resolution']
        snaps        = meta_data['snaps']
        iterations   = meta_data['iterations']
        dt           = meta_data['time_step']
        mu           = meta_data['viscosity']
        rho          = meta_data['mass_density']
        Re           = meta_data['reynolds_number']

        # Create domain
        self.xval,self.yval,self.interlist=self.__getxy__(layers,ddd*scale,d*scale,width*scale,gamma,origin)
        domain = self.interlist[0].domain
        for k in range(1, len(self.interlist)):
            domain = domain + self.interlist[k].domain

        self.domain = domain
        if isfile(join(self.path,self.name)+'.xml'):
            self.mesh = Mesh(join(self.path,self.name)+'.xml')
            self.mesh_isgenerated = True;
            if self.verbose:
                print("")
                print("Mesh loaded from file: ", join(self.path,self.name)+'.xml')
        else:
            self.generate_mesh(res)
            File(join(self.path,self.name)+'.xml') << self.mesh
            if self.verbose:
                print("")
                print("Mesh newly generated!")


    def __update_meta_data_to_file__(self):

        # Print data to screen
        if self.verbose==True:
            meta_data    = self.meta_data
            width        = meta_data['inlet_width']
            d            = meta_data['inlet_length']
            ddd          = meta_data['device_to_device_distance']
            layers       = meta_data['number_of_layers']
            scale        = meta_data['scale']
            gamma        = meta_data['gamma']
            origin       = meta_data['origin']
            mtype        = meta_data['mesh_type']
            res          = meta_data['mesh_resolution']
            snaps        = meta_data['snaps']
            iterations   = meta_data['iterations']
            dt           = meta_data['time_step']
            mu           = meta_data['viscosity']
            rho          = meta_data['mass_density']
            Re           = meta_data['reynolds_number']

            print("")
            print("inlet_width (m):\t\t\t",width)
            print("inlet_length (m):\t\t\t",d)
            print("device_to_device_distance (m):\t\t",ddd)
            print("layers:\t\t\t\t\t",layers)
            print("scale: \t\t\t\t\t",scale)
            print("")
            print("origin:\t\t\t\t\t", origin)
            print("gamma:\t\t\t\t\t", gamma)
            print("")
            print("mesh_type:\t\t\t\t",mtype)
            print("mesh_resolution:\t\t\t",res)
            print("")
            print("sanps:\t\t\t\t\t", snaps)
            print("iterations:\t\t\t\t",  iterations)
            print("time_step (seconds):\t\t\t", dt)
            print("viscosity (Pa):\t\t\t\t", mu)
            print("mass density (kg/cubic meter):\t\t", rho)
            print("Reynolds number:\t\t\t", Re)

        # Write data to file
        with open(join(self.path,self.name)+'.yaml', 'w') as file:
            documents = yaml.dump(self.meta_data, file)

    def __getxy__(self,layers,ddd,d,width,gamma,origin=[0,0]):
        xval = [origin[0]]
        yval = [origin[1]]
        interlist = []
        ycor = 2**(layers-1)*ddd
        mesh_type = self.meta_data['mesh_type']

        if mesh_type=="rectangular":
            interlist.append(RI([0,0],li=d/2,lo=d/2,lm=ycor+width*gamma,wi=width,wo=width*gamma,wm=width*gamma))
        elif mesh_type=="triangular":
            interlist.append(TI([0,0],li=d/2,lo=d/2,lm=ycor+width*gamma,wi=width,wo=width*gamma,wm=width*gamma))
        elif mesh_type=="digitized":
            interlist.append(DI([0,0],li=d/2,lo=d/2,lm=ycor+width*gamma,wi=width,wo=width*gamma,wm=width*gamma))
        elif mesh_type=="curved":
            interlist.append(CI([0,0],li=d/2,lo=d/2,lm=ycor+width*gamma,wi=width,wo=width*gamma,wm=width*gamma))
        else:
            raise("Interection type not found")

        def __recursiveformula__(n, x, y, ddd, d, width, gamma):

            if n == layers:
                return None

            ycor = 2**(layers-n-1)*ddd

            xval.append(x+d)
            xval.append(x+d)
            yval.append(y+ycor)
            yval.append(y-ycor)

            if mesh_type=="rectangular":
                new_inter = RI([x+d,y+ycor],li=d/2,lo=d/2,lm=ycor+width*gamma,wi=width,wo=width*gamma,wm=width*gamma)
                interlist.append(new_inter)
                new_inter = RI([x+d,y-ycor],li=d/2,lo=d/2,lm=ycor+width*gamma,wi=width,wo=width*gamma,wm=width*gamma)
                interlist.append(new_inter)
            elif mesh_type=="triangular":
                new_inter = TI([x+d,y+ycor],li=d/2,lo=d/2,lm=ycor+width*gamma,wi=width,wo=width*gamma,wm=width*gamma)
                interlist.append(new_inter)
                new_inter = TI([x+d,y-ycor],li=d/2,lo=d/2,lm=ycor+width*gamma,wi=width,wo=width*gamma,wm=width*gamma)
                interlist.append(new_inter)
            elif mesh_type=="digitized":
                new_inter = DI([x+d,y+ycor],li=d/2,lo=d/2,lm=ycor+width*gamma,wi=width,wo=width*gamma,wm=width*gamma)
                interlist.append(new_inter)
                new_inter = DI([x+d,y-ycor],li=d/2,lo=d/2,lm=ycor+width*gamma,wi=width,wo=width*gamma,wm=width*gamma)
                interlist.append(new_inter)
            elif mesh_type=="curved":
                new_inter = CI([x+d,y+ycor],li=d/2,lo=d/2,lm=ycor+width*gamma,wi=width,wo=width*gamma,wm=width*gamma)
                interlist.append(new_inter)
                new_inter = CI([x+d,y-ycor],li=d/2,lo=d/2,lm=ycor+width*gamma,wi=width,wo=width*gamma,wm=width*gamma)
                interlist.append(new_inter)

            __recursiveformula__(n+1,x+d,y+ycor, ddd, d, width*gamma, gamma)
            __recursiveformula__(n+1,x+d,y-ycor, ddd, d, width*gamma, gamma)

        __recursiveformula__(1, 0, 0, ddd, d, width*gamma, gamma)
        return np.array(xval), np.array(yval), interlist

    def walls(self,x):
        """
        Returns True if x is a point on the walls
        """
        if any(inter.walls(x) == True for inter in self.interlist):
            return True
        else:
            return False

    def inflow(self, x):
        """
        Returns True if x is a point on the inlet
        """

        return self.interlist[0].inflow(x)

    def outflow(self, x):
        """
        Returns True if x is a point on the outlet
        """

        outflow_inter = []
        for inter in self.interlist:
            if inter.origin[0] == np.max(self.xval):
                outflow_inter.append(inter)

        if any(inter.outflow(x) == True for inter in outflow_inter):
            return True
        else:
            return False

    def outlets(self):
        """
        Solving time dependent Navier Stokes equation for a given number of snap shots 
        """

        outflow_inter = []
        for inter in self.interlist:
            if inter.origin[0] == np.max(self.xval):
                outflow_inter.append(inter)

        outlets = []

        for x, [lower, upper] in [inter.outlets() for inter in outflow_inter]:
            outlets.append(lower)
            outlets.append(upper)

        return x, outlets

    def solve(self,snaps=None,overwrite=False):
        """
        Solving time dependent Navier Stokes equation for a given number of snap shots 
        """

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
        type   = meta_data['mesh_type']
        res    = int(meta_data['mesh_resolution'])
        snapsF = int(meta_data['snaps'])
        iter   = int(meta_data['iterations'])
        dt     = float(meta_data['time_step'])
        mu     = float(meta_data['viscosity'])
        rho    = float(meta_data['mass_density'])
        Re     = float(meta_data['reynolds_number'])

        if(snaps != None):
            if(snapsF < snaps):
                self.meta_data['snaps'] = snaps
                self.__update_meta_data_to_file__()
        snaps = self.meta_data['snaps']

        v_mean = Re*mu/rho/width
        inflow_profile=(str(v_mean*scale*3/2/(width*scale/2)**2)+'*(x[1]+'+str(width*scale/2)+')*('+str(width*scale/2)+'-x[1])', '0')

        self.solver(self.path,self.name,snaps,iter,mu/scale,rho/scale**3,dt,inflow_profile,overwrite,verbose=self.verbose,consecutive_run_allowed=True)

    def hdf2pvd(self):
        """
        Converts HDF files in the working folder to PVD files
        """

        if self.initialized == False:
            raise RuntimeError("Manifold is not initialized!")

        destination = join(self.path,self.name)
        source      = self.path

        if not isdir(destination):
            mkdir(destination)

        self.hdf_to_pvd(destination,source,self.name,verbose=self.verbose)


    def get_Q(self, series=False):
        """
        Returns an array of the exit flow rate at each exit
        """

        if self.initialized == False:
            raise RuntimeError("Manifold is not initialized!")

        V = VectorFunctionSpace(self.mesh, 'P', 2)
        Q = FunctionSpace(self.mesh, 'P', 1)

        u_  = Function(V)
        p_  = Function(Q)

        hdf = HDF5File(self.mesh.mpi_comm(),join(self.path,self.name)+'.h5', "r")
        attr = hdf.attributes("p")
        nsteps = attr['count']

        start = nsteps-1
        if series: start = 0
        x, outlets = self.outlets()
        x_cor = np.linspace(x ,x ,100)[1:-1]

        series = []
        for step in range(start, nsteps):

            averages = []

            for [y1, y2] in outlets:

                dataset = "u/vector_%d"%step
                hdf.read(u_,dataset)
                u = interpolate(u_, V)

                y_cor = np.linspace(y2,y1,100)[1:-1]
                velocities = []
                for x,y in zip(x_cor,y_cor):
                    velocities.append(u(x,y))
                velocities = np.array(velocities)

                scale = self.meta_data['scale']
                averages.append(np.trapz(velocities[:,0],-1*y_cor)/scale/scale)

            series.append(averages) #Divide by scale^2 (because of 2D)

        hdf.close()
        return np.array(series)

    def get_Qin(self):
        """
        Returns the inlet flow rate rate of the bifurcated manifold
        """

        if(self.initialized == False):
            raise RuntimeError("Manifold not initialized!")

        width           = self.meta_data['inlet_width']
        reynolds_number = self.meta_data['reynolds_number']
        mass_density    = self.meta_data['mass_density']
        viscosity       = self.meta_data['viscosity']
        v_mean = Re*mu/rho/width
        return np.array([width*v_mean])
