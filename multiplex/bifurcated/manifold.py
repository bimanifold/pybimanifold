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
    Bifurcated Manifold class 

    Initialziation paramters:
    cwd (str):       current working directory
    name (str):      unique name of the manifold (default: "day_month_year_hour_minute_second")
    verbose (bool):  prints information to stdout

    Examples:
        mymanifold = BifurcatedManifold()
        mymanifold = BifurcatedManifold(cwd="working_directory")
        mymanifold = BifurcatedManifold(cwd="unique_name")
        mymanifold = BifurcatedManifold(cwd="working_directory",name="unique_name")
        mymanifold = BifurcatedManifold(cwd="working_directory",name="unique_name",verbose=False)

    """

    def __init__(self,cwd,name=datetime.now().strftime("%d_%m_%Y_%H_%M_%S"),verbose=False):
        super().__init__(cwd=cwd,name=name,verbose=verbose)


    def __initialize_domain__(self):

        width        = self.meta_data['inlet_width']
        d            = self.meta_data['inlet_length']
        ddd          = self.meta_data['device_to_device_distance']
        layers       = self.meta_data['number_of_layers']
        scale        = self.meta_data['scale']
        gamma        = self.meta_data['gamma']
        origin       = self.meta_data['origin']
        mtype        = self.meta_data['mesh_type']
        res          = self.meta_data['mesh_resolution']

        # Create domain
        self.xval,self.yval,self.interlist=self.__getxy__(layers,ddd*scale,d*scale,width*scale,gamma,origin)
        domain = self.interlist[0].domain
        for k in range(1, len(self.interlist)):
            domain = domain + self.interlist[k].domain

        self.domain = domain


    def __print_domain_meta_data_to_screen__(self):

        # Print data to screen
        if self.verbose==True:
            width        = self.meta_data['inlet_width']
            d            = self.meta_data['inlet_length']
            ddd          = self.meta_data['device_to_device_distance']
            layers       = self.meta_data['number_of_layers']
            scale        = self.meta_data['scale']
            gamma        = self.meta_data['gamma']
            origin       = self.meta_data['origin']
            mtype        = self.meta_data['mesh_type']

            print("")
            print("inlet_width (m):\t\t\t",width)
            print("inlet_length (m):\t\t\t",d)
            print("device_to_device_distance (m):\t\t",ddd)
            print("layers:\t\t\t\t\t",layers)
            print("scale: \t\t\t\t\t",scale)
            print("origin:\t\t\t\t\t", origin)
            print("gamma:\t\t\t\t\t", gamma)
            print("mesh_type:\t\t\t\t",mtype)
            print("")


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
            raise Exception("Interection type not found")

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

        Prameters:
        x (dolfin.Point):   point in the x-y plane

        Example:
        >>> import dolfin
        >>> mymanifold.BifurcatedManifold()
        >>> mymanifold.load("sample.yaml")
        >>> x = dolfin.Point(4,5)
        >>> mymanifold.walls(x)
        False

        """

        if any(inter.walls(x) == True for inter in self.interlist):
            return True
        else:
            return False


    def inflow(self, x):
        """
        Returns True if x is a point on the inlet

        Prameters:
        x (dolfin.Point):   point in the x-y plane

        Example:
        >>> import dolfin
        >>> mymanifold.BifurcatedManifold()
        >>> mymanifold.load("sample.yaml")
        >>> x = dolfin.Point(4,5)
        >>> mymanifold.walls(x)
        False

        """

        return self.interlist[0].inflow(x)

    def outflow(self, x):
        """
        Returns True if x is a point on the outlet

        Prameters:
        x (dolfin.Point):   point in the x-y plane

        Example:
        >>> import dolfin
        >>> mymanifold.BifurcatedManifold()
        >>> mymanifold.load("sample.yaml")
        >>> x = dolfin.Point(4,5)
        >>> mymanifold.walls(x)
        False

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


    def get_Q(self, series=False):
        """
        Returns a numpy array of the exit flow rate at each exit.

        Parameters:
        series (bool):  if true return the flow rate for each snap shot

        Return:
        Q (np.array):   flow rate for each outlet

        Examples:
        >>> Qflows = mymanifold.get_Q()
        >>> Qflows = mymanifold.get_Q(series=True)

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
                averages.append(np.trapz(velocities[:,0],-1*y_cor)/scale/scale) #Divide by scale^2 (because of 2D)

            series.append(averages) 

        hdf.close()
        return np.array(series)



    def get_Qin(self):
        """
        Returns the inlet flow rate rate of the bifurcated manifold

        Parameters:
        None

        Returns:
        Qin (int):  inlet flow rate

        Example:
        >> Qin = mymanifold.get_Qin()

        """
        
        if(self.initialized == False):
            raise RuntimeError("Manifold not initialized!")

        width           = self.meta_data['inlet_width']
        reynolds_number = self.meta_data['reynolds_number']
        mass_density    = self.meta_data['mass_density']
        viscosity       = self.meta_data['viscosity']
        Re              = self.meta_data['reynolds_number']
        v_mean = Re*viscosity/mass_density/width
        return width*v_mean
