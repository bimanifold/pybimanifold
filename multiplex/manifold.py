from dolfin import *
from mshr import *
import numpy as np
from os.path import isfile, isdir, split
from os import mkdir, makedirs
from re import findall
from multiplex.intersection_rectangular import RectangularIntersection as RI
from multiplex.intersection_triangular import TriangularIntersection as TI
from multiplex.intersection_digitized import DigitizedIntersection as DI
from multiplex.intersection_curved import CurvedIntersection as CI
from multiplex.meshelement import MeshElement

class Manifold(MeshElement):

    def __init__(self,path,name,geoparams=[],solverparams=[],gamma=0.79370052598,origin=[0,0],verbose=True):
        super().__init__()
        self.path = path
        self.name = name
        self.verbose = verbose

        # Check if file already exists
        new_analysis = False
        if isfile(path+name+".txt"):

            # Read from file
            if verbose==True:
                print(path+name+".txt"+" exits!")
            file1 = open(path+name+".txt","r")
            f1 = file1.readlines()

            def get_number(string):
                temp=findall(r'[+-]? *(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?', string)
                if temp != []:
                    temp=temp[0]
                    temp=float(temp)
                    if temp==int(temp):
                        temp = int(temp)
                    return temp

            width  = get_number(f1[0])
            d      = get_number(f1[1])
            ddd    = get_number(f1[2])
            layers = get_number(f1[3])
            scale  = get_number(f1[4])
            res    = get_number(f1[6])
            cells  = get_number(f1[7])
            type   = f1[8].split('\n')[0].split('\t')[-1]
            snaps  = get_number(f1[10])
            iter   = get_number(f1[11])
            dt     = get_number(f1[12])
            mu     = get_number(f1[13])
            rho    = get_number(f1[14])
            Re     = get_number(f1[15])
            if(snaps < solverparams[0]):
                snaps = solverparams[0]
                new_analysis = True
            solverparams = [snaps,iter,dt,mu,rho,Re]
            geoparams    = [type,res,scale,layers,ddd,d,width]
            self.num_cells = cells

        else:

            if isdir(path):
                if verbose==True:
                    print(path+" exits!")
            else:
                makedirs(path)
                if verbose==True:
                    print(path+" created!")
            [type,res,scale,layers,ddd,d,width] = geoparams
            new_analysis = True

        self.solverparams = solverparams
        self.geoparams    = geoparams
        [snaps,iter,dt,mu,rho,Re] = self.solverparams
        self.type = type
        self.resolution = res


        # Create domain
        self.xval,self.yval,self.interlist=self.getxy(layers,ddd*scale,d*scale,width*scale,gamma,origin)
        domain = self.interlist[0].domain
        for k in range(1, len(self.interlist)):
            domain = domain + self.interlist[k].domain

        self.domain = domain
        if isfile(path+name+'.xml'):
            self.mesh = Mesh(path+name+'.xml')
            print("Mesh loaded from file: ", path+name+'.xml')
        else:
            self.generate_mesh(self.resolution)
            File(path+name+'.xml') << self.mesh
            if self.verbose:
                print("Mesh newly generated!")
            if (new_analysis==False):
                raise Exception("MESH ERROR: could not load previous mesh")

        self.v_mean = Re*mu/rho/width
        self.inflow_profile=(str(self.v_mean*scale*3/2/(width*scale/2)**2)+'*(x[1]+'+str(width*scale/2)+')*('+str(width*scale/2)+'-x[1])', '0')


        # Print data to screen
        if verbose==True:
            print("")
            print("wi (m):\t\t\t",width)
            print("d (m):\t\t\t",d)
            print("ddd (m):\t\t",ddd)
            print("layers:\t\t\t",layers)
            print("scale:\t\t\t",scale)
            print("")
            print("mesh res:\t\t",res)
            print("num cells:\t\t",self.mesh.num_cells())
            print("type:\t\t\t",type)
            print("")
            print("num_snaps:\t\t", snaps)
            print("num_iter:\t\t",  iter)
            print("dt (seconds):\t\t", dt)
            print("mu (Pa):\t\t", mu)
            print("rho (kg/cubic meter):\t", rho)
            print("Reynolds number:\t", Re)
            print("mean velocity (m/s):\t", self.v_mean)
            print("")
            print("inflow_profile:\t", self.inflow_profile)
            print("")

        # Write data to file
        if new_analysis==True:
            file1 = open(path+name+".txt","a+")
            file1.write("wi (m):\t\t\t"+str(width)+"\n")
            file1.write("d (m):\t\t\t"+str(d)+"\n")
            file1.write("ddd (m):\t\t"+str(ddd)+"\n")
            file1.write("layers:\t\t\t"+str(layers)+"\n")
            file1.write("scale:\t\t\t"+str(scale)+"\n")
            file1.write("\n")
            file1.write("mesh res:\t\t"+str(res)+"\n")
            file1.write("num cells:\t\t"+str(self.mesh.num_cells())+"\n")
            file1.write("type:\t\t\t"+str(type)+"\n")
            file1.write("\n")
            file1.write("num_snaps:\t\t"+str(snaps)+"\n")
            file1.write("num_iter:\t\t"+str(iter)+"\n")
            file1.write("dt (seconds):\t\t"+str(dt)+"\n")
            file1.write("mu (Pa):\t\t"+str(mu)+"\n")
            file1.write("rho (kg/cubic meter):\t"+str(rho)+"\n")
            file1.write("Reynolds number:\t"+str(Re)+"\n")
            file1.write("mean velocity (m/s):\t"+str(self.v_mean)+"\n")
            file1.write("\n")
            file1.write("inflow_profile:\t"+str(self.inflow_profile)+"\n")
            file1.write("\n")
            file1.close()

    def getxy(self,layers,ddd,d,width,gamma,origin=[0,0]):
        xval = [origin[0]]
        yval = [origin[1]]
        interlist = []
        ycor = 2**(layers-1)*ddd

        if self.type=="rectangular":
            interlist.append(RI([0,0],li=d/2,lo=d/2,lm=ycor+width*gamma,wi=width,wo=width*gamma,wm=width*gamma))
        elif self.type=="triangular":
            interlist.append(TI([0,0],li=d/2,lo=d/2,lm=ycor+width*gamma,wi=width,wo=width*gamma,wm=width*gamma))
        elif self.type=="digitized":
            interlist.append(DI([0,0],li=d/2,lo=d/2,lm=ycor+width*gamma,wi=width,wo=width*gamma,wm=width*gamma))
        elif self.type=="curved":
            interlist.append(CI([0,0],li=d/2,lo=d/2,lm=ycor+width*gamma,wi=width,wo=width*gamma,wm=width*gamma))
        else:
            raise("Interection type not found")

        def recursiveformula(n, x, y, ddd, d, width, gamma):

            if n == layers:
                return None

            ycor = 2**(layers-n-1)*ddd

            xval.append(x+d)
            xval.append(x+d)
            yval.append(y+ycor)
            yval.append(y-ycor)

            if self.type=="rectangular":
                new_inter = RI([x+d,y+ycor],li=d/2,lo=d/2,lm=ycor+width*gamma,wi=width,wo=width*gamma,wm=width*gamma)
                interlist.append(new_inter)
                new_inter = RI([x+d,y-ycor],li=d/2,lo=d/2,lm=ycor+width*gamma,wi=width,wo=width*gamma,wm=width*gamma)
                interlist.append(new_inter)
            elif self.type=="triangular":
                new_inter = TI([x+d,y+ycor],li=d/2,lo=d/2,lm=ycor+width*gamma,wi=width,wo=width*gamma,wm=width*gamma)
                interlist.append(new_inter)
                new_inter = TI([x+d,y-ycor],li=d/2,lo=d/2,lm=ycor+width*gamma,wi=width,wo=width*gamma,wm=width*gamma)
                interlist.append(new_inter)
            elif self.type=="digitized":
                new_inter = DI([x+d,y+ycor],li=d/2,lo=d/2,lm=ycor+width*gamma,wi=width,wo=width*gamma,wm=width*gamma)
                interlist.append(new_inter)
                new_inter = DI([x+d,y-ycor],li=d/2,lo=d/2,lm=ycor+width*gamma,wi=width,wo=width*gamma,wm=width*gamma)
                interlist.append(new_inter)
            elif self.type=="curved":
                new_inter = CI([x+d,y+ycor],li=d/2,lo=d/2,lm=ycor+width*gamma,wi=width,wo=width*gamma,wm=width*gamma)
                interlist.append(new_inter)
                new_inter = CI([x+d,y-ycor],li=d/2,lo=d/2,lm=ycor+width*gamma,wi=width,wo=width*gamma,wm=width*gamma)
                interlist.append(new_inter)

            recursiveformula(n+1,x+d,y+ycor, ddd, d, width*gamma, gamma)
            recursiveformula(n+1,x+d,y-ycor, ddd, d, width*gamma, gamma)

        recursiveformula(1, 0, 0, ddd, d, width*gamma, gamma)
        return np.array(xval), np.array(yval), interlist

    def walls(self,x):
        if any(inter.walls(x) == True for inter in self.interlist):
            return True
        else:
            return False

    def inflow(self, x):
        return self.interlist[0].inflow(x)

    def outflow(self, x):

        outflow_inter = []
        for inter in self.interlist:
            if inter.origin[0] == np.max(self.xval):
                outflow_inter.append(inter)

        if any(inter.outflow(x) == True for inter in outflow_inter):
            return True
        else:
            return False

    def outlets(self):

        outflow_inter = []
        for inter in self.interlist:
            if inter.origin[0] == np.max(self.xval):
                outflow_inter.append(inter)

        outlets = []

        for x, [lower, upper] in [inter.outlets() for inter in outflow_inter]:
            outlets.append(lower)
            outlets.append(upper)

        return x, outlets

    def solve(self,overwrite=False):
        [snaps,iter,dt,mu,rho,Re]           = self.solverparams
        [type,res,scale,layers,ddd,d,width] = self.geoparams
        self.solver(self.path,self.name, snaps, iter, mu/scale, rho/scale**3, dt, self.inflow_profile, overwrite, verbose=self.verbose, consecutive_run_allowed=True)

    def hdf2pvd(self):
        destination = self.path+self.name+'/'
        source      = self.path

        if not isdir(destination):
            mkdir(destination)

        self.hdf_to_pvd(destination,source,self.name,verbose=self.verbose)


    def get_Q(self, series=False):

        V = VectorFunctionSpace(self.mesh, 'P', 2)
        Q = FunctionSpace(self.mesh, 'P', 1)

        u_  = Function(V)
        p_  = Function(Q)

        hdf = HDF5File(self.mesh.mpi_comm(),self.path+self.name+".h5", "r")
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

                averages.append(np.trapz(velocities[:,0],-1*y_cor)/self.geoparams[2]/self.geoparams[2])

            series.append(averages) #Divide by scale^2 (because of 2D)

        hdf.close()

        return np.array(series)

    def get_Qin(self):
        return np.array([self.geoparams[-1]*self.v_mean])
