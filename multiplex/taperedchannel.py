import numpy as np
from dolfin import *
from mshr import *
from multiplex.meshelement import MeshElement


class TaperedChannel(MeshElement):

    def __init__(self,channeltype,startwidth,endwith,length, lengthtapered, radius=1):
        super().__init__()

        self.h1 = startwidth
        self.h2 = endwith
        self.L  = lengthtapered/2
        self.R  = radius
        self.type = channeltype
        self.length = length/2
        self.mesh

        b = np.array([         self.h2,           0,   self.h1,     0])
        A = np.array([[      self.L**3,   self.L**2,    self.L,     1],
                      [    3*self.L**2,    2*self.L,         1,     0],
                      [   (-self.L)**3,(-self.L)**2,   -self.L,     1],
                      [ 3*(-self.L)**2, 2*(-self.L),         1,     0]])
        self.coefspline1 = np.linalg.solve(A,b)

        b = np.array([self.h2,0,0,self.h1,0,0])
        A = np.array([[      self.L**5,       self.L**4,      self.L**3,    self.L**2,   self.L,  1],
                      [    5*self.L**4,     4*self.L**3,    3*self.L**2,     2*self.L,        1,  0],
                      [   20*self.L**3,    12*self.L**2,    6*self.L**1,            2,        0,  0],
                      [   (-self.L)**5,    (-self.L)**4,   (-self.L)**3, (-self.L)**2,  -self.L,  1],
                      [ 5*(-self.L)**4,  4*(-self.L)**3, 3*(-self.L)**2,  2*(-self.L),        1,  0],
                      [20*(-self.L)**3, 12*(-self.L)**2, 6*(-self.L)**1,            2,        0,  0]])
        self.coefspline2 = np.linalg.solve(A,b)


    def spline0(self,x):
        if x < -self.L:
            return self.h1
        if x < self.L:
            gamma = np.arctan((self.h2-self.h1)/2/self.L)
            return self.h1 +(self.h2-self.h1)/2+x*np.tan(gamma)
        return self.h2

    def curvspline0(self,x):
        R        = 1
        gamma    = np.arctan((self.h2-self.h1)/2/self.L)
        alpha    = (np.pi - gamma)/2
        b        = R/np.tan(alpha)
        if x < -b-self.L:
            return self.h1
        if x < -self.L+b*np.cos(gamma):
            return -np.sqrt(R**2-(x+self.L+b)**2)+self.h1+R
        if x < self.L-b*np.cos(gamma):
            gamma = np.arctan((self.h2-self.h1)/2/self.L)
            return self.h1 +(self.h2-self.h1)/2+x*np.tan(gamma)
        if x < self.L+b:
            return  np.sqrt(R**2-(x-self.L-b)**2)+self.h2-R
        return self.h2

    def spline1(self,x):
        if x < -self.L:
            return self.h1
        if x < self.L:
            return np.dot(self.coefspline1, np.array([x**3,x**2,x,1]))
        return self.h2

    def spline2(self,x):
        if x < -self.L:
            return self.h1
        if x < self.L:
            return np.dot(self.coefspline2, np.array([x**5,x**4,x**3,x**2,x,1]))
        return self.h2

    def gety(self,x):
        if self.type == 'spline0':
            try:
                return np.array([self.spline0(xval) for xval in x])
            except TypeError:
                return self.spline0(x)
        if self.type == 'spline0curv':
            try:
                return np.array([self.curvspline0(xval) for xval in x])
            except TypeError:
                return self.curvspline0(x)
        if self.type == 'spline1':
            try:
                return np.array([self.spline1(xval) for xval in x])
            except TypeError:
                return self.spline1(x)
        if self.type == 'spline2':
            try:
                return np.array([self.spline2(xval) for xval in x])
            except TypeError:
                return self.spline2(x)

    def getxy(self, resolution):
        if self.type == 'curvspline0':
            R        = 1
            gamma    = np.arctan((self.h2-self.h1)/2/self.L)
            alpha    = (np.pi - gamma)/2
            b        = R/np.tan(alpha)
            x = np.linspace(-self.L-b,self.L+b,resolution)
            y = self.gety(x)
            x = np.insert(x, 0, -self.length)
            y = np.insert(y, 0, self.h1)
            x = np.append(x, self.length)
            y = np.append(y, self.h2)
        else:
            x = np.linspace(-self.L,self.L,resolution)
            y = self.gety(x)
            x = np.insert(x, 0, -self.length)
            y = np.insert(y, 0, self.h1)
            x = np.append(x, self.length)
            y = np.append(y, self.h2)
        return x,y

    def domain(self, resolution):
        x, y = self.getxy(2*resolution)
        domain_vertices = []
        for xval,yval in zip(x,y):
            domain_vertices.append(Point(xval,-yval))
        for xval,yval in zip(np.flip(x, axis=0),np.flip(y, axis=0)):
            domain_vertices.append(Point(xval,yval))
        return Polygon(domain_vertices)

    def inflow(self, x):
        return near(x[0], -self.length)

    def outflow(self, x):
        return near(x[0],  self.length)

    def walls(self, x):
        return near(x[1], self.gety(x[0])) or near(x[1], -self.gety(x[0]))

    def generate_mesh(self, resolution, meshtype='unstructured'):
        if meshtype=='unstructured':
            self.mesh = generate_mesh(self.domain(resolution), resolution)
        if meshtype=='structured':
            self.mesh = RectangleMesh(Point(-self.length,-self.h2),Point(self.length,self.h2), resolution[0],resolution[1])
            x = self.mesh.coordinates()
            x[:,1] = np.abs(self.gety(x[:,0]) / np.max(x[:,1])) * x[:,1]
        self.mesh_isgenerated = True

    def getxcor(self):
        x = self.mesh.coordinates()
        return np.unique(x[:,0])
