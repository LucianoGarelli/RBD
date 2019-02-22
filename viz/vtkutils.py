## Import some modules
import vtk, sys
from numpy import *
from math import pi,cos,sin,tan,atan
from distutils.version import LooseVersion, StrictVersion
DEG = pi/180.0

def vtkvers():
    vers = vtk.vtkVersion().GetVTKVersion()
    return LooseVersion(vers) >= LooseVersion("6.0.0")

def setinputdata_wrp(vtkobject,data):
    if vtkvers():
        vtkobject.SetInputData(data)
    else:
        vtkobject.SetInput(data)

def setinputconn_wrp(vtkobject,data):
    if vtkvers():
        vtkobject.SetInputConnection(data)
    else:
        vtkobject.SetInput(data)

def setsourcedata_wrp(vtkobject,data):
    if vtkvers():
        vtkobject.SetSourceData(data)
    else:
        vtkobject.SetSource(data)

def setsourceconn_wrp(vtkobject,source):
    if vtkvers():
        vtkobject.SetSourceConnection(source.GetOutputPort())
    else:
        vtkobject.SetSource(source.GetOutput())

## Set vtk `nodes' from coords array x
def make_nodes(x,nodes):
    ids = vtk.vtkIdList()                   # ids of particles
    nnod = x.shape[0]
    ndim = x.shape[1]
    xx = zeros(3)
    for k in range(0,nnod):
        if ndim<3:
            xx[0:2] = x[k]
        else:
            xx = x[k]
        nodes.InsertNextPoint(xx)
        ids.InsertNextId(k)

## Set the connectivity of a nesh from an integer array
## `conec' and a typical element `elem'
def set_mesh_connectivity(mesh,conec,elem):
    nelem = conec.shape[0]
    nel = conec.shape[1]
    mesh.Allocate(nelem,0)
    for k in range(0,nelem):
        tids = elem.GetPointIds()
        # tids = vtk.vtkIdList()
        for l in range(0,nel):
            node = int(conec[k,l])
            tids.SetId(l,node)
        mesh.InsertNextCell(elem.GetCellType(),tids)

## Make a FEM mesh grid from cords and connectivities
def make_actor(xnod,icone,elem,**kwargs):
    color = kwargs.get("color",[0,1,0])
    opacity = kwargs.get("opacity",1.0)
    nodes = vtk.vtkPoints()
    make_nodes(xnod,nodes)
    nnod = xnod.shape[0]
    ndim = xnod.shape[1]
    nel = icone.shape[1]

    grid = vtk.vtkPolyData()
    grid.SetPoints(nodes)

    ## if elem is None: Guess element type from nel and dim
    set_mesh_connectivity(grid,icone,elem)

    # Create mapper
    mapper = vtk.vtkDataSetMapper()
    setinputdata_wrp(mapper,grid)

    # Create actor
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    prop = actor.GetProperty()
    prop.SetDiffuseColor(color)
    prop.SetOpacity(opacity)
    return (actor,mapper,dict(grid=grid,nodes=nodes))

## Make a tube from coordinates of the center of the tubexs
def make_tube(xnod,**kwargs):
    color = kwargs.get("color",[0,1,0])
    radius = kwargs.get("radius",0.1)
    sides = kwargs.get("sides",6)
    
    nodes = vtk.vtkPoints()
    make_nodes(xnod,nodes)
    nnod = xnod.shape[0]
    ndim = xnod.shape[1]

    icone = zeros((1,nnod))
    for l in range(0,nnod):
        icone[0,l] = l
    
    grid = vtk.vtkPolyData()
    grid.SetPoints(nodes)

    elem = vtk.vtkPolyLine()
    elem.GetPointIds().SetNumberOfIds(nnod)
    set_mesh_connectivity(grid,icone,elem)

    tuber = vtk.vtkTubeFilter()
    print radius
    tuber.SetRadius(radius)
    tuber.SetNumberOfSides(sides)
    if vtkvers():
        tuber.SetInputData(grid)
    else:
        tuber.SetInput(grid)

    # Create mapper
    mapper = vtk.vtkDataSetMapper()
    mapper.SetInputConnection(tuber.GetOutputPort())

    # Create actor
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetDiffuseColor(color)
    return (actor,mapper)

## Make a patch of quads. The coordinates are defined
## by the mapping function `xfun'.
## usage: (actor,mapper,data)
##        = make_quad_patch(Nxi,Neta,color=,fundata=)
def make_quad_patch(Nxi,Neta,xfun,**kwargs):
    color = kwargs.get("color",[0,1,0])
    fundata = kwargs.get("fundata",None)
    
    Nxi1 = Nxi+1
    Neta1 = Neta+1
    nnod = Nxi1*Neta1
    xnod = zeros((nnod,3))
    xiv = zeros((nnod,2))
    node = 0
    for jxi in range(0,Nxi1):
        xi = float(jxi)/Nxi
        for jeta in range(0,Neta1):
            eta = float(jeta)/Neta
            xnod[node] = xfun(xi,eta,fundata)
            xiv[node] = [xi,eta]
            node += 1

    icone = zeros((Nxi*Neta,4))
    elem = 0
    for jxi in range(0,Nxi):
        for jeta in range(0,Neta):
            base = Neta1*jxi+jeta
            icone[elem] = base + array([0,Neta1,Neta1+1,1])
            elem += 1
    (actor,mapper,data) = make_actor(xnod,icone,vtk.vtkQuad(),**kwargs)
    data['xi'] = xiv
    return (actor,mapper,data)

def draw_axes(xc=array([0,0,0]),L=1,D=0.1):
    cyl = vtk.vtkCylinderSource()
    cyl.SetHeight(L)
    #cyl.SetCenter(xc)
    cyl.SetRadius(D/2.0)
    cyl.SetResolution(30)

    cylMapper = vtk.vtkPolyDataMapper()
    cylMapper.SetInputConnection(cyl.GetOutputPort() )

    ## Create the cylinder actor
    cylY = vtk.vtkActor()
    cylY.SetMapper(cylMapper)
    cylY.SetPosition(xc)

    cylX = vtk.vtkActor()
    cylX.SetMapper(cylMapper)
    cylX.RotateZ(90.0)
    cylX.SetPosition(xc)

    cylZ = vtk.vtkActor()
    cylZ.SetMapper(cylMapper)
    cylZ.RotateX(90.0)
    cylZ.SetPosition(xc)

    cylX.GetProperty().SetColor(1.0,0.0,0.0)
    cylY.GetProperty().SetColor(0.0,1.0,0.0)
    cylZ.GetProperty().SetColor(0.0,0.0,1.0)
    return [cylX,cylY,cylZ]

## Patch of quads, usage:
##    qp = quad_patch(Nxi,Neta,color=
class quad_patch():
    def __init__(self,Nxi,Neta,**kwargs):
        color = kwargs.get("color",[0,1,0])
        self.pre_reset_coords(**kwargs)

        Nxi1 = Nxi+1
        Neta1 = Neta+1
        nnod = Nxi1*Neta1
        xnod = zeros((nnod,3))
        self.xi = zeros((nnod,2))
        self.nnod = nnod
        node = 0
        for jxi in range(0,Nxi1):
            xi = float(jxi)/Nxi
            for jeta in range(0,Neta1):
                eta = float(jeta)/Neta
                xnod[node] = self.xfun(xi,eta)
                self.xi[node] = [xi,eta]
                node += 1

        icone = zeros((Nxi*Neta,4))
        elem = 0
        for jxi in range(0,Nxi):
            for jeta in range(0,Neta):
                base = Neta1*jxi+jeta
                icone[elem] = base + array([0,Neta1,Neta1+1,1])
                elem += 1
        (self.actor,self.mapper,data) \
          = make_actor(xnod,icone,vtk.vtkQuad(),**kwargs)
        self.nodes = data['nodes']
        self.grid = data['grid']

    def reset_coords(self):
        self.pre_reset_coords()
        nnod = self.nnod
        for k in range(0,nnod):
            x = self.xfun(self.xi[k,0],self.xi[k,1])
            self.nodes.SetPoint(k,x[0],x[1],x[2])
        self.nodes.Modified()
          
##   s = patch_t(N=1,color=)
##   s.xcorners = (array of 4x3) corners of the patch (counter cockwise)
class patch_t(quad_patch):
    def __init__(self,xcorners,N=1,**kwargs):
        self.xcorners = xcorners
        quad_patch.__init__(self,N,N,**kwargs)
        
    def pre_reset_coords(self,**kwargs):
        None

    def xfun(self,xi,eta):
        xc = self.xcorners
        x = xc[0]*(1-xi)*(1-eta) \
          + xc[1]*xi*(1-eta) \
          + xc[2]*xi*eta \
          + xc[3]*(1-xi)*eta
        return x

## Sphere object, usage
##   s = sphere(N,color=)
##   s.R = radius,
##   s.xc = sphere center
##   s.reset_coords()
class sphere(quad_patch):
    def __init__(self,N=20,**kwargs):
        quad_patch.__init__(self,N,N,**kwargs)
        
    def pre_reset_coords(self,**kwargs):
        (self.xc,self.R) = self.params()

    def xfun(self,xi,eta):
        phi = (xi-0.5)*2.0*pi
        theta = (eta-0.5)*pi
        rho = cos(theta)
        z = sin(theta)
        x = array([rho*cos(phi),rho*sin(phi),z])
        x = x*self.R + self.xc
        return x

class sphere2_t(sphere):
    def __init__(self,N,**kwargs):
        self.xc = zeros((3))
        self.color = kwargs.get("color",[0,1,0])
        self.R = kwargs.get("R",0.5)
        sphere.__init__(self,N,**kwargs)
        
    def params(self):
        return (self.xc,self.R)

## Vector product of 2 3D vectors
def pvec(x,y):
    z = zeros((3))
    z[0] = x[1]*y[2]-x[2]*y[1]
    z[1] = x[2]*y[0]-x[0]*y[2]
    z[2] = x[0]*y[1]-x[1]*y[0]
    return z

## Norm of a vector
def norm(x):
    sum = 0.0
    n = x.shape[0]
    for k in range(0,n):
        sum += x[k]**2
    return sqrt(sum)

## Normalize a vector, returns a normalized copy
def normalize(x):
    xx = copy(x)
    ax = norm(xx)
    if ax>0.0:
        xx /= ax
    return xx

## Given a 3D vector `x' returns a 3x3 matrix where S[0]
## is `x' and the S[1] S[2] are unit length and orthogonal to `x'
## also S[0],S[1],S[2] form a right-handed basis
def orth_basis(x):
    k = 0
    if abs(x[1])<abs(x[0]):
        k = 1
    if abs(x[2])<abs(x[1]):
        k = 2
    aux = zeros((3))
    aux[k] = 1
    t1 = pvec(x,aux)
    t2 = pvec(x,t1)
    S = zeros((3,3))
    S[0] = normalize(x)
    S[1] = normalize(t1)
    S[2] = normalize(t2)
    return S

## Makes a tube.
## usage: t = tube(Nphi,Nz,color=)
##  self.xends = array 2x3
##  self.R = 
class tube(quad_patch):
    def __init__(self,Nphi=20,Nz=10,**kwargs):
        quad_patch.__init__(self,Nphi,Nz,**kwargs)
        
    def pre_reset_coords(self,**kwargs):
        (self.xends,self.R) = self.params()
        self.x0 = self.xends[0]
        self.dx = self.xends[1]-self.x0
        S = orth_basis(self.dx)
        self.t1 = S[1]
        self.t2 = S[2]
        self.L = norm(self.dx)

    def xfun(self,xi,eta):
        phi = (xi-0.5)*2.0*pi
        x = self.x0 + \
            self.R * cos(phi) * self.t1 + \
            self.R * sin(phi) * self.t2 + \
            eta * self.dx
        return x

## Makes a sector (as a crash dummy Secci disk)
## We use this to mark a lab reference system. 
class cross_mark_sector_t(quad_patch):
    def __init__(self,indx,color,Nphi=5,Nsector=4,**kwargs):
        self.indx = indx
        self.Nsector = 4
        self.xc = array([0,0,0])
        self.t1 = array([1,0,0])
        self.t2 = array([0,1,0])
        self.R = 0.04
        quad_patch.__init__(self,Nphi,1,color=color,**kwargs)

    def pre_reset_coords(self,**kwargs):
        None
        
    def xfun(self,xi,eta):
        dphi = 2*pi/float(self.Nsector)
        phi = (self.indx+xi)*dphi
        x = self.xc + \
            eta*(self.R * cos(phi) * self.t1 + \
            self.R * sin(phi) * self.t2)
        return x

class sectors_t():
    def __init__(self,**kwargs):
        self.Nsector = 4
        self.sectors = [None] * self.Nsector
        xc = array(kwargs.get("xc",[0,0,0]))
        normal = array(kwargs.get("normal",[0,0,1]))
        ren = kwargs.get("ren",None)
        R = kwargs.get("R",0.04)
        S = orth_basis(normal)
        for k in range(0,self.Nsector):
            c = k%2
            self.sectors[k] = cross_mark_sector_t(k,[c,c,c])
            ren.AddActor(self.sectors[k].actor)
            cms = self.sectors[k]
            cms.xc = xc
            cms.t1 = S[1]
            cms.t2 = S[2]
            cms.R = R
        self.xc = zeros((3))

    def reset_coords(self):
        for k in range(0,self.Nsector):
            self.sectors[k].xc = self.xc
            ## copyto(self.sectors[k].xc,self.xc)
            self.sectors[k].reset_coords()

class chord_t(tube):
    def __init__(self,xe,**kwargs):
        self.xe = copy(xe)
        self.R = kwargs.get("R",0.05)
        tube.__init__(self,20,2,**kwargs)
        
    def params(self):
        xends = copy(self.xe)
        return (xends,self.R)

class tube_array_t():
    ## It's a topology of bars, defined by
    ## xnod,icone. Only xnod must change
    ## at each time step
    def __init__(self,xnod,icone,ren,**kwargs):
        self.xnod = copy(xnod)
        self.icone = copy(icone)
        self.Nnod = xnod.shape[0]
        self.Ntubes = icone.shape[0]
        self.tubes = [None] * self.Ntubes
        ndim = 3
        for k in range(0,self.Ntubes):
            xe = zeros((2,ndim))
            xe[0] = self.xnod[icone[k,0]]
            xe[1] = self.xnod[icone[k,1]]
            self.tubes[k] = chord_t(copy(xe),**kwargs)
            ren.AddActor(self.tubes[k].actor)

    def reset_coords(self):
        for k in range(0,self.Ntubes):
            xe = self.tubes[k].xe
            xe[0] = self.xnod[self.icone[k,0]]
            xe[1] = self.xnod[self.icone[k,1]]
            self.tubes[k].reset_coords()
            
class box_t(tube_array_t):
    def __init__(self,ren,**kwargs):
        x0 = array(kwargs.get("x0",[0.0,0.0,0.0]))
        L = array(kwargs.get("L",[1.0,1.0,1.0]))
        self.xnod0 = array([[0.0, 0.0, 0.0],
                            [1.0, 0.0, 0.0],
                            [1.0, 1.0, 0.0],
                            [0.0, 1.0, 0.0],
                            [0.0, 0.0, 1.0],
                            [1.0, 0.0, 1.0],
                            [1.0, 1.0, 1.0],
                            [0.0, 1.0, 1.0]])
        xn0 = self.xnod0
        for k in range(0,3):
            xn0[:,k] = x0[k] + L[k]*xn0[:,k]
        self.xnod = copy(self.xnod0)
        self.icone = array([[0, 1],
                 [1, 2],
                 [2, 3],
                 [3, 0],
                 [0, 4],
                 [1, 5],
                 [2, 6],
                 [3, 7],
                 [4, 5],
                 [5, 6],
                 [6, 7],
                 [7, 4]])
        tube_array_t.__init__(self,self.xnod0,self.icone,ren,**kwargs)
        self.dx = 0.0
        
    def reset_coords(self):
        dx = array([self.dx,0,0])
        for k in range(0,8):
            self.xnod[k] = self.xnod0[k] + dx
        tube_array_t.reset_coords(self)

## Makes an arrow. Usage:
## arrow = arrow_t(color=,tipl=,tipr=,shaftr=)
## tipl= tip length to arrow length ratio
## tipr= tip radius to arrow length ratio
## shaftr= shaft radius to arrow length ratio
## At each time step:
## arrow.x = arrow start position
## arrow.v = arrow vector
class arrow_t():
    def __init__(self,x=None,v=None,scale=1.0,**kwargs):
        self.color = [1.0,1.0,1.0]
        if 'color' in kwargs:
            self.color = kwargs['color']

        ## Geometry shape when vector has length reflen
        self.tipl = kwargs.get("tipl",0.1)
        self.tipr = kwargs.get("tipr",0.04)
        self.shaftr = kwargs.get("shaftr",0.02)

        arw = vtk.vtkArrowSource()
        self.arw = arw
        arw.SetTipResolution(20)
        arw.SetShaftResolution(20)
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(arw.GetOutputPort())
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        self.actor = actor
        self.mapper = mapper
        self.x = x
        self.v = v
        self.reflen = 0.2
        self.reset_coords()

    def reset_coords(self):
        x = self.x
        v = self.v
        if x is None:
            x = array([0,0,0])
        if v is None:
            v = array([1,0,0])
        a = self.actor
        a.SetPosition(x)
        vv = copy(v)
        ## print vv 
        ## avoid singularity when x,y components are 0
        if vv[0]==0.0:
            vv[0] = 1e-5
        rho = sqrt(vv[0]**2+vv[1]**2)
        theta = math.atan2(vv[2],rho)*180.0/math.pi
        phi = math.atan2(vv[1],vv[0])*180.0/math.pi
        a.SetOrientation(0.0,0.0,0.0)
        # a.RotateY(-theta)
        # a.RotateZ(phi)
        a.RotateWXYZ(-theta,0,1,0)
        a.RotateWXYZ(phi,0,0,1)
        a.SetScale(norm(v))
        a.GetProperty().SetDiffuseColor(self.color)
        av = norm(vv)
        s = 1.0/max(0.01,av/self.reflen)
        self.arw.SetTipLength(self.tipl*s)
        self.arw.SetTipRadius(self.tipr*s)
        self.arw.SetShaftRadius(self.shaftr*s)

class torus_t(quad_patch):
    ## Nphi along big cirumference
    ## Ntheta along small cirumference
    def __init__(self,Nphi=40,Ntheta=10,**kwargs):
        ndim = 3
        self.xc = kwargs.get("xc",zeros((ndim)))
        self.R1 = kwargs.get("R1",1.0)
        self.R2 = kwargs.get("R2",0.1)
        self.normal = kwargs.get("normal",array([0.0,0.0,1.0]))
        quad_patch.__init__(self,Nphi,Ntheta,**kwargs)
        
    def pre_reset_coords(self,**kwargs):
        ## (self.xc,self.R1,self.R2,self.normal) = self.params()
        S = orth_basis(self.normal)
        self.t1 = S[1]
        self.t2 = S[2]

    def xfun(self,xi,eta):
        phi = (xi-0.5)*2.0*pi
        theta = (eta-0.5)*2.0*pi
        rho = self.R1 + self.R2*cos(theta)
        Z = self.R2*sin(theta)
        X = rho*cos(phi)
        Y = rho*sin(phi)
        x = self.xc + Z*self.normal + X*self.t1 + Y*self.t2
        return x

class curved_tube_t():
    def __init__(self,ren=None,npoints=80,**kwargs):
        self.npoints = 80 ## Nbr of points in tail
        tubgrid = vtk.vtkPolyData()
        nodes = vtk.vtkPoints()                # positions of particles
        color = kwargs.get("color",[1,1,1])
        rtube = kwargs.get("rtube",0.03)
        ndim = 3
        self.x = zeros((npoints,ndim))
        ids = vtk.vtkIdList()                   # ids of particles
        ## Init to [0,1] on x-axis
        for k in xrange(0,npoints):
            xi = float(k)/npoints
            self.x[k] = array([xi,0,0])
            nodes.InsertNextPoint(self.x[k])
            ids.InsertNextId(k)
        tubgrid.SetPoints(nodes)

        elem = vtk.vtkPolyLine()
        elem.GetPointIds().SetNumberOfIds(npoints)
        tids = elem.GetPointIds()

        tubgrid.Allocate(1,0)
        for k in xrange(0,npoints):
            tids.SetId(k,k)
        tubgrid.InsertNextCell(elem.GetCellType(),tids)
        tuber = vtk.vtkTubeFilter()
        tuber.SetRadius(0.4*rtube)
        tuber.SetNumberOfSides(7)

        tuber.SetInputData(tubgrid)
        tbmapper = vtk.vtkDataSetMapper()
        tbmapper.SetInputConnection(tuber.GetOutputPort())
        tbactor = vtk.vtkActor()
        tbactor.SetMapper(tbmapper)
        prop = tbactor.GetProperty()
        prop.SetDiffuseColor(color)
        self.actor = tbactor
        self.nodes = nodes
        if ren:
            ren.AddActor(tbactor)

    def pre_reset_coords(self,**kwargs):
        None

    def reset_coords(self):
        self.pre_reset_coords()
        for k in range(0,self.npoints):
            self.nodes.SetPoint(k,self.x[k])
        self.nodes.Modified()
    
## makes a text 
class text_t():
    def __init__(self,**kwargs):
        text = vtk.vtkVectorText()
        text.SetText("Hello")
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(text.GetOutputPort())
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        actor.SetScale(0.02)
        actor.RotateY(180.0)
        self.actor = actor
        self.mapper = mappe

def vtk2array(vtk_array):
    at = vtk_array.GetDataType()
    if at == 11:
        #vtkDoubleArray
        pt='d'
    elif at == 12:
        #vtkIdTypeArray
        pt='l'
    #this is slow. numpy.zeros would be faster.
    # r = array(pt, [0]*vtk_array.GetSize())
    r = zeros((vtk_array.GetSize()),dtype=double)
    vtk_array.ExportToVoidPointer(r)
    return r

## Makes each sector of the rocket. Usually we use 4 sectors
## YBYB (Yellow-Black...) so as to visualize better the rotation
class rocket_sector_t(quad_patch):
    def __init__(self,rocket,indx,color,Nphi=5,Nlen=10,
                 Nsector=4,**kwargs):
        ## Each sector contains only its index and a
        ## ptr to the rocket (the father)
        self.indx = indx
        self.rocket = rocket
        quad_patch.__init__(self,Nphi,Nlen,color=color,**kwargs)

    def pre_reset_coords(self,**kwargs):
        None

    def xfun(self,xi,eta):
        ## The (father) rocket
        r = self.rocket
        ## The width of the angular sector in radians.
        ## Each sector spans azimutally dphi*[indx,indx+1]
        dphi = 2*pi/float(r.Nsector)
        ## Azimutal coordinate
        phi = (self.indx+xi)*dphi
        ## Versors
        vers = r.versors
        ## The body of the rocket spans from eta=[0,1-head]
        ## and the head is eta=[1-head,1]
        head = 0.3                        # Ratio of head length
        ## the radius at this `eta' is fac*R
        if eta<1-head:
            ## This is the cylindrical part
            fac = 1
        else:
            ## This is the semispehrical part
            fac = sqrt((1-eta)/head)
        ## Compute the coordinate for (xi,eta) in [0,1]*x[0,1]
        ## planar coords
        x = r.xc + \
            fac*r.R*(cos(phi)*vers[1]+sin(phi)*vers[2])+r.L*(eta-0.5)*vers[0]
        return x

## Makes a rocket colored by strips so as to
## show the rotation
class rocket_t():
    def __init__(self,**kwargs):
        ## Number of sectors, usually 4
        self.Nsector = 4
        ## Here we store each of the sectors who are
        ## quad_patch_t objects
        self.sectors = [None] * self.Nsector
        ## This is the center of the basis odf the rocket
        self.xc = array(kwargs.get("xc",[0,0,0]))
        ## This is the radius of the cylindrical section
        self.R = array(kwargs.get("R",0.2))
        ## Length of the rocket
        self.L = array(kwargs.get("L",1))
        ## Versors that define the orientation of the rocket
        self.versors = [array(kwargs.get("ex",[1,0,0])), \
                        array(kwargs.get("ey",[0,1,0])), \
                        array(kwargs.get("ez",[0,0,1]))]
        ## The VTK render
        self.ren = kwargs.get("ren",None)
        ## Colors of strips, normally Black/Yellow
        colors = [[0,0,0],[1,1,0]]
        colors = array(kwargs.get("colors",colors))
        ## Create each sector and set the color alternatively
        for k in range(0,self.Nsector):
            self.sectors[k] = rocket_sector_t(self,k,colors[k%2],**kwargs)
            if self.ren:
                self.ren.AddActor(self.sectors[k].actor)
        ## Set the coords of the base
        self.xc = zeros((3))

    def reset_coords(self):
        for k in range(0,self.Nsector):
            ## copyto(self.sectors[k].xc,self.xc)
            self.sectors[k].reset_coords()

## Define versors along body axes in order to see the orientation of the
## body. Another option is to use a `rocket_t' object
class versors_t():
    def __init__(self,**kwargs):
        ## The VTK render
        self.ren = kwargs.get("ren",None)
        ## The origin of the versors
        self.xc = kwargs.get("xc",array([0,0,0]))
        ## The vectors defining the versors
        w = kwargs.get("vers",None)
        if not w:
            w = zeros((3,3))
            for l in range(0,3):
                w[l,l] = 1.0
        self.vers = w
        ## These are the actual objects, which are three
        ## arrow_t objects
        self.versors = list()
        ## Scale width and head
        self.escal1 = kwargs.get("escal1",1)
        self.escal2 = kwargs.get("escal2",3)
        s1 = self.escal1
        self.ndim = 3
        for k in range(0,self.ndim):
            ## Scale each of the versors
            ## Set the colors of the versors (RGB)
            c = [0,0,0]
            c[k] = 1
            ## Add the vector to an array of vectors
            self.versors.append(arrow_t(color=c,shaftr=0.1*s1,
                                tipl=0.5*s1,tipr=0.2*s1))
            self.versors[k].x = self.xc
            if self.ren:
                self.ren.AddActor(self.versors[k].actor)

    def reset_coords(self):
        for k in range(0,self.ndim):
            self.versors[k].v = self.escal2*self.vers[k]
            self.versors[k].x = self.xc
            self.versors[k].reset_coords()
