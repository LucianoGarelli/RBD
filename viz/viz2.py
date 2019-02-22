#! /usr/bin/env python

## Avoid warnings from h5py
import warnings
with warnings.catch_warnings():
    warnings.filterwarnings("ignore",category=FutureWarning)
    import h5py

## Import modules
import vtk, time, os, sys, signal, glob, re, json
import random as rndm
from numpy import *
from math import pi,cos,sin,tan,atan
from vtkutils import *
from utils import *

## Rotation around X axis
def X(alpha):
  R = array([[1,0,0], \
             [0,cos(alpha),sin(alpha)], \
             [0,-sin(alpha),cos(alpha)]])
  return R

## Rotation around Z axis
def Z(alpha):
    R = array([[cos(alpha),sin(alpha),0], \
            [-sin(alpha),cos(alpha),0], \
            [0,0,1]])
    return R

## Not used anymore, based on an expression
## product of matrices
def make_versors2(eulang,omegahbdy):
    nt = xc.shape[0]
    Oh = zeros((nt,9))
    omegah = zeros((nt,3))
    for k in range(0,nt):
        O = X(eulang[k,0])
        O = O.dot(Z(eulang[k,1]))
        O = O.dot(X(eulang[k,2]))
        Oh[k,0:3] = O[:,0]
        Oh[k,3:6] = O[:,1]
        Oh[k,6:9] = O[:,2]
        omegah[k,:] = O.dot(omegahbdy[k,:])
    return (Oh,omegah)

def make_versors(eulang,omegahbdy):
    ## Number of time steps
    nt = xc.shape[0]
    ## Each row of Oh will store the versors
    Oh = zeros((nt,9))
    ## Angular velocity in earth axes
    omegah = zeros((nt,3))
    for k in range(0,nt):
        ## Euler angles
        phi = eulang[k,0]
        theta = eulang[k,1]
        psi = eulang[k,2]
        ## Cos and Sin of Euler angles
        ct = cos(theta)
        cp = cos(psi)
        cf = cos(phi)
        st = sin(theta)
        sp = sin(psi)
        sf = sin(phi)
        ## Compute orthogonal transformation matrix from
        ## body system to earth system
        O = array([[ct*cp,sf*st*cp-cf*sp,cf*st*cp+sf*sp], \
                   [ct*sp,sf*st*sp+cf*cp,cf*st*sp-sf*cp], \
                   [-st,  sf*ct,        cf*ct]])
        ## Transform the velocity in body axes to earth
        omegah[k,:] = O.dot(omegahbdy[k,:])
        ## In Oh we store in each row the versors. As the versors
        ## are columns in O we have to transpose the matrix
        Oh[k,:] = O.transpose().reshape(9)
    return (Oh,omegah)

if None:
    h5f = h5py.File('tempo.h5','r')
    xc = h5f['gdata/value/xc/value'][:]
    xc = transpose(xc)
    Oh = h5f['gdata/value/Oh/value'][:]
    Oh = transpose(Oh)
    omegah = h5f['gdata/value/omegah/value'][:]
    omegah = transpose(omegah)
    h5f.close()
else:
    ## Body_ang_vel', -> Velocidad angular ejes cuerpo
    ## Body_vel',  -> Velocidad translacional ejes cuerpo
    ## Euler_ang', -> Anguler de Euler
    ## Inertial_coord', -> Coord cg ejes tierra
    ## Inertial_vel'] -> Velocidad ejes tierra
    ## h5f = h5py.File('Data.hdf5','r')
    ## h5f = h5py.File('Data-planoxz.hdf5','r')
    h5f = h5py.File('../Data.hdf5','r')
    xc = h5f['/Inertial_coord'][:]
    ## Third component is elevation, so we have to reverse sign
    xc[:,2] = -xc[:,2]
    eulang = h5f['/Euler_ang'][:]
    omegahbdy = h5f['/Body_ang_vel'][:]
    [Oh,omegah] = make_versors(eulang,omegahbdy)
    h5f.close()

## Nbr of time steps
nt = xc.shape[0]
## min corner of the bounding box
x0 = amin(xc,0)
## Compute array of lengths
Lbox = amax(xc,0)
Lbox = Lbox-x0
## Maximum total length, to be a scale of the visualization
L = amax(Lbox)

## Ndims and time step
ndim = 3
Dt = 0.004

## Visualize from frame1 to frame2, if frame2<0 proceed
## until exhausting all frames in the input vectors
frame1 = 0
frame2 = -1
## Save the frames so as to make a video or not
mkvideo = 1

## Increment in frames. Accelerates the visualization
frame_inc = 5
## If true
interactive = False

## Used for a flying camera
## Initial azimutal position of camera
phi0 = -65
phi0Rcam = -90
static_cam = 0

## Speed of rotation camera
wcam = 0
wRcam = 0.0
Rcammax = 30
Rcammin = Rcammax
## Elevation of camera
camelev = 30.0
if mkvideo:
    frame_inc = 2

# Create the usual VTK rendering stuff.
ren = vtk.vtkRenderer()
renWin = vtk.vtkRenderWindow()
renWin.AddRenderer(ren)
renWin.SetSize(800,600)
## Color of background
ren.SetBackground(.05, .1, .2)

## Scales a rocket representation of the Omega vector
wscal1 = 6 ## Scale the width and head of the Omega vector
wscal2 = 3 ## Scale the length of the Omega vector
omega = arrow_t(color=[1,0,1],shaftr=0.5*wscal1,
                tipl=2.5*wscal1,tipr=1*wscal1)
ren.AddActor(omega.actor)
frame=0
## Set the initial position and vector
omega.x = xc[frame][0:3]
omega.v = 0.5*omegah[frame]

## Define versors along body axes
versors = versors_t(escal1=25,escal2=10,ren=ren)

## A box that is computed from the bounding box of the movement of the
## body
box = box_t(x0=x0,L=Lbox,color=[0,1,1],ren=ren,R=0.5)

## A rocket for the orientation of the body (you can use the versors or the
## rocket).
## ref:= Refines the mesh that defines the surface of the rocket
ref=5
rocket = rocket_t(ren=ren,R=2,L=12,Nphi=5*ref,Nlen=10*ref)

## Set the camera
cam = ren.GetActiveCamera()
## Radius of flying camera. If Rcammax!=Rcammin
## then we have a camera with varying distances to
## the view point
Rcammax = 1.8*L
Rcammin = Rcammax

## Focal point of the camera
xto = x0+0.5*Lbox
## Position of the camera
xfrom = xto+Rcammax*array([-0.5,1,-0.2])
cam.SetFocalPoint(xto)
cam.SetPosition(xfrom)
## In this case we reverse the ViewUp vector
## sue to the fact that ni the eqs. the z coordinate
## is directed downwards
cam.SetViewUp(0.,0.,-1.);

## Lights are automatically created by VTK.
## Warning: If you create lights then the auto lights are turned off.
# Create a vtkLight, and set the light parameters.
if 0:
    light1 = vtk.vtkLight()
    light1.SetFocalPoint(xto)
    xl = xto+array([-0.5*L,-L,0])
    light1.SetPosition(xl)
    ren.AddLight(light1)
    light2 = vtk.vtkLight()
    light2.SetFocalPoint(xto)
    xl = xto+array([+0.5*L,-L,0])
    light2.SetPosition(xl)
    ren.AddLight(light2)
    lightKit = vtk.vtkLightKit()
    lightKit.SetKeyLightIntensity(0.35)
    lightKit.AddLightsToRenderer(ren)

if 0:
    lightKit = vtk.vtkLightKit()
    lightKit.SetKeyLightIntensity(0.35)
    lightKit.AddLightsToRenderer(ren)

## Create axes RGB at the origin
axes = draw_axes([0.,0.,0.],2,0.05)
ren.AddActor(axes[0])
ren.AddActor(axes[1])
ren.AddActor(axes[2])

## If interactive you can manipulate the scene
if interactive:
    iren = vtk.vtkRenderWindowInteractor()
    iren.SetRenderWindow(renWin)

    style = vtk.vtkInteractorStyleTrackballCamera()
    iren.SetInteractorStyle(style)
    iren.Initialize()
    iren.Start()
    sys.exit()

# Render the scene
renWin.Render()

## This is used to store the frames
## for creating a movie
w2i = vtk.vtkWindowToImageFilter()
w2i.SetInput(renWin)
w2i.Update()

## The writer
writer = vtk.vtkPNGWriter()
writer.SetInputConnection(w2i.GetOutputPort())
# writer.SetCompressionToJPEG()

## This is added so that it gives time to set no border in the OpenGL window
## and other stuff like minimizing other windows.
if mkvideo:
    renWin.Render()
    ## raw_input("Enter something to continue: ")
    renWin.BordersOff()
    ## Create a clean dir "./frames" in order to store
    ## the frames there
    if os.path.isdir("./frames"):
        if query_yes_no("./frames dir exist, delete? ","no"):
            os.system("/bin/rm -rf ./frames")
        else:
            print "ABORT!"
            sys.exit()
    os.mkdir("./frames")

# Frame loop: rotating the camera and modify object coordinates
last = time.time()
frame = frame1
zframe = 0
while 1:
    now = time.time()
    if now>last+2:
        print "frame ",frame,time.asctime()
        last = now

    t = frame*Dt
    if frame>=nt:
        print "Exhausted frames. Exiting. nframes ",nt
        break
    ## Set the coords for the Omega arrow
    omega.x = xc[frame][0:3]
    omega.v = wscal2*omegah[frame]
    omega.reset_coords()

    ## Set the coords and orientation of the rocket
    rocket.versors = Oh[frame].reshape(3,3)
    rocket.xc = xc[frame][0:3]
    rocket.reset_coords()

    ## Set coords and orientation of the versors
    versors.xc = xc[frame][0:3]
    for k in range(0,ndim):
        k0 = k*ndim
        versors.vers[k] = Oh[frame][k0:k0+3]
    versors.reset_coords()

    ## Update camera (if fyling camera is used)
    if 0:
        phi = phi0*pi/180.0 + wcam*frame
        Rcamav =(Rcammin+Rcammax)/2.0
        DRcam = (Rcammax-Rcammin)/2.0
        phiR = phi0Rcam*pi/180.0 + wRcam*frame
        cer = camelev*pi/180.0
        Rcam = (Rcamav + DRcam*sin(phiR))
        Hcam = sin(cer)*Rcam
        Rcam *= cos(cer)
        view_dir = array([Rcam*cos(phi),Rcam*sin(phi),Hcam])
        pfrom = xto + view_dir
        cam.SetPosition(pfrom)

    renWin.Render()
    ren.ResetCameraClippingRange()

    ## Save current frame
    if mkvideo:
        w2i.Modified()
        png = "./frames/frame.%04d.png" % zframe
        writer.SetFileName(png)
        writer.Write()
        zframe += 1
    else:
        time.sleep(0.05)

    if frame2>=0 and frame >= frame2:
        break

    ## Update frame counter
    frame += frame_inc
