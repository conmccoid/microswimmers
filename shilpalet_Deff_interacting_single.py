import numpy as np
# import matplotlib as mpl
# from mpl_toolkits.mplot3d import Axes3D
# import matplotlib.pyplot as plt
# from test_wrapper import stokesforcedoubletmodule
#  from test_wrapper import brinkmanforcedoubletmodule
# from test_wrapper import stericforcesmodule
from forcedoubletnetwork_wrapper import stericforcesmodule
from forcedoubletnetwork_wrapper import stokesletsvelmodule
# from forcedoubletnetwork_wrapper import networkforcesmodule
from forcedoubletnetwork_wrapper import velbetafromforcesmodule
from forcedoubletnetwork_wrapper import forcedoubletvelmodule
from forcedoubletnetwork_wrapper import forcedoubletbetamodule
import time
import math
# from mpl_toolkits.mplot3d import axes3d
import os
# import vtk
from pylab import *
import vtk
from vtk import *
import gc


# from vtk.numpy_interface import dataset_adapter as dsa
# from vtk.numpy_interface import algorithms as algs
# import vtk 


ppfac=-1.0

if(ppfac==-1.0):
	mover="pull"
else:
	mover="push"

# n = 2
# nv=n*n
trial="1"
spacing="2"
pts_file="initPos_shilpalet_cube_spacing_"+spacing+"_trial_"+trial+".txt"
initposfile=open(pts_file,"r")
lines=initposfile.readlines()
result=[]
for x in lines:
    result.append(x.split(' ')[0])
initposfile.close()
nv=float(result[0])
nv=int(nv)
print('Number of swimmers:', nv)
x0=np.zeros(3*nv)
beta=np.zeros(3*nv)
u=np.zeros(nv)
for j in range(0,nv):
	x0[j]=float(result[j+1])
	x0[j+nv]=float(result[j+nv+1])
	x0[j+2*nv]=float(result[j+2*nv+1])
	beta[j]=float(result[j+3*nv+1])
	beta[j+nv]=float(result[j+4*nv+1])
	beta[j+2*nv]=float(result[j+5*nv+1])
x=x0.copy()
pi = math.pi
#Regularization parameter for delta function
dh = 0.05#0.15/4.71
# dh_Network= 0.025
#Push=-1, Pull=1 I think

# DIPOLE MAGNITUDE
b0=(-53./3.)*dh*np.pi
#set to have U0=1.0

# STRESSLET MAGNITUDE
q0=1.0*ppfac 
# q0=+/-1 results in stagnation point at s0=22.5*dh^2*b0= -/+ 0.15609788497 if U0=1.0
# ppfac=+1 is pusher, ppfac=-1 is puller



#Brinkman Parameter
a =  0.00000000
#Self driven swimming speed
U_0=1.0
# U_0=C*p0*(45/(840*dh*np.pi))

#Poisson value
lam=2*dh/U_0


dt = 1.0e-4 #TIME STEP
Tf = 10#10.0*dt #FINAL TIME
Nt = int(Tf/dt)
Nskip = 100#int(.01/dt)
t = linspace(0.,Tf,Nt+1)
kappa =2000.0
tolfac=2.0
tol =tolfac*dh
print(tol)
print('Number of timesteps '+str(Nt))
#Number of organisms
C=.1*ppfac
Nt_store=int(Nt/Nskip)+1
x_store=np.zeros((3*nv,Nt_store))
t_store=np.zeros(Nt_store)
p=0

# beta1 = np.ones(nv)
# beta2 = np.zeros(nv)
# beta3 = np.zeros(nv)
# beta = np.zeros(3*nv)
# beta[0:nv]=beta1
# beta[nv:2*nv]=beta2
# beta[2*nv:3*nv]=beta3

tumble=np.zeros(nv)
for i in range(0,nv):
		tmbl_smpl=np.random.uniform(0,1)
		tum=-lam*np.log(1-tmbl_smpl)
		tumble[i]=tum


################################################################
######   Setting up VTK visualization for Swimmers   #########
################################################################
forcedoubletPoints   = vtk.vtkPoints()        # vtk array of points
forcedoubletVertices = vtk.vtkCellArray()     # vtk array for vertices (I'm not really sure what the diffference between points and vertices or if we even need this)
forcedoubletBeta    = vtk.vtkDoubleArray()   # just for example, I'm adding "forces" to each of the points so we need an array to hold it
forcedoubletPointId = vtk.vtkDoubleArray()   # demonstrate saving scalar data at each point, I'm just saving a id
forcedoubletSpeed = vtk.vtkDoubleArray()   # demonstrate saving scalar data at each point, I'm just saving a id
forcedoubletNonDimSpeed = vtk.vtkDoubleArray()   # demonstrate saving scalar data at each point, I'm just saving a id
# forcedoubletVelocities = vtk.vtkDoubleArray()     # vtk array for the edges


# We need to tell vtk that force has 3 components, it defaults to 1 component so if we're saving a scalar as in point_id,
# we don't need to do this step.
forcedoubletBeta.SetNumberOfComponents(3)
# forcedoubletVelocities.SetNumberOfComponents(3)
# next we can name the data we're saving
# forcedoubletVelocities.SetName("Velocity")   # name the data in the force array Force
forcedoubletBeta.SetName("Beta")   # name the data in the force array Force
forcedoubletPointId.SetName("ID")   # name the data in the point_id array ID
forcedoubletSpeed.SetName("Speed")
forcedoubletNonDimSpeed.SetName("NonDimSpeed")

# fill in the arrays. The points will be at (i,j,k). All the forces will be set to ((i+1)/nx, (j+1)/ny, (k+1)/nz) and the
# id of the point will be (i+j*nx+k*nx*ny)
for k in range(0,nv):
    pid = forcedoubletPoints.InsertNextPoint(x[k],x[k+nv],x[k+2*nv]) # set a point to (i,j,k) and get whatever ID vtk assigns to it
    forcedoubletVertices.InsertNextCell(1)                      # create a new cell for a vertex (I don't really understand why ...)
    forcedoubletVertices.InsertCellPoint(pid)                   # set the vertex to be at the point we just made ...
    beta_data = [beta[k],beta[k+nv],beta[k+2*nv]]  # create the force vector
    # beta_data = [beta(k),beta(k+nv),beta(k+2*nv)]  # create the force vector
    forcedoubletBeta.InsertNextTuple(beta_data)               # save the current force in the array
    forcedoubletPointId.InsertNextValue(k)     # save an id for the current point 
    forcedoubletSpeed.InsertNextValue(U_0) 
    forcedoubletNonDimSpeed.InsertNextValue(1.) 


# create a polydata object and add the points and data to it
polydata = vtk.vtkPolyData()                  # create polydata object
polydata.SetPoints(forcedoubletPoints)                    # associate the points with the polydata object
polydata.SetVerts(forcedoubletVertices)                   # save the vertices
polydata.GetPointData().AddArray(forcedoubletBeta)       # add array of forces 
# polydata.GetPointData().AddArray(forcedoubletVelocities)       # add array of velocities
polydata.GetPointData().SetScalars(forcedoubletPointId)  # add array of ids
polydata.GetPointData().SetScalars(forcedoubletSpeed)  # add array of ids
polydata.GetPointData().SetScalars(forcedoubletNonDimSpeed)  # add array of ids

# now write the file
polydata.Modified()            # no idea what this does
if vtk.VTK_MAJOR_VERSION <= 5:
    polydata.Update()

writer = vtk.vtkXMLPolyDataWriter();
writer.SetFileName('./T'+trial+'_'+mover+'/points_0.vtp');
if vtk.VTK_MAJOR_VERSION <= 5:
    writer.SetInput(polydata)
else:
    writer.SetInputData(polydata)
writer.Write()


#fig = plt.figure()
#plt.plot(xv, yv, 'o')

# fig = plt.figure()
# plt.quiver(xv, zv, beta1, beta3, units='width')
# # plt.quiver(xf, yf, f1, f2, color='blue')
# plt.show()
paramNum=int(7)
my_path = os.path.abspath(__file__)
param_net=np.zeros(paramNum)
param_net[0]=dh
param_net[1]=C
param_net[2]=q0
param_net[3]=b0
param_net[4]=kappa
param_net[5]=tol
param_net[6]=a



# yv=yv.astype(float)
fig = plt.figure()
	    
# start = time.time()
# fs=np.zeros(3*nv)
# vs = np.zeros(3*nv)
# vq = np.zeros(3*nv)
# ubetas = np.zeros(3*nv)
# ubetaq = np.zeros(3*nv)

# x_output=np.zeros((3*nv,Nt))
# outctr=0
# x_output[:,outctr]=x
# outctr=outctr+1
Nstore=int(Nt/Nskip)
# xstore=np.zeros((nv*3,Nstore))
tstore=np.zeros(Nstore)
# Ncounter=0

Deff=np.zeros(Nstore)
MSD=np.zeros(Nstore)
MSDx=np.zeros(Nstore)
MSDy=np.zeros(Nstore)
MSDz=np.zeros(Nstore)
Bx=np.zeros(Nstore)
By=np.zeros(Nstore)
Bz=np.zeros(Nstore)

print(param_net)
# if (nv<2):
# 	xstore=x
# else:
# 	xstore[:,Ncounter]=x
# Ncounter=Ncounter+1
for j in range(0,Nt):	
	# print(t[j])
	vs= np.zeros(3*nv)
	vq= np.zeros(3*nv)
	ubetas= np.zeros(3*nv)
	ubetaq= np.zeros(3*nv)
	fs= np.zeros(3*nv)
	
	# Computing steric forces
	# print('Calculating steric forces')
	stericforcesmodule.getstericforces( nv, x, fs, param_net)
	#Computing forcedoublet velocities and beta rates due to the forcedoublets
	# print('Calculating velocity from forcedoublets')
	forcedoubletvelmodule.getforcedoubletvel(x,x,beta,vq,param_net,nv,nv)
	# print('Calculating beta from forcedoublets')
	forcedoubletbetamodule.getforcedoubletbeta(x,ubetaq,beta,param_net,nv)
	#Computing forcedoublet velocities and beta rates due to steric forces
	# print('Calculating beta from forces')
	velbetafromforcesmodule.getvelbetafromforces(nv,nv,dh,beta,x,x,fs,vs,ubetas)

	for i in range(0,nv):


		if(j*dt>tumble[i]):
			tmbl_smpl=np.random.uniform(0,1)
			tum=-lam*np.log(1-tmbl_smpl)
			tumble[i]=j*dt+tum
			theta=2*np.pi*np.random.uniform(0,1)
			cosphi=np.random.uniform(-1,1)
			phi=np.arccos(cosphi)
			beta[i]=np.sin(phi)*np.cos(theta)
			beta[i+nv]=np.sin(phi)*np.sin(theta)
			beta[i+2*nv]=np.cos(phi)
		else:
			beta[i] = beta[i] + dt*(ubetaq[i]+ubetas[i])
			beta[i+nv] = beta[i+nv] + dt*(ubetaq[i+nv]+ubetas[i+nv])
			beta[i+2*nv] = beta[i+2*nv] + dt*(ubetaq[i+2*nv]+ubetas[i+2*nv])
			betan=0.0	
			betan=pow(np.sqrt( beta[i]*beta[i] + beta[i+nv]*beta[i+nv] + beta[i+2*nv]*beta[i+2*nv]),-1)
			beta[i]=beta[i]*betan
			beta[i+nv]=beta[i+nv]*betan
			beta[i+2*nv]=beta[i+2*nv]*betan

		u[i]=sqrt((vs[i]+vq[i])**2+(vs[i+nv]+vq[i+nv])**2+(vs[i+2*nv]+vq[i+2*nv])**2)
		x[i] = x[i] + dt*(vs[i]+vq[i])
		x[i+nv] = x[i+nv] + dt*(vs[i+nv]+vq[i+nv])
		x[i+2*nv] = x[i+2*nv] + dt*(vs[i+2*nv]+vq[i+2*nv])
			
	
	# Computing network forces
	# networkforcesmodule.getnetworkforces( numNetworkPoints, numSpringConnect, dt, springData, springConnect, networkPos, fn)
	# # Computing steric forces
	# stericforcesmodule.getstericforces( nv, x, fs, param_net)
	
	# #Computing network velocity due to force doublets 
	# forcedoubletvelmodule.getforcedoubletvel( networkPos, x, beta, vns, param_net, numNetworkPoints, nv)	
	# # Computing network velocity due to network forces
	# stokesletsvelmodule.getstokesletsvel( networkPos, networkPos, vnf, fn, dh_Network, numNetworkPoints)

	#Computing forcedoublet velocities and beta rates due to the forcedoublets
	# forcedoubletvelmodule.getforcedoubletvel(x,x,beta,vs,param_net,nv,nv)
	# forcedoubletbetamodule.getforcedoubletbeta(x,ubetas,beta,param_net,nv)
	#Computing forcedoublet velocities and beta rates due to network forces
	# velbetafromforcesmodule.getvelbetafromforces(nv,numNetworkPoints,dh_Network,beta,x,networkPos,fn,vf,ubetaf)
	# velbetafromforcesmodule.getvelbetafromforces(nv,numNetworkPoints,dh_Network,beta,x,networkPos,fn,vsf,ubetasf)
	
	
	

	if ((j+1)%(Nskip)==0):
		print(t[j+1])
		print('Starting to save')
		for l in range(0,nv):
			dx=(x[l]-x0[l])
			dy=(x[l+nv]-x0[l+nv])
			dz=(x[l+2*nv]-x0[l+2*nv])
			MSD[p]=MSD[p]+(dx*dx+dy*dy+dz*dz)
			MSDx[p]=MSDx[p]+(dx*dx)
			MSDy[p]=MSDy[p]+(dy*dy)
			MSDz[p]=MSDz[p]+(dz*dz)
			Bx[p]=Bx[p]+np.abs(beta[l])
			By[p]=By[p]+np.abs(beta[l+nv])
			Bz[p]=Bz[p]+np.abs(beta[l+2*nv])
		MSD[p]=MSD[p]/nv
		MSDx[p]=MSDx[p]/nv
		MSDy[p]=MSDy[p]/nv
		MSDz[p]=MSDz[p]/nv
		Bx[p]=Bx[p]/nv
		By[p]=By[p]/nv
		Bz[p]=Bz[p]/nv


		Deff[p]=MSD[p]/(6*t[j+1])
		tstore[p]=t[j+1]
		# xstore[:,p]=x
		p=p+1
		################################################################
		######   Setting up VTK visualization for Swimmers   #########
		################################################################
		print('Starting to allocate')
		forcedoubletPoints   = vtk.vtkPoints()        # vtk array of points
		forcedoubletVertices = vtk.vtkCellArray()     # vtk array for vertices (I'm not really sure what the diffference between points and vertices or if we even need this)
		forcedoubletBeta    = vtk.vtkDoubleArray()   # just for example, I'm adding "forces" to each of the points so we need an array to hold it
		forcedoubletPointId = vtk.vtkDoubleArray()   # demonstrate saving scalar data at each point, I'm just saving a id
		forcedoubletSpeed = vtk.vtkDoubleArray()   # saving speed
		forcedoubletNonDimSpeed = vtk.vtkDoubleArray()   # saving speed
		# forcedoubletVelocities = vtk.vtkDoubleArray()     # vtk array for the edges


		# We need to tell vtk that force has 3 components, it defaults to 1 component so if we're saving a scalar as in point_id,
		# we don't need to do this step.
		forcedoubletBeta.SetNumberOfComponents(3)
		# forcedoubletVelocities.SetNumberOfComponents(3)
		# next we can name the data we're saving
		# forcedoubletVelocities.SetName("Velocity")   # name the data in the force array Force
		forcedoubletBeta.SetName("Beta")   # name the data in the force array Force
		forcedoubletPointId.SetName("ID")   # name the data in the point_id array ID
		forcedoubletSpeed.SetName("Speed")   # name the data in the point_id array ID
		forcedoubletNonDimSpeed.SetName("NonDimSpeed")   # name the data in the point_id array ID

		# fill in the arrays. The points will be at (i,j,k). All the forces will be set to ((i+1)/nx, (j+1)/ny, (k+1)/nz) and the
		# id of the point will be (i+j*nx+k*nx*ny)
		print('Setting each point')
		for g in range(0,nv):
		    pid = forcedoubletPoints.InsertNextPoint(x[g],x[g+nv],x[g+2*nv]) # set a point to (i,j,k) and get whatever ID vtk assigns to it
		    forcedoubletVertices.InsertNextCell(1)                      # create a new cell for a vertex (I don't really understand why ...)
		    forcedoubletVertices.InsertCellPoint(pid)                   # set the vertex to be at the point we just made ...
		    beta_data = [beta[g], beta[g+nv], beta[g+2*nv]]  # create the force vector
		    # beta_data = [beta(k),beta(k+nv),beta(k+2*nv)]  # create the force vector
		    forcedoubletBeta.InsertNextTuple(beta_data)               # save the current force in the array
		    forcedoubletPointId.InsertNextValue(g)     # save an id for the current point 
		    forcedoubletSpeed.InsertNextValue(u[g])     # save an id for the current point 
		    forcedoubletNonDimSpeed.InsertNextValue(u[g]/U_0)     # save an id for the current point 



		# create a polydata object and add the points and data to it
		print('Starting to polydata')
		polydata = vtk.vtkPolyData()                  # create polydata object
		polydata.SetPoints(forcedoubletPoints)                    # associate the points with the polydata object
		polydata.SetVerts(forcedoubletVertices)                   # save the vertices
		polydata.GetPointData().AddArray(forcedoubletBeta)       # add array of forces 
		# polydata.GetPointData().AddArray(forcedoubletVelocities)       # add array of velocities
		polydata.GetPointData().SetScalars(forcedoubletPointId)  # add array of ids
		polydata.GetPointData().SetScalars(forcedoubletSpeed)  # add array of ids
		polydata.GetPointData().SetScalars(forcedoubletNonDimSpeed)  # add array of ids

		print('Starting to write')
		# now write the file
		polydata.Modified()            # no idea what this does
		if vtk.VTK_MAJOR_VERSION <= 5:
		    polydata.Update()
		print('Starting to finish')
		writer = vtk.vtkXMLPolyDataWriter();
		name_points= './T'+trial+'_'+mover+'/points_'+str(j+1) +'.vtp'
		writer.SetFileName(name_points);
		if vtk.VTK_MAJOR_VERSION <= 5:
		    writer.SetInput(polydata)
		else:
		    writer.SetInputData(polydata)
		print('Starting to finish')    
		writer.Write()
		# # plt.plot(x[0:nv], x[nv:2*nv],'bo',markersize=2.0)
		# plt.hold(True)
		# plt.plot(x[0:nv], x[nv:2*nv],'ro',markersize=1.0)
		# mpl.rcParams['font.sans-serif'] = "Arial"
		# # Then, "ALWAYS use sans-serif fonts"
		# mpl.rcParams['font.family'] = "sans-serif"
		# plt.plot(x[0:nv], x[nv:2*nv],'bo',markersize=0.50)
		# plt.plot(x[0:nv], x[nv:2*nv],'ro',markersize=1.0)
		# axes = plt.gca()
		# for tick in axes.get_xticklabels():
		# 	tick.set_fontname("Arial")
		# for tick in axes.get_yticklabels():
		# 	tick.set_fontname("Arial")
		# plt.savefig("/Users/ahoover/Dropbox/CollectiveActiveSwimmers/mycode/graph" + str(j+1) + ".png")
		# # xstore[:,Ncounter]=x
		# Ncounter=Ncounter+1
		
	# 	# print('Steric forces are '+str(f))	
	# 	# print('Velocity associated with forces '+str(vf))	
	# 	# print('Change in Beta associated with forces '+str(ubetaf))	
	# 	# print('Beta is '+str(beta))	
	# 	# print('Position is '+str(x))	
	# 	# print('Time is '+str(t[j]))
	# 	# if(np.sum(vf)>0):
		# 	break
		# 	# print(vs)
	del vs,vq,fs,ubetaq,ubetas
	gc.collect()
np.savetxt('./T'+trial+'_'+mover+'/Deff.curve',  np.c_[tstore,Deff], fmt='%1.4e')   # use exponential notation
np.savetxt('./T'+trial+'_'+mover+'/Deff_nondim.curve',  np.c_[tstore,Deff/(U_0*U_0*lam/3.0)], fmt='%1.4e')   # use exponential notation
np.savetxt('./T'+trial+'_'+mover+'/MSD.curve', np.c_[tstore,MSD], fmt='%1.4e')   # use exponential notation
np.savetxt('./T'+trial+'_'+mover+'/MSDx.curve', np.c_[tstore,MSDx], fmt='%1.4e')   # use exponential notation
np.savetxt('./T'+trial+'_'+mover+'/MSDy.curve', np.c_[tstore,MSDy], fmt='%1.4e')   # use exponential notation
np.savetxt('./T'+trial+'_'+mover+'/MSDz.curve', np.c_[tstore,MSDz], fmt='%1.4e')   # use exponential notation
np.savetxt('./T'+trial+'_'+mover+'/Bx.curve',  np.c_[tstore,Bx], fmt='%1.4e')   # use exponential notation
np.savetxt('./T'+trial+'_'+mover+'/By.curve',  np.c_[tstore,By], fmt='%1.4e')   # use exponential notation
np.savetxt('./T'+trial+'_'+mover+'/Bz.curve',  np.c_[tstore,Bz], fmt='%1.4e')   # use exponential notation


# ################################################################
# ######    Setting up VTK visualization for Trajectories     #########
# ################################################################

# trajectoryPoints   = vtk.vtkPoints()        # vtk array of points
# trajectoryVertices = vtk.vtkCellArray()     # vtk array for vertices (I'm not really sure what the diffference between points and vertices or if we even need this)
# trajectoryMSD    = vtk.vtkDoubleArray()   # just for example, I'm adding "forces" to each of the points so we need an array to hold it
# trajectoryPointId = vtk.vtkDoubleArray()   # demonstrate saving scalar data at each point, I'm just saving a id
# trajectoryEdges    = vtk.vtkCellArray()     # vtk array for the edges
# trajectorMSD.SetNumberOfComponents(1)

# # next we can name the data we're saving
# networkForces.SetName("Force")   # name the data in the force array Force
# trajectoryPointId.SetName("ID")   # name the data in the point_id array ID


# for k in range(0,nv_store):
#     pid = networkPoints.InsertNextPoint(x_store[],networkPos[k+numNetworkPoints],networkPos[k+2*numNetworkPoints]  )           # set a point to (i,j,k) and get whatever ID vtk assigns to it
#     networkVertices.InsertNextCell(1)                      # create a new cell for a vertex (I don't really understand why ...)
#     networkVertices.InsertCellPoint(pid)                   # set the vertex to be at the point we just made ...
#     force_data = [0.0, 0.0, 0.0]  # create the force vector
#     networkForces.InsertNextTuple(force_data)               # save the current force in the array
#     networkPointId.InsertNextValue(k)     # save an id for the current point 



# for k in range(0,nV):
#     line = vtk.vtkLine()                           # an edge is a line
#     line.GetPointIds().SetId(0,springConnect[k])      # set the start point to be the point (i,j,k), which is the (i+j*nx+k*nx*ny) point 
#     line.GetPointIds().SetId(1,springConnect[k+numSpringConnect])    # set the end point to be the point (i+1,j,k), which is the (1+i+j*nx+k*nx*ny) point 

# networkEdges.InsertNextCell(line)

# # create a polydata object and add the points and data to it
# polydata = vtk.vtkPolyData()                  # create polydata object
# polydata.SetPoints(networkPoints)                    # associate the points with the polydata object
# polydata.SetVerts(networkVertices)                   # save the vertices
# polydata.SetLines(networkEdges)                      # save the edges
# polydata.GetPointData().AddArray(networkForces)       # add array of forces 
# polydata.GetPointData().SetScalars(networkPointId)  # add array of ids


# # now write the file
# polydata.Modified()            # no idea what this does
# if vtk.VTK_MAJOR_VERSION <= 5:
#     polydata.Update()

# writer = vtk.vtkXMLPolyDataWriter();
# writer.SetFileName("./visualization/network_0.vtp");
# if vtk.VTK_MAJOR_VERSION <= 5:
#     writer.SetInput(polydata)
# else:
#     writer.SetInputData(polydata)
# writer.Write()

# 	if ((j+1)%(Nskip)==0):
# 		print('vf='+str(vf[0]))
# 		print('vs='+str(vs[0]))
# 		print('dt='+str(dt))

# 		print('x='+str(x[0]))
# 		print('y='+str(x[nv]))
# 		print('z='+str(x[nv*2]))
# 		print('fx='+str(f[0]))
# 		print('fx='+str(f[nv]))
# 		print('fz='+str(f[nv*2]))
# 		print('betax='+str(beta[0]/np.sqrt(beta[0]*beta[0]+beta[nv]*beta[nv]+beta[2*nv]*beta[2*nv])))
# 		print('betay='+str(beta[nv]/np.sqrt(beta[0]*beta[0]+beta[nv]*beta[nv]+beta[2*nv]*beta[2*nv])))
# 		print('betaz='+str(beta[nv*2]/np.sqrt(beta[0]*beta[0]+beta[nv]*beta[nv]+beta[2*nv]*beta[2*nv])))
# 		# print(x[]/dh)
# 		# print(beta[0]/np.sqrt(beta[1]*beta[1]+beta[3]*beta[3]+))
# 		# print(beta[0]/np.sqrt(beta[1]*beta[1]+beta[3]*beta[3]))
# 		# plt.plot(x[0:nv], x[nv:2*nv],'bo',markersize=2.0)
# 		plt.hold(True)
# 		# plt.plot(x[0:nv], x[nv:2*nv],'ro',markersize=1.0)
# 		mpl.rcParams['font.sans-serif'] = "Arial"
# 		# Then, "ALWAYS use sans-serif fonts"
# 		mpl.rcParams['font.family'] = "sans-serif"
# 		plt.plot(x[0:nv], x[nv:2*nv],'bo',markersize=0.50)
# 		plt.plot(x[0:nv], x[nv:2*nv],'ro',markersize=1.0)
# 		axes = plt.gca()
# 		for tick in axes.get_xticklabels():
# 			tick.set_fontname("Arial")
# 		for tick in axes.get_yticklabels():
# 			tick.set_fontname("Arial")
# 		plt.savefig("/Users/ahoover/Dropbox/CollectiveActiveSwimmers/mycode/graph" + str(j+1) + ".png")
# 		xstore[:,Ncounter]=x
# 		Ncounter=Ncounter+1
		
# 		# print('Steric forces are '+str(f))	
# 		# print('Velocity associated with forces '+str(vf))	
# 		# print('Change in Beta associated with forces '+str(ubetaf))	
# 		# print('Beta is '+str(beta))	
# 		# print('Position is '+str(x))	
# 		print('Time is '+str(t[j]))
# 		# if(np.sum(vf)>0):
# 		# 	break
# 		# 	# print(vs)
		
# 	# if ((j+1)%Nt==0):

# 	# 	plt.figure(dpi=1200)
# 	# 	mpl.rcParams['font.sans-serif'] = "Arial"
# 	# 	# Then, "ALWAYS use sans-serif fonts"
# 	# 	mpl.rcParams['font.family'] = "sans-serif"		
# 	# 	plt.plot(xstore[0,0:Ncounter], xstore[2,0:Ncounter],'b-',linewidth=1.0)
# 	# 	plt.plot(xstore[1,0:Ncounter], xstore[3,0:Ncounter],'b-',linewidth=1.0)
# 	# 	plt.plot(x[0:nv], x[nv:2*nv],'ro',markersize=2.0)
# 	# 	plt.rc('text', usetex=True)
# 	# 	plt.rc('font', family='serif')
# 	# 	axes = plt.gca()
# 	# 	axes.set_xlim([-.25,5])
# 	# 	axes.set_ylim([-1.0,1.0])	
# 	# 	axes.grid(color='k', linestyle='--', linewidth=1,alpha=.5, )
# 	# 	axes.set_aspect('equal', adjustable='box')
# 	# 	plt.xlabel(r'x')
# 	# 	plt.ylabel(r'y')#,fontsize=16)
# 	# 	plt.title(r"Brinkman Fluid $\alpha = $"+str(a)+", Time = "+ str(Tf) + r", tol = "+ str(tolfac)+r" $\delta$",fontsize=16,horizontalalignment='center')#,fontname="Arial")#, color='gray')
# 	# 	for tick in axes.get_xticklabels():
# 	# 		tick.set_fontname("Arial")
# 	# 	for tick in axes.get_yticklabels():
# 	# 		tick.set_fontname("Arial")
# 	# 	plt.savefig("/Users/ahoover/Dropbox/CollectiveActiveSwimmers/mycode/StokesFluid_Time"+str(Tf)+"_tol20e-1.png")

# 	# 	# plt.savefig("/Users/hooverap/Dropbox/CollectiveActiveSwimmers/mycode/graph_lasttime.png")	
# 		# plt.quiver(x[0:nv], x[nv:2*nv], vs[0:nv],vs[nv:2*nv], units='width')
# 		# plt.savefig("/Users/ahoover/Dropbox/CollectiveActiveSwimmers/mycode/graph" + str(j) + ".png")
# 		# plt.savefig("/Users/hooverap/Dropbox/CollectiveActiveSwimmers/mycode/graph" + str(j) + ".png")
# 		# plt.show()
# 		# print(ubetas)
# 		# print(ubetaf)
# 	# if(j==Nt-1):
# 	# 	print('Velocity associated with forces '+str(vf))	
# 	# 	print('Change in Beta associated with forces '+str(ubetaf))	
# 	# 	print(vf)	
# print(x)

# plt.figure()

# print(vs)
# print(vf)
# print(ubetas)
# print(ubetaf)
# fig = plt.figure()
	

# print(my_path)
# fig = plt.figure()
# plt.quiver(xv, yv, u, v, units='width')
# plt.quiver(xf, yf, f1, f2, color='blue')
# plt.show()
