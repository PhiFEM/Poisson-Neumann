from __future__ import print_function
import numpy as np
from dolfin import *
import sympy
import matplotlib.pyplot as plt
parameters['allow_extrapolation'] = True
parameters["form_compiler"]["representation"] = 'uflacs'

# Polynome of the spaces
polV = 1
polPhi = polV+1
polQuadr = 2*(polV+polPhi)-2

# test case
test_case = 2

# Size of the domain in the test case 1
R = 0.47

# size of the domain in the test case 2
R1 = 1
R2 = 2

# Rotation of the mesh
theta0=0.0


# Construction of phi
class phi_expr(UserExpression):
	def eval(self, value, x):
		xxx = x[0]
		yyy = x[1]
		if test_case ==1:
			rrr = (xxx**2+yyy**2)**(0.5)
			if xxx!=0:
				theta = np.arctan2(yyy,xxx) - theta0
			if xxx==0 and yyy>0:
				theta = 0.5*np.pi - theta0
			if xxx==0 and yyy<0:
				theta = -0.5*np.pi - theta0
			if xxx==0 and yyy==0:
				theta = 0 
			value[0] = rrr**4*(5.0 + 3.0*np.sin(7.0*theta + 7*np.pi/36.0))/2.0 - R**4
		if test_case==2:
			rot_x = np.cos(theta0-np.pi/8)*xxx-np.sin(theta0-np.pi/8)*yyy
			rot_y = np.sin(theta0-np.pi/8)*xxx+np.cos(theta0-np.pi/8)*yyy
			value[0] = np.max([abs(rot_x)/R1,abs(rot_y)/R2])-1

	def value_shape(self):
		return (2,)


# Construction of the mesh
i=6
N = int(10*2**((i)))
if test_case ==1:
	mesh_macro = RectangleMesh(Point(-1.0, -1.0), Point(1.0, 1.0), N, N)
if test_case ==2:
	size = float(1.1*(R1**2+R2**2)**(0.5))
	mesh_macro = RectangleMesh(Point(-size, -size), Point(size, size), N, N)
V_phi = FunctionSpace(mesh_macro, 'CG',polPhi)
phi = phi_expr(element=V_phi.ufl_element())
phi = interpolate(phi,V_phi)
domains = MeshFunction("size_t", mesh_macro, mesh_macro.topology().dim())
domains.set_all(0)
for ind in range(mesh_macro.num_cells()):
	mycell = Cell(mesh_macro,ind)
	v1x,v1y,v2x,v2y,v3x,v3y = mycell.get_vertex_coordinates()
	if phi(v1x,v1y)<=0.0 or phi(v2x,v2y)<=0.0 or phi(v3x,v3y)<=0.0 or near(phi(v1x,v1y),0.0) or near(phi(v2x,v2y),0.0) or near(phi(v3x,v3y),0.0):
		domains[ind] = 1
mesh = SubMesh(mesh_macro, domains, 1)


# Computation of the source term
V=FunctionSpace(mesh,'CG',1)
f= Function(V)
plot_sol1 = plot(f)


# Construction of the mesh
i=1
N = int(10*2**((i)))
if test_case ==1:
	mesh_macro = RectangleMesh(Point(-1.0, -1.0), Point(1.0, 1.0), N, N)
if test_case ==2:
	size = float(1.1*(R1**2+R2**2)**(0.5))
	mesh_macro = RectangleMesh(Point(-size, -size), Point(size, size), N, N)
V_phi = FunctionSpace(mesh_macro, 'CG',polPhi)
phi = phi_expr(element=V_phi.ufl_element())
phi = interpolate(phi,V_phi)
domains = MeshFunction("size_t", mesh_macro, mesh_macro.topology().dim())
domains.set_all(0)
for ind in range(mesh_macro.num_cells()):
	mycell = Cell(mesh_macro,ind)
	v1x,v1y,v2x,v2y,v3x,v3y = mycell.get_vertex_coordinates()
	if phi(v1x,v1y)<=0.0 or phi(v2x,v2y)<=0.0 or phi(v3x,v3y)<=0.0 or near(phi(v1x,v1y),0.0) or near(phi(v2x,v2y),0.0) or near(phi(v3x,v3y),0.0):
		domains[ind] = 1
mesh = SubMesh(mesh_macro, domains, 1)


# Facets and cells where we apply the ghost penalty
mesh.init(1,2)
vertex_ghost = MeshFunction("size_t", mesh, mesh.topology().dim()-2)
facet_ghost = MeshFunction("size_t", mesh, mesh.topology().dim()-1)
cell_ghost = MeshFunction("size_t", mesh, mesh.topology().dim())
vertex_ghost.set_all(0)
facet_ghost.set_all(0)
cell_ghost.set_all(0)

for mycell in cells(mesh):
	for myfacet in facets(mycell):
		v1, v2 = vertices(myfacet)
		if phi(v1.point().x(),v1.point().y())*phi(v2.point().x(),v2.point().y())<=0 or near(phi(v1.point().x(),v1.point().y())*phi(v2.point().x(),v2.point().y()),0.0):
			cell_ghost[mycell] = 1
			for myfacet2 in facets(mycell):
				facet_ghost[myfacet2] = 1
				v1, v2 = vertices(myfacet2)
				vertex_ghost[v1] = 1
				vertex_ghost[v2] = 1

for mycell in cells(mesh):
	if cell_ghost[mycell] == 0:
		for myfacet in facets(mycell):
			v1, v2 = vertices(myfacet)
			if facet_ghost[myfacet] == 1:
				facet_ghost[myfacet] = 2
			if vertex_ghost[v1] ==0 and vertex_ghost[v2] ==0:
				facet_ghost[myfacet] = 3
	v1, v2, v3 = vertices(mycell)
	if vertex_ghost[v1]==0 and vertex_ghost[v2]==0 and vertex_ghost[v3]==0:
		cell_ghost[mycell] = 2
submesh = SubMesh(mesh, domains, 1)


# Plot and save
plot_sol2 = plot(mesh)
if test_case ==1:
	plt.savefig('test_case_1_domain.png')
if test_case ==2:
	plt.savefig('test_case_2_domain.png')
