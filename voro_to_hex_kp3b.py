# generate pure hex mesh for random pebble bed in a cylinder container
# using voronoi cell to split domain for each pebble.
# each voronoi cell is then divided into hex mesh that will conformal pebbles  
#=================================================================================
#
# one voro quad to four nek quads
# add chamfers
# 
#=================================================================================
# 
# use new method for top and bot ghost pebbles
# direct extract top and bot ghost pebbles from the DEM data
# but offset in z direction to avoid perpendicular wall.
#

# python module
import numpy as np
import math
from scipy.spatial import Voronoi, voronoi_plot_2d
import os
import time
# not working to include this, not know why
#from matplotlib.pyplot import *

#=================================================================================
# user defined files
from data_class import *
from merge_vertices import *
from splitting_subroutines import *
from polygon_divide import *
from write_rea import *
from dump_vtk import *
from fix_elements import *
from check_mesh import *
from top_quads_smooth import *
#=================================================================================

starttime = time.time()

#=================================================================================
# part 1
# read pebble coordinates, cylinder infos

scaling=1.0

#case = '18pebbles' 	# 18 pebbles case
#case = 'tamu146'	# tamu 146 pebbes case
#case = '1568'	# tamu 146 pebbes case
#case = 'kp1'
case = 'kp3'

ghost_pebble_tol_side = 1.5
ghost_pebble_tol = 2.0

#r_chamfer = 0.2 # chamfer radius relative to pebble radius
r_chamfer = 0.1 # chamfer radius relative to pebble radius
#r_chamfer = 0.15 # chamfer radius relative to pebble radius


if case == '18pebbles':
	pebble_file_path = 'pebbles.txt'
	pebble_file = open(pebble_file_path) 
	cyl_radius = float(pebble_file.readline().split()[0])*scaling
	line = pebble_file.readline()
	cyl_bot = float(line.split()[0])*scaling
	cyl_top = float(line.split()[1])*scaling

	pebble_diameter = float(pebble_file.readline().split()[0])*scaling
	pebble_radius = pebble_diameter/2.0

	number_of_pebbles = int(pebble_file.readline().split()[0])

	random_pebbles = pebbles('random_pebbles') # instance of class pebbles

	for ipebble in range(0,number_of_pebbles):
		line = pebble_file.readline()
		p_x = float(line.split()[0])
		p_y = float(line.split()[1])
		p_z = float(line.split()[2])
		ifghost = False
		tag = ''
		random_pebbles.new_pebble(p_x,p_y,p_z,ifghost,tag)    # add new pebble

elif case == 'tamu146':
	#pebble_file_path = 'tamu_pebbles.txt' # non-touching pebbles
	pebble_file_path = 'tamu_pebbles_touching.txt' # touching pebbles
	pebble_file = open(pebble_file_path) 

	cyl_radius = float(pebble_file.readline().split()[0])

	line = pebble_file.readline()
	cyl_bot = float(line.split()[0])
	cyl_top = float(line.split()[1])

	pebble_diameter = float(pebble_file.readline().split()[0])
	
	# ensure pebble diameter is one
	scaling = 1.0/pebble_diameter

	pebble_diameter = 1.0
	pebble_radius = pebble_diameter/2.0
	
	cyl_radius = cyl_radius*scaling
	cyl_bot = cyl_bot*scaling
	cyl_top = cyl_top*scaling

	number_of_pebbles = int(pebble_file.readline().split()[0])

	random_pebbles = pebbles('random_pebbles') # instance of class pebbles

	for ipebble in range(0,number_of_pebbles):
		line = pebble_file.readline()
		p_x = float(line.split(',')[0])*0.0254*scaling
		p_y = float(line.split(',')[1])*0.0254*scaling
		p_z = float(line.split(',')[2])*0.0254*scaling
		ifghost = False
		tag = ''
		random_pebbles.new_pebble(p_x,p_y,p_z,ifghost,tag)    # add new pebble
	
elif case == '1568':
	pebble_file_path = '1568_pebble.txt' # touching pebbles
	pebble_file = open(pebble_file_path) 

	cyl_radius = float(pebble_file.readline().split()[0])

	line = pebble_file.readline()
	cyl_bot = float(line.split()[0])
	cyl_top = float(line.split()[1])

	pebble_diameter = float(pebble_file.readline().split()[0])
	
	# ensure pebble diameter is one
	scaling = 1.0/pebble_diameter

	pebble_diameter = 1.0
	pebble_radius = pebble_diameter/2.0
	
	cyl_radius = cyl_radius*scaling
	cyl_bot = cyl_bot*scaling
	cyl_top = cyl_top*scaling

	number_of_pebbles = int(pebble_file.readline().split()[0])

	random_pebbles = pebbles('random_pebbles') # instance of class pebbles

	mx_rr = 0.0
	for ipebble in range(0,number_of_pebbles):
		line = pebble_file.readline()
		p_x = float(line.split()[0])*scaling/1.5 # 1.5 is the needed as inconsistance
		p_y = float(line.split()[1])*scaling/1.5
		p_z = float(line.split()[2])*scaling/1.5
		
		rr = (p_x**2.0+p_y**2.0)**0.5
		
		if rr > mx_rr: mx_rr = rr
		
		ifghost = False
		tag = ''
		random_pebbles.new_pebble(p_x,p_y,p_z,ifghost,tag)    # add new pebble
	
	print 'mx_rr is ',mx_rr
	print 'actual cyl_radius should be ',mx_rr+0.5
	print 'given cyl_radius is ',cyl_radius
	#cyl_radius = mx_rr + 0.5
	
elif case == 'kp1':
	pebble_file_path = 'tsk1_pebble_kp.dat' # touching pebbles
	pebble_file = open(pebble_file_path) 

	cyl_radius = float(pebble_file.readline().split()[0])

	line = pebble_file.readline()
	cyl_bot = float(line.split()[0])
	cyl_top = float(line.split()[1])

	pebble_diameter = float(pebble_file.readline().split()[0])
	
	# ensure pebble diameter is one
	scaling = 1.0/pebble_diameter

	pebble_diameter = 1.0
	pebble_radius = pebble_diameter/2.0
	
	cyl_radius = cyl_radius*scaling
	cyl_bot = cyl_bot*scaling
	cyl_top = cyl_top*scaling

	number_of_pebbles = int(pebble_file.readline().split()[0])

	random_pebbles = pebbles('random_pebbles') # instance of class pebbles

	mx_rr = 0.0
	for ipebble in range(0,number_of_pebbles):
		line = pebble_file.readline()
		p_x = float(line.split()[0])*scaling
		p_y = float(line.split()[1])*scaling
		p_z = float(line.split()[2])*scaling
		
		rr = (p_x**2.0+p_y**2.0)**0.5
		
		if rr > mx_rr: mx_rr = rr
		
		ifghost = False
		tag = ''
		random_pebbles.new_pebble(p_x,p_y,p_z,ifghost,tag)    # add new pebble
	
	print 'mx_rr is ',mx_rr
	print 'actual cyl_radius should be ',mx_rr+0.5
	print 'given cyl_radius is ',cyl_radius
	#cyl_radius = mx_rr + 0.5
	
elif case == 'kp3':
	pebble_file_path = 'tsk3_pebble_kp_new.dat' # touching pebbles
	pebble_file = open(pebble_file_path) 

	cyl_radius = float(pebble_file.readline().split()[0])

	line = pebble_file.readline()
	cyl_bot = float(line.split()[0])
	cyl_top = float(line.split()[1])

	pebble_diameter = float(pebble_file.readline().split()[0])
	
	# ensure pebble diameter is one
	scaling = 1.0/pebble_diameter

	pebble_diameter = 1.0
	pebble_radius = pebble_diameter/2.0
	
	cyl_radius = cyl_radius*scaling
	cyl_bot = cyl_bot*scaling
	cyl_top = cyl_top*scaling

	number_of_pebbles = int(pebble_file.readline().split()[0])

	random_pebbles = pebbles('random_pebbles') # instance of class pebbles

	mx_rr = 0.0
	for ipebble in range(0,number_of_pebbles):
		line = pebble_file.readline()
		p_x = float(line.split()[0])*scaling
		p_y = float(line.split()[1])*scaling
		p_z = float(line.split()[2])*scaling
		
		rr = (p_x**2.0+p_y**2.0)**0.5
		
		if rr > mx_rr: mx_rr = rr
		
		ifghost = False
		tag = ''
		random_pebbles.new_pebble(p_x,p_y,p_z,ifghost,tag)    # add new pebble
	
	print 'mx_rr is ',mx_rr
	print 'actual cyl_radius should be ',mx_rr+0.5
	print 'given cyl_radius is ',cyl_radius
	#cyl_radius = mx_rr + 0.5	

number_of_pebbles = int(random_pebbles.number_of_pebbles())
print '========================================='
print 'number of real pebbles: ' + str(number_of_pebbles)

# construct ghost pebbles
# 
# if one pebble is 2 pebble-diameters away from boundary, 
# then define one ghost pebble
#
# two types of ghost pebble 
# 1st type, mirror of cylinder boundary
# 2nd type, mirror of top and bot planes 
#
# however, in the new method, top and bot ghost pebbles are read from dem data.
#
# 1. top ghost pebbles
#top_gpebble_file_path = 'tsk1_top_ghost.dat' # touching pebbles
top_gpebble_file_path = 'tsk3_top_ghost.dat' # touching pebbles
top_gpebble_file = open(top_gpebble_file_path) 
number_of_gpebbles = int(top_gpebble_file.readline().split()[0])
for ipebble in range(0,number_of_gpebbles):
	line = top_gpebble_file.readline()
	p_x = float(line.split()[0])*scaling
	p_y = float(line.split()[1])*scaling
	p_z = float(line.split()[2])*scaling+0.6
	ifghost = True
	tag = 'TOP'
	random_pebbles.new_pebble(p_x,p_y,p_z,ifghost,tag)    # add new pebble
top_gpebble_file.close()

# 2. bot ghost pebbles
#bot_gpebble_file_path = 'tsk1_bot_ghost.dat' # touching pebbles
bot_gpebble_file_path = 'tsk3_bot_ghost.dat' # touching pebbles
bot_gpebble_file = open(bot_gpebble_file_path) 
number_of_gpebbles = int(bot_gpebble_file.readline().split()[0])
for ipebble in range(0,number_of_gpebbles):
	line = bot_gpebble_file.readline()
	p_x = float(line.split()[0])*scaling
	p_y = float(line.split()[1])*scaling
	p_z = float(line.split()[2])*scaling-0.6
	ifghost = True
	tag = 'BOT'
	random_pebbles.new_pebble(p_x,p_y,p_z,ifghost,tag)    # add new pebble
bot_gpebble_file.close()

number_of_pebbles = int(random_pebbles.number_of_pebbles())
print 'number of all (real+ghost)pebbles: ' + str(number_of_pebbles)
print '==========================================================='


# 3. cylinder sidewall ghost pebbles
number_of_pebbles = int(random_pebbles.number_of_pebbles())
for ipebble in range(0,number_of_pebbles):
	p_x = random_pebbles.pebble_xyz(ipebble)[0]
	p_y = random_pebbles.pebble_xyz(ipebble)[1]
	p_z = random_pebbles.pebble_xyz(ipebble)[2]
	p_r = (p_x**2.0 + p_y**2.0)**0.5
	d_r = cyl_radius - p_r

	if p_r > 0.5*pebble_diameter:
		if d_r <= pebble_diameter*ghost_pebble_tol_side: 
			# construct ghost pebble
			nv_x = p_x/p_r
			nv_y = p_y/p_r
			g_r = d_r + cyl_radius
			g_x = g_r*nv_x
			g_y = g_r*nv_y
			g_z = p_z
			ifghost = True
			tag = 'SW '
			random_pebbles.new_pebble(g_x,g_y,g_z,ifghost,tag)

#print 'all pebbles'
#print random_pebbles.all_pebbles_xyz()

number_of_pebbles = int(random_pebbles.number_of_pebbles())
print 'number of all (real+ghost)pebbles: ' + str(number_of_pebbles)
print '==========================================================='

#=================================================================================
# part 2
# construct Voronoi cell from pebble coordinates
# and extract information from it

vor = Voronoi(random_pebbles.all_pebbles_xyz())

vor_vert = voro_vertices('voro_vertices') # instance of class voro_vertices

for ivert in range(0,len(vor.vertices)):
	v_x = vor.vertices[ivert][0]
	v_y = vor.vertices[ivert][1]
	v_z = vor.vertices[ivert][2]
	vor_vert.new_vertice(v_x,v_y,v_z)
		
#extract   vor.ridge_vertices, voronoi ridge becomes face in 3d
vor_face = voro_faces('voro_faces')
for iface in range(0,len(vor.ridge_vertices)):
	vlist = list(vor.ridge_vertices[iface])
	plist = list(vor.ridge_points[iface])
	vor_face.new_face(vlist,plist)
	for v in vlist:
		vor_vert.v_to_f[v].append(vor_face.nf-1) # info from vertice to voro face
		
for iface in range(0,vor_face.nf):
	p1 = vor_face.f_to_p[iface][0]
	p1ghost = random_pebbles.p_ghost[p1]
	p2 = vor_face.f_to_p[iface][1]
	p2ghost = random_pebbles.p_ghost[p2]
	ifBothSideGhostPebble = p1ghost and p2ghost
	vor_face.ifghost[iface] = ifBothSideGhostPebble

vtkFileName = 'polygon_0.vtk'
dump_poly_vtk(vor_vert,vor_face,vtkFileName)


#tol = [0.05,0.15,0.175,0.05]
#tol = [0.05,0.06,0.0725,0.06,0.05] # tol smaller will reduce multiplicity

tol = [0.05,0.07,0.09,0.07,0.05]# tol smaller will reduce multiplicity
#tol = [0.05,0.08,0.1,0.05] # tol smaller will reduce multiplicity
nm = 5
minAngle = 20
minArea = 0.04

# update cylinder radius based on max pebble_to_center radius
#cyl_radius = mx_rr + 0.5

find_nbr_vertices(vor_vert,vor_face)
# a quicker way to merge vertice
# and store merge information, which will be used when vertices are moved

for im in range(0,nm):
	merge_vertices(vor_vert,tol[im])

vtkFileName = 'polygon_1.vtk'
print ' dumping '+ vtkFileName
dump_poly_vtk(vor_vert,vor_face,vtkFileName)

#==============================================================

for iface in range(0,vor_face.nf):# pack edge list in a face for all faces
	vor_face.find_edge_info(iface)
	vor_face.calculate_edge_length(iface,vor_vert)
	vor_face.clean_up_zero_length_edge(iface,vor_vert)

for iface in range(0,vor_face.nf):
	p1 = vor_face.f_to_p[iface][0]
	p1xyz = random_pebbles.pebble_xyz(p1)
	p2 = vor_face.f_to_p[iface][1]
	p2xyz = random_pebbles.pebble_xyz(p2)
	vor_face.add_f_pp_center(iface,p1xyz,p2xyz,vor_vert)
	vor_face.add_f_center(iface,vor_vert,random_pebbles,pebble_diameter)
	vor_face.calculate_vertice_angles(iface,vor_vert)

# fix sharp triangle
for iface in range(0,vor_face.nf):	
	vor_face.fix_sharp_triangle(iface,vor_vert,minAngle)

# update face info after sharp triangle is fixed
vtkFileName = 'polygon_2.vtk'
print ' dumping '+ vtkFileName
dump_poly_vtk(vor_vert,vor_face,vtkFileName)

for iface in range(0,vor_face.nf):
	vor_face.find_edge_info(iface)
	vor_face.calculate_edge_length(iface,vor_vert)
	vor_face.clean_up_zero_length_edge(iface,vor_vert)

for iface in range(0,vor_face.nf):
	p1 = vor_face.f_to_p[iface][0]
	p1xyz = random_pebbles.pebble_xyz(p1)
	p2 = vor_face.f_to_p[iface][1]
	p2xyz = random_pebbles.pebble_xyz(p2)
	vor_face.add_f_pp_center(iface,p1xyz,p2xyz,vor_vert)
	vor_face.add_f_center(iface,vor_vert,random_pebbles,pebble_diameter)
	vor_face.calculate_vertice_angles(iface,vor_vert)


## fix sharp triangle
## remove tiny face, ???
#for iface in range(0,vor_face.nf):
#	vor_face.remove_tiny_face(iface,vor_vert,minArea)
#
#for iface in range(0,vor_face.nf):
#	vor_face.find_edge_info(iface)
#	vor_face.calculate_edge_length(iface,vor_vert)
#	vor_face.clean_up_zero_length_edge(iface,vor_vert)
#
#	p1 = vor_face.f_to_p[iface][0]
#	p1xyz = random_pebbles.pebble_xyz(p1)
#	p2 = vor_face.f_to_p[iface][1]
#	p2xyz = random_pebbles.pebble_xyz(p2)
#	vor_face.add_f_pp_center(iface,p1xyz,p2xyz,vor_vert)
#	vor_face.add_f_center(iface,vor_vert,random_pebbles,pebble_diameter)
#	vor_face.calculate_vertice_angles(iface,vor_vert)
#
#for iface in range(0,vor_face.nf):	
#	vor_face.fix_sharp_triangle(iface,vor_vert,minAngle)
#
#for iface in range(0,vor_face.nf):
#	vor_face.find_edge_info(iface)
#	vor_face.calculate_edge_length(iface,vor_vert)
#	vor_face.clean_up_zero_length_edge(iface,vor_vert)
#	
#	p1 = vor_face.f_to_p[iface][0]
#	p1xyz = random_pebbles.pebble_xyz(p1)
#	p2 = vor_face.f_to_p[iface][1]
#	p2xyz = random_pebbles.pebble_xyz(p2)
#	vor_face.add_f_pp_center(iface,p1xyz,p2xyz,vor_vert)
#	vor_face.add_f_center(iface,vor_vert,random_pebbles,pebble_diameter)
#	vor_face.calculate_vertice_angles(iface,vor_vert)


# fix concave angle

for iface in range(0,vor_face.nf):
	vor_face.fix_concave_angle(iface,vor_vert)

for iface in range(0,vor_face.nf):
	vor_face.find_edge_info(iface)
	vor_face.calculate_edge_length(iface,vor_vert)
	vor_face.clean_up_zero_length_edge(iface,vor_vert)
	
	p1 = vor_face.f_to_p[iface][0]
	p1xyz = random_pebbles.pebble_xyz(p1)
	p2 = vor_face.f_to_p[iface][1]
	p2xyz = random_pebbles.pebble_xyz(p2)
	vor_face.add_f_pp_center(iface,p1xyz,p2xyz,vor_vert)
	vor_face.add_f_center(iface,vor_vert,random_pebbles,pebble_diameter)
	vor_face.calculate_vertice_angles(iface,vor_vert)

for iface in range(0,vor_face.nf):	
	vor_face.fix_sharp_triangle(iface,vor_vert,minAngle)

for iface in range(0,vor_face.nf):
	vor_face.find_edge_info(iface)
	vor_face.calculate_edge_length(iface,vor_vert)
	vor_face.clean_up_zero_length_edge(iface,vor_vert)
	
	p1 = vor_face.f_to_p[iface][0]
	p1xyz = random_pebbles.pebble_xyz(p1)
	p2 = vor_face.f_to_p[iface][1]
	p2xyz = random_pebbles.pebble_xyz(p2)
	vor_face.add_f_pp_center(iface,p1xyz,p2xyz,vor_vert)
	vor_face.add_f_center(iface,vor_vert,random_pebbles,pebble_diameter)
	vor_face.calculate_vertice_angles(iface,vor_vert)

# bend corner vertice to reduce non-right-hand elements
beta = 0.2  # distance of move
delta = 20.0 # angles for search
for iface in range(0,vor_face.nf):
	vor_face.bend_corner_vertice(iface,vor_vert,random_pebbles,beta,delta)

beta = 0.1  # distance of move
delta = 10.0 # angles for search
for iface in range(0,vor_face.nf):
	vor_face.bend_corner_vertice(iface,vor_vert,random_pebbles,beta,delta)
	
	
vtkFileName = 'polygon_3.vtk'
print ' dumping '+ vtkFileName
dump_poly_vtk(vor_vert,vor_face,vtkFileName)

#======================================================
# before start dividing the polygon, we detect which polygon shold be add chamfer.
# if pebble-pebble distance is too close for this polygon.
#
for iface in range(0,vor_face.nf):
	p1 = vor_face.f_to_p[iface][0]
	p1xyz = random_pebbles.pebble_xyz(p1)
	p2 = vor_face.f_to_p[iface][1]
	p2xyz = random_pebbles.pebble_xyz(p2)
	vor_face.detect_chamfer(iface,p1xyz,p2xyz,pebble_diameter,1.025) # if pebble-pebble distance is less than 1.02*diameter, then consider them touching

vtkFileName = 'polygon_with_chamfer.vtk'
print ' dumping '+ vtkFileName
dump_poly_vtk_for_chamfer(vor_vert,vor_face,vtkFileName)

# redistribute polygon to avoid contact with pebbles
# not only vertices, but also detect edges. if edge is too close to a pebble, move its vertice away from it.
for iface in range(0,vor_face.nf):
	vor_face.redistribute(iface,vor_vert,random_pebbles,pebble_diameter,r_chamfer)

for iface in range(0,vor_face.nf):
	#tag = 'BOT'
	#vor_face.project_to_plane(iface,vor_vert,random_pebbles,cyl_bot,tag)
	#tag = 'TOP'
	#vor_face.project_to_plane(iface,vor_vert,random_pebbles,cyl_top,tag)
	cyl_radius_no_bl = cyl_radius - 0.01
	vor_face.project_to_sidewall(iface,vor_vert,random_pebbles,cyl_radius_no_bl)

# recalculate edge length after redistribution
for iface in range(0,vor_face.nf):
	vor_face.find_edge_info(iface)
	vor_face.calculate_edge_length(iface,vor_vert)
	vor_face.clean_up_zero_length_edge(iface,vor_vert)
	
	p1 = vor_face.f_to_p[iface][0]
	p1xyz = random_pebbles.pebble_xyz(p1)
	p2 = vor_face.f_to_p[iface][1]
	p2xyz = random_pebbles.pebble_xyz(p2)
	vor_face.add_f_pp_center(iface,p1xyz,p2xyz,vor_vert)
	vor_face.add_f_center(iface,vor_vert,random_pebbles,pebble_diameter)
	vor_face.calculate_vertice_angles(iface,vor_vert)
	
vtkFileName = 'polygon_4.vtk'
print ' dumping '+ vtkFileName
dump_poly_vtk(vor_vert,vor_face,vtkFileName)


# add mid edge vertice to super long edge 
longEdge = 0.6*pebble_diameter #0.6*pebble_diameter
for iface in range(0,vor_face.nf):
	vor_face.add_mid_vertice(iface,vor_vert,longEdge)

for iface in range(0,vor_face.nf):
	vor_face.find_edge_info(iface)
	vor_face.calculate_edge_length(iface,vor_vert)
	vor_face.clean_up_zero_length_edge(iface,vor_vert)
	
	p1 = vor_face.f_to_p[iface][0]
	p1xyz = random_pebbles.pebble_xyz(p1)
	p2 = vor_face.f_to_p[iface][1]
	p2xyz = random_pebbles.pebble_xyz(p2)
	vor_face.add_f_pp_center(iface,p1xyz,p2xyz,vor_vert)
	vor_face.add_f_center(iface,vor_vert,random_pebbles,pebble_diameter)
	vor_face.calculate_vertice_angles(iface,vor_vert)

vtkFileName = 'polygon_5.vtk'
print ' dumping '+ vtkFileName
dump_poly_vtk(vor_vert,vor_face,vtkFileName)
# from now on, cannot move vertice anymore.
#
# fix concave quad:
for iface in range(0,vor_face.nf):
	concave_quad_fix2(iface,vor_face,vor_vert,140)
	#concave_quad_fix3(iface,vor_face,vor_vert,140)
	
for iface in range(0,vor_face.nf):
	vor_face.find_edge_info(iface)
	vor_face.calculate_edge_length(iface,vor_vert)
	vor_face.clean_up_zero_length_edge(iface,vor_vert)
	
	p1 = vor_face.f_to_p[iface][0]
	p1xyz = random_pebbles.pebble_xyz(p1)
	p2 = vor_face.f_to_p[iface][1]
	p2xyz = random_pebbles.pebble_xyz(p2)
	vor_face.add_f_pp_center(iface,p1xyz,p2xyz,vor_vert)
	vor_face.add_f_center(iface,vor_vert,random_pebbles,pebble_diameter)
	vor_face.calculate_vertice_angles(iface,vor_vert)
	
vtkFileName = 'polygon_6.vtk'
print ' dumping '+ vtkFileName
dump_poly_vtk(vor_vert,vor_face,vtkFileName)

maxAngle = 150
for iface in range(0,vor_face.nf):
	polygon_divide_new(iface,vor_face,vor_vert,maxAngle)

for iface in range(0,vor_face.nf):
	vor_face.find_edge_info(iface)
	vor_face.calculate_edge_length(iface,vor_vert)
	vor_face.clean_up_zero_length_edge(iface,vor_vert)

for iface in range(0,vor_face.nf):
	vor_face.find_edge_info(iface)
	vor_face.calculate_edge_length(iface,vor_vert)
	vor_face.add_f_center(iface,vor_vert,random_pebbles,pebble_diameter)
	vor_face.calculate_vertice_angles(iface,vor_vert)

vtkFileName = 'polygon_7.vtk'
print ' dumping '+ vtkFileName
dump_poly_vtk(vor_vert,vor_face,vtkFileName)

maxAngle = 150
for iface in range(0,vor_face.nf):
	polygon_divide_new(iface,vor_face,vor_vert,maxAngle)

for iface in range(0,vor_face.nf):
	vor_face.find_edge_info(iface)
	vor_face.calculate_edge_length(iface,vor_vert)
	vor_face.clean_up_zero_length_edge(iface,vor_vert)

for iface in range(0,vor_face.nf):
	vor_face.find_edge_info(iface)
	vor_face.calculate_edge_length(iface,vor_vert)
	vor_face.add_f_center(iface,vor_vert,random_pebbles,pebble_diameter)
	vor_face.calculate_vertice_angles(iface,vor_vert)

vtkFileName = 'polygon_8.vtk'
print ' dumping '+ vtkFileName
dump_poly_vtk(vor_vert,vor_face,vtkFileName)

#================


# hexagon divide
for iface in range(0,vor_face.nf):
	polygon_divide_for_hexagon(iface,vor_face,vor_vert)

for iface in range(0,vor_face.nf):
	vor_face.find_edge_info(iface)
	vor_face.calculate_edge_length(iface,vor_vert)
	vor_face.clean_up_zero_length_edge(iface,vor_vert)

for iface in range(0,vor_face.nf):
	vor_face.find_edge_info(iface)
	vor_face.calculate_edge_length(iface,vor_vert)
	vor_face.add_f_center(iface,vor_vert,random_pebbles,pebble_diameter)
	vor_face.calculate_vertice_angles(iface,vor_vert)

vtkFileName = 'polygon_9.vtk'
print ' dumping '+ vtkFileName
dump_poly_vtk(vor_vert,vor_face,vtkFileName)


## pentagon divide
#for iface in range(0,vor_face.nf):
#	polygon_divide_for_pentagon(iface,vor_face,vor_vert)
#	#concave_quad_fix2(iface,vor_face,vor_vert,150)
#
#for iface in range(0,vor_face.nf):
#	vor_face.find_edge_info(iface)
#	vor_face.calculate_edge_length(iface,vor_vert)
#	vor_face.clean_up_zero_length_edge(iface,vor_vert)
#
#for iface in range(0,vor_face.nf):
#	vor_face.find_edge_info(iface)
#	vor_face.calculate_edge_length(iface,vor_vert)
#	vor_face.add_f_center(iface,vor_vert,random_pebbles,pebble_diameter)
#	vor_face.calculate_vertice_angles(iface,vor_vert)


#=================================================================================
# new divide for octagon,heptagon,hexagon, pentagon
for iface in range(0,vor_face.nf):
	polygon_divide_for_octagon_new(iface,vor_face,vor_vert)
	polygon_divide_for_heptagon_new(iface,vor_face,vor_vert)
	polygon_divide_for_hexagon_new(iface,vor_face,vor_vert)
	#polygon_divide_for_pentagon_new(iface,vor_face,vor_vert)

for iface in range(0,vor_face.nf):
	vor_face.find_edge_info(iface)
	vor_face.calculate_edge_length(iface,vor_vert)
	vor_face.clean_up_zero_length_edge(iface,vor_vert)

for iface in range(0,vor_face.nf):
	vor_face.find_edge_info(iface)
	vor_face.calculate_edge_length(iface,vor_vert)
	vor_face.add_f_center(iface,vor_vert,random_pebbles,pebble_diameter)
	vor_face.calculate_vertice_angles(iface,vor_vert)
#=================================================================================

# final fix concave quad
for iface in range(0,vor_face.nf):
	#polygon_divide_for_pentagon(iface,vor_face,vor_vert)
	concave_quad_fix2(iface,vor_face,vor_vert,120)

for iface in range(0,vor_face.nf):
	vor_face.find_edge_info(iface)
	vor_face.calculate_edge_length(iface,vor_vert)
	vor_face.clean_up_zero_length_edge(iface,vor_vert)

for iface in range(0,vor_face.nf):
	vor_face.find_edge_info(iface)
	vor_face.calculate_edge_length(iface,vor_vert)
	vor_face.add_f_center(iface,vor_vert,random_pebbles,pebble_diameter)
	vor_face.calculate_vertice_angles(iface,vor_vert)

vtkFileName = 'polygon_10.vtk'
print ' dumping '+ vtkFileName
dump_poly_vtk(vor_vert,vor_face,vtkFileName)
#==================================================================================

# check, edge number of polygon, avoid chamfer polygon
#
ng8 = 0
n8 = 0
n7 = 0
n6 = 0
n5 = 0
nmax = 0
for iface in range(0,vor_face.nf):
	ifghost = vor_face.ifghost[iface]
	ifcollapse = vor_face.ifcollapse[iface]
	ifchamfer = vor_face.ifchamfer[iface]
	if (not ifghost) and (not ifcollapse) and (not ifchamfer):

		nedges = len(vor_face.f_to_v[iface])
	
		if nedges > nmax:
			nmax = nedges
		
		if nedges > 8: 
			ng8 = ng8+1
		if nedges == 8:
			n8 = n8+1
		if nedges == 7:
			n7 = n7+1
		if nedges == 6:
			n6 = n6+1
		if nedges == 5:
			n5 = n5+1

print 'ng8: ' + str(ng8)
print 'n8: ' + str(n8)
print 'n7: ' + str(n7)
print 'n6: ' + str(n6)
print 'n5: ' + str(n5)

#=================================================================================
# part 3
# split polygon into quads
print 'starting generating quads'

vor_quads = voro_quads('voro_quads')

for iface in range(0,vor_face.nf):
	# only do this if not ghost face and non-collapsed face
	ifghost = vor_face.ifghost[iface]
	ifcollapse = vor_face.ifcollapse[iface]
	ifchamfer = vor_face.ifchamfer[iface]
	if (not ifghost) and (not ifcollapse) and (not ifchamfer):

		pxyz = []
		p1 = vor_face.f_to_p[iface][0]
		px = random_pebbles.xyz[p1][0]
		py = random_pebbles.xyz[p1][1]
		pz = random_pebbles.xyz[p1][2]
		pxyz.append([px,py,pz])

		p2 = vor_face.f_to_p[iface][1]
		px = random_pebbles.xyz[p2][0]
		py = random_pebbles.xyz[p2][1]
		pz = random_pebbles.xyz[p2][2]
		pxyz.append([px,py,pz])
					
		p1tag = random_pebbles.tag[p1]
		p2tag =	random_pebbles.tag[p2]
	
		ptag = [p1tag,p2tag]
	
		if p1tag == '' and p2tag == '': 
			ifInternal = True
		else:
			ifInternal = False
		
		
		nedges = len(vor_face.f_to_v[iface])
		
		if nedges == 3:
		# if triangle
			tri = []
			for iv in  range(0,3):
				v = vor_face.f_to_v[iface][iv]
				trix = vor_vert.xyz[v][0]
				triy = vor_vert.xyz[v][1]
				triz = vor_vert.xyz[v][2]
				tri.append([trix,triy,triz])
			if_has_ghost_pebble = random_pebbles.p_ghost[p1] or random_pebbles.p_ghost[p2]
			quads = tri_to_quad1(tri,pxyz,vor_face.f_center[iface],pebble_diameter,ptag,random_pebbles)
			for iquad in range(0,3):
				qxyz = quads[iquad]
				vor_quads.new_quad(qxyz,p1,p2)
		elif nedges == 4:
		# if rectangle
			rec = []
			for iv in  range(0,4):
				v = vor_face.f_to_v[iface][iv]
				recx = vor_vert.xyz[v][0]
				recy = vor_vert.xyz[v][1]
				recz = vor_vert.xyz[v][2]
				rec.append([recx,recy,recz])	
			quads = rec_to_quad(rec,pxyz,vor_face.f_center[iface],pebble_diameter,ptag,random_pebbles)
				
			for iquad in range(0,4):
				qxyz = quads[iquad]
				vor_quads.new_quad(qxyz,p1,p2)
		
		else:
			
			maxAngle = max(vor_face.angles[iface])
			if maxAngle <= 120: #130:
				# new method. better than old methd
				face_center = vor_face.f_center[iface]
				for iedge in range(0,nedges):
					ie1 = iedge
					ie2 = (iedge+1)%nedges
					
					v1 = vor_face.f_to_e[iface][ie1][0]
					v2 = vor_face.f_to_e[iface][ie1][1]
					
					xyz1 = vor_vert.xyz[v1]
					xyz2 = vor_vert.xyz[v2]
		
					mid1 = line_split(xyz1,xyz2,0.5)
									
					v1 = vor_face.f_to_e[iface][ie2][0]
					v2 = vor_face.f_to_e[iface][ie2][1]
					
					xyz1 = vor_vert.xyz[v1]
					xyz2 = vor_vert.xyz[v2]
		
					mid2 = line_split(xyz1,xyz2,0.5)
					
					quad = [face_center,mid2,xyz1,mid1]	
					vor_quads.new_quad(quad,p1,p2)
				
			else:
				# old method. 
				for iedge in range(0,nedges):
					tri = []
					# 1st point of tri is always the face center
					trix= vor_face.f_center[iface][0]
					triy= vor_face.f_center[iface][1]
					triz= vor_face.f_center[iface][2]
					tri.append([trix,triy,triz])
					
					v1 = vor_face.f_to_e[iface][iedge][0]
					if v1 < 0:
						print 'EORROR: ghost face used'
					trix = vor_vert.xyz[v1][0]
					triy = vor_vert.xyz[v1][1]
					triz = vor_vert.xyz[v1][2]
					tri.append([trix,triy,triz])
					
					v2 = vor_face.f_to_e[iface][iedge][1]
					if v2 < 0:
						print 'EORROR: ghost face used'
					trix = vor_vert.xyz[v2][0]
					triy = vor_vert.xyz[v2][1]
					triz = vor_vert.xyz[v2][2]
					tri.append([trix,triy,triz])
		
					if_has_ghost_pebble = random_pebbles.p_ghost[p1] or random_pebbles.p_ghost[p2]
		
					fc = [(tri[0][0]+tri[1][0]+tri[2][0])/3.0,(tri[0][1]+tri[1][1]+tri[2][1])/3.0,(tri[0][2]+tri[1][2]+tri[2][2])/3.0]
					fc,moved = move_away_from_pebbles(fc,pxyz,pebble_diameter,1.01,ptag)
					
					# now tri and pxyz are packed
					quads = tri_to_quad2(tri,pxyz,fc,pebble_diameter,ifInternal,ptag,random_pebbles)
		
					# 3 quads are generated
					for iquad in range(0,3):
						qxyz = quads[iquad]
						vor_quads.new_quad(qxyz,p1,p2)
		
		
	if (not ifghost) and (not ifcollapse) and (ifchamfer): # special process to deal with polygon with chamfers

		p1 = vor_face.f_to_p[iface][0]
		p1xyz = random_pebbles.xyz[p1]
		p2 = vor_face.f_to_p[iface][1]
		p2xyz = random_pebbles.xyz[p2]
		pxyz = [p1xyz,p2xyz]
		
		p1tag = random_pebbles.tag[p1]
		p2tag =	random_pebbles.tag[p2]
		ptag = [p1tag,p2tag]
	
		if p1tag == '' and p2tag == '': 
			ifInternal = True
		else:
			ifInternal = False
		
		
		nedges = len(vor_face.f_to_v[iface])
		
		for iedge in range(0,nedges):
			v1 = vor_face.f_to_e[iface][iedge][0]
			v1xyz = vor_vert.xyz[v1]
			v2 = vor_face.f_to_e[iface][iedge][1]
			v2xyz = vor_vert.xyz[v2]
			
			pmiddle = line_split(pxyz[0],pxyz[1],0.5)
			
			if p1tag == 'SW ' or p2tag == 'SW ':
				# if sidewall pebble,  project to siddewall
				# this is the pebble-pebble center
				pmiddle = proj_vertice_to_cyl(pmiddle,cyl_radius)
			
			
			d1 = distance(v1xyz,pmiddle)
			d2 = distance(v2xyz,pmiddle)
			r_c = r_chamfer*(pebble_diameter/2.0)
			vc1xyz =  line_split(pmiddle,v1xyz,r_c/d1)
			vc2xyz =  line_split(pmiddle,v2xyz,r_c/d2)
			
			vc1xyz,moved = move_away_from_pebbles(vc1xyz,pxyz,pebble_diameter,1.005,ptag)
			vc2xyz,moved = move_away_from_pebbles(vc2xyz,pxyz,pebble_diameter,1.005,ptag)
			
			vec1 = vector_minus(vc1xyz,pmiddle)
			vec2 = vector_minus(vc2xyz,pmiddle)
			vec3 = line_split(vec1,vec2,0.5)
			vec3norm= np.linalg.norm(vec3)
			vec3 = [vec3[0]/vec3norm,vec3[1]/vec3norm,vec3[2]/vec3norm]
			
			#vcmxyz = [pmiddle[0]+vec3[0]*r_c,pmiddle[1]+vec3[1]*r_c,pmiddle[2]+vec3[2]*r_c]
						
			m12 = line_split(v1xyz,v2xyz,0.5)
			d_m12 = distance(m12,pmiddle)
			
			if d_m12 <= (r_c +0.04): r_c = d_m12 - 0.04
			
			vcmxyz =  line_split(pmiddle,m12,r_c/d_m12)
			vcmxyz,moved = move_away_from_pebbles(vcmxyz,pxyz,pebble_diameter,1.005,ptag)
			
			
			
			rec = [v1xyz,vc1xyz,vc2xyz,v2xyz]
			
			
			#quads = rec_to_quad_for_chamfer(rec,vcmxyz,pxyz,pebble_diameter,ifInternal,ptag,random_pebbles)
			#	
			#for iquad in range(0,4):
			#	qxyz = quads[iquad]
			#	vor_quads.new_quad(qxyz,p1,p2)
			#	if iquad == 1:
			#		vor_quads.edge_tag[vor_quads.nquads-1][0] = 'C  '
			#	if iquad == 2:
			#		vor_quads.edge_tag[vor_quads.nquads-1][3] = 'C  '
			
			
			#quads = rec_to_quad_for_chamfer2(rec,vcmxyz,pxyz,pebble_diameter,ifInternal,ptag,random_pebbles)
			#
			## add linear flag to near chamfer element
			#for iquad in range(0,8):
			#	qxyz = quads[iquad]
			#	vor_quads.new_quad(qxyz,p1,p2)
			#	if iquad == 2:
			#		vor_quads.edge_tag[vor_quads.nquads-1][0] = 'C  '
			#		vor_quads.linearFlag[vor_quads.nquads-1] = True
			#	if iquad == 4:
			#		vor_quads.edge_tag[vor_quads.nquads-1][3] = 'C  '
			#		vor_quads.linearFlag[vor_quads.nquads-1] = True
			
			quads = rec_to_quad_for_chamfer3(rec,vcmxyz,pxyz,pebble_diameter,ifInternal,ptag,random_pebbles)
			
			# add linear flag to near chamfer element
			for iquad in range(0,10):
				qxyz = quads[iquad]
				vor_quads.new_quad(qxyz,p1,p2)
				if iquad == 2:
					vor_quads.edge_tag[vor_quads.nquads-1][0] = 'C  '
					vor_quads.linearFlag[vor_quads.nquads-1] = True
				if iquad == 5:
					vor_quads.edge_tag[vor_quads.nquads-1][3] = 'C  '
					vor_quads.linearFlag[vor_quads.nquads-1] = True
			
			
#vtkFileName = 'voro_quad.vtk'
#dump_quad_vtk(vor_quads,vtkFileName)

bot_quads = voro_quads('bot_quads')			
top_quads = voro_quads('top_quads')
vor_hex8s= voro_hex8s('voro_hex8s')

#
# 1. store the list of top quads that are on the same pebble
# 2. fix concave top quads to avoid non right hand element
# 

voro_quad_to_four_quads = False

if voro_quad_to_four_quads:
	for iquad in range (0,vor_quads.nquads):
		#vor_quads.xyz[iquad]
		p1 = vor_quads.quad_to_p[iquad][0]
		p1xyz = random_pebbles.pebble_xyz(p1)
		p2 = vor_quads.quad_to_p[iquad][1]
		p2xyz = random_pebbles.pebble_xyz(p2)
		ppxyz = [p1xyz,p2xyz]
		
		p1tag = random_pebbles.tag[p1]
		p2tag =	random_pebbles.tag[p2]
		ptag = [p1tag,p2tag]
		
		# split on voro quad to four bot quads
		mid = line_split(p1xyz,p2xyz,0.5)
		if 'C  ' in vor_quads.edge_tag[iquad]:
			four_quads,four_edge_tags = quad_to_4quads_for_chamfer(vor_quads.xyz[iquad],vor_quads.edge_tag[iquad],ppxyz,pebble_diameter,r_chamfer,ptag)
		else:
			four_quads,four_edge_tags = quad_to_4quads(vor_quads.xyz[iquad],vor_quads.edge_tag[iquad],ppxyz,pebble_diameter,ptag)
			
			
		for i in range(0,4):
			quad_bot = four_quads[i]
			bot_quads.new_quad(quad_bot,vor_quads.quad_to_p[iquad][0],vor_quads.quad_to_p[iquad][1])
			bot_quads.edge_tag[bot_quads.nquads-1] = four_edge_tags[i]
			bot_quads.add_linear_mid(bot_quads.nquads-1,ppxyz,pebble_diameter,ptag,random_pebbles)
			
			
			for ipebble in range(0,2):
				pebble_id = vor_quads.quad_to_p[iquad][ipebble]
				ifghost = random_pebbles.p_ghost[pebble_id]
				if not ifghost:
					# only construct hex for non-ghost pebble
					pxyz = random_pebbles.xyz[pebble_id]
					quad_top = project_to_pebble(quad_bot,pxyz,pebble_diameter)
					top_quads.new_quad(quad_top,pebble_id,pebble_id)
					top_quads.edge_tag[top_quads.nquads-1] = four_edge_tags[i]
					top_quads.add_linear_mid(top_quads.nquads-1,ppxyz,pebble_diameter,ptag,random_pebbles)

					qbot_id = bot_quads.nquads-1
					qtop_id = top_quads.nquads-1
					vor_hex8s.new_hex8(qbot_id,qtop_id,pebble_id)
					
					
					#
					# adjust top quad mid so it won't touch bot quad mid
					#
					for imid in range(0,4):
						bot_m = bot_quads.mid_xyz[qbot_id][imid]
						top_m = top_quads.mid_xyz[qtop_id][imid]
					
						r_bot_m = distance(pxyz,bot_m)
						r_top_m = distance(pxyz,top_m)
					
						nv = vector_minus(bot_m,pxyz)
						nv[0] = nv[0]/r_bot_m
						nv[1] = nv[1]/r_bot_m
						nv[2] = nv[2]/r_bot_m
					
						iv1 = imid
						iv2 = imid + 1
						if iv2 > 3: iv2 = 0
					
						pradius = pebble_diameter/2.0
						if r_bot_m > (pradius+0.015):
							r_top_m = pradius
							new_top_m = [0,0,0]
							new_top_m[0] = pxyz[0]+r_top_m*nv[0]
							new_top_m[1] = pxyz[1]+r_top_m*nv[1]
							new_top_m[2] = pxyz[2]+r_top_m*nv[2]
			
							top_quads.mid_xyz[qtop_id][imid] = new_top_m
						else:
							r_top_m = r_bot_m - 0.015
						
							new_top_m = [0,0,0]
							new_top_m[0] = pxyz[0]+r_top_m*nv[0]
							new_top_m[1] = pxyz[1]+r_top_m*nv[1]
							new_top_m[2] = pxyz[2]+r_top_m*nv[2]
						
							top_quads.mid_xyz[qtop_id][imid] = new_top_m
					
else:
	for iquad in range (0,vor_quads.nquads):
		#vor_quads.xyz[iquad]
		
		linearFlag = vor_quads.linearFlag[iquad]
		
		p1 = vor_quads.quad_to_p[iquad][0]
		p1xyz = random_pebbles.pebble_xyz(p1)
		p2 = vor_quads.quad_to_p[iquad][1]
		p2xyz = random_pebbles.pebble_xyz(p2)
		ppxyz = [p1xyz,p2xyz]
		mid = line_split(p1xyz,p2xyz,0.5)
		
		p1tag = random_pebbles.tag[p1]
		p2tag =	random_pebbles.tag[p2]
		ptag = [p1tag,p2tag]

		bot_quads.new_quad(vor_quads.xyz[iquad],vor_quads.quad_to_p[iquad][0],vor_quads.quad_to_p[iquad][1])
		bot_quads.edge_tag[bot_quads.nquads-1] = vor_quads.edge_tag[iquad]
		
		# bot quad is linear.
		# however, if mid point is too close to pebble, move it.
		bot_quads.add_linear_mid(bot_quads.nquads-1)
		
		
		for ipebble in range(0,2):
			pebble_id = vor_quads.quad_to_p[iquad][ipebble]
			ifghost = random_pebbles.p_ghost[pebble_id]
			if not ifghost:
				# only construct hex for non-ghost pebble
				pxyz = random_pebbles.xyz[pebble_id]
				
				tol = 0.01 #0.0075
				quad_top = project_to_pebble(bot_quads.xyz[bot_quads.nquads-1],pxyz,pebble_diameter,tol)
				
				top_quads.new_quad(quad_top,pebble_id,pebble_id)
				top_quads.edge_tag[top_quads.nquads-1] = vor_quads.edge_tag[iquad]
				
				#top_quads.add_curve_mid(top_quads.nquads-1,pxyz,pebble_diameter)
				
				# use linear now, but will be reset later
				top_quads.add_linear_mid(top_quads.nquads-1)

					
					# adjust this top_quads for chamfer quads.
					# cannot direct project to pebble surface.
					# should just project to 
					#if 'C  ' in four_edge_tags:
						#adjusted_top_quad = adjust_top_quad_for_chamfer(quad_bot,quad_top,four_edge_tags,mid,pxyz)
						#top_quads.xyz[top_quads.nquads-1] = adjusted_top_quad
					
				qbot_id = bot_quads.nquads-1
				qtop_id = top_quads.nquads-1
				vor_hex8s.new_hex8(qbot_id,qtop_id,pebble_id)
				
				random_pebbles.top_quads[pebble_id].append(qtop_id) # store topquads on the same pebble
				top_quads.quad_to_quad[qtop_id].append(qbot_id)
				bot_quads.quad_to_quad[qbot_id].append(qtop_id)
		
				if (not linearFlag):
					bot_m = bot_quads.mid_xyz[qbot_id]
					top_quads.mid_xyz[qtop_id] = project_to_pebble(bot_m,pxyz,pebble_diameter,tol)


for iquad in range (0,bot_quads.nquads):
	p1 = bot_quads.quad_to_p[iquad][0]
	p2 = bot_quads.quad_to_p[iquad][1]
	p1tag = random_pebbles.tag[p1]
	p2tag =	random_pebbles.tag[p2]
	qtag = 'E  '
	if p1tag != '':  qtag = p1tag
	if p2tag != '':  qtag = p2tag
	bot_quads.set_tag(iquad,qtag)

	
#vtkFileName = 'top_quad.vtk'
#dump_quad_2nd_vtk(top_quads,vtkFileName)
			
#vtkFileName = 'bot_quad.vtk'
#dump_quad_2nd_vtk(bot_quads,vtkFileName)

# now, fix concave top quads
# this is the 'ultimate solution' to fix non-rith-hand elements
#fix_concave_top_quads(top_quads,random_pebbles)

#vtkFileName = 'top_quad_fixed.vtk'
#dump_quad_2nd_vtk(top_quads,vtkFileName)


print 'generating hex20 elements '

# from voro_hex8 convert to nek hex20 elements
nlayers = 2 							# number of layers, including boundary layer
bl_ratio = 0.1      # boundary layer thickness
bl_thickness = 0.015
# boundary layer must on pebble surface.
# if this voro_hex8 is on cylinder side, need boundary layer as well
# not consider this right now.

hex20 = nek_hex20s('nek_hex20s')

# for internal hex20 elements.
for ihex8 in range(0,vor_hex8s.nhex8s):
	qbot_id = vor_hex8s.quad_bot[ihex8]
	qtop_id = vor_hex8s.quad_top[ihex8]
	
	ipebble = top_quads.quad_to_p[qtop_id][0]
	pxyz = random_pebbles.xyz[ipebble]
		
	pz_ratio = (pxyz[2]-cyl_bot)/(cyl_top-cyl_bot)
	add_addition_layer = True
	add_addition_layer_to_entrance_and_exit_pebbles = True

	# construct boundary layer
	quad_m = []
	quad_mm = []
	for iv in range(0,4):
		xyz1 = top_quads.xyz[qtop_id][iv]
		xyz2 = bot_quads.xyz[qbot_id][iv]
		#new_vert = line_split(xyz1,xyz2,bl_ratio)
		new_vert = line_split_for_boundary_layer(xyz1,xyz2,bl_thickness)
		quad_m.append(new_vert)
    
	v8 = []
	for iv in range(0,4):
		v8.append(quad_m[iv])
	for iv in range(0,4):
		v8.append(top_quads.xyz[qtop_id][iv])
		
	e12 = []
	for ie in range(0,4):
		xyz1 = top_quads.mid_xyz[qtop_id][ie]
		xyz2 = bot_quads.mid_xyz[qbot_id][ie]
		#new_vert = line_split(xyz1,xyz2,bl_ratio)
		new_vert = line_split_for_boundary_layer(xyz1,xyz2,bl_thickness)
		quad_mm.append(new_vert)
		e12.append(new_vert)
	for ie in range(0,4):
		e12.append(top_quads.mid_xyz[qtop_id][ie])
	for ie in range(0,4):
		mx = (quad_m[ie][0]+top_quads.xyz[qtop_id][ie][0])/2.0
		my = (quad_m[ie][1]+top_quads.xyz[qtop_id][ie][1])/2.0
		mz = (quad_m[ie][2]+top_quads.xyz[qtop_id][ie][2])/2.0
		e12.append([mx,my,mz])
    
	s6 = []
	for iface in range(0,4):
		s6.append(top_quads.edge_tag[qtop_id][iface])
	s6.append('E  ')
	s6.append('PW ')
	
	hex20.new_hex(v8,e12,s6)

	#
	if add_addition_layer:
		if (pz_ratio>0.0) and (pz_ratio<1.0):
			#quad_m,quad_mm is known.
			#quad_m1, quad_mm1  are mid of quad_m and quat_bot
			quad_m1 = []
			quad_mm1 = []
			for iv in range(0,4):
				xyz1 = quad_m[iv]
				xyz2 = bot_quads.xyz[qbot_id][iv]
				new_vert = line_split(xyz1,xyz2,0.5)
				quad_m1.append(new_vert)
			
			v8 = []
			for iv in range(0,4):
				v8.append(quad_m1[iv])
			for iv in range(0,4):
				v8.append(quad_m[iv])
			
			e12 = []
			for ie in range(0,4):
				xyz1 = quad_mm[ie]
				xyz2 = bot_quads.mid_xyz[qbot_id][ie]
				new_vert = line_split(xyz1,xyz2,0.5)
				quad_mm1.append(new_vert)
				e12.append(new_vert)
			for ie in range(0,4):
				e12.append(quad_mm[ie])
			for ie in range(0,4):
				mxyz = line_split(quad_m[ie],quad_m1[ie],0.5)
				e12.append(mxyz)
		
			# no need to update s6, but s6[5] change to E
			s6 = []
			for iface in range(0,4):
				s6.append(bot_quads.edge_tag[qbot_id][iface])
			s6.append('E  ')
			s6.append('E  ')
			
			hex20.new_hex(v8,e12,s6)
			quad_m = quad_m1
			quad_mm = quad_mm1
		else:
			# add additional 3 layers to entrance and exit pebbles
		
			# 1st layer
			#quad_m,quad_mm is known.
			#quad_m1, quad_mm1 
			quad_m1 = []
			quad_mm1 = []
			for iv in range(0,4):
				xyz1 = quad_m[iv]
				xyz2 = bot_quads.xyz[qbot_id][iv]
				new_vert = line_split(xyz1,xyz2,0.25)
				quad_m1.append(new_vert)
			
			v8 = []
			for iv in range(0,4):
				v8.append(quad_m1[iv])
			for iv in range(0,4):
				v8.append(quad_m[iv])
			
			e12 = []
			for ie in range(0,4):
				xyz1 = quad_mm[ie]
				xyz2 = bot_quads.mid_xyz[qbot_id][ie]
				new_vert = line_split(xyz1,xyz2,0.25)
				quad_mm1.append(new_vert)
				e12.append(new_vert)
			for ie in range(0,4):
				e12.append(quad_mm[ie])
			for ie in range(0,4):
				mxyz = line_split(quad_m[ie],quad_m1[ie],0.5)
				e12.append(mxyz)
		
			# no need to update s6, but s6[5] change to E
			s6 = []
			for iface in range(0,4):
				s6.append(bot_quads.edge_tag[qbot_id][iface])
			s6.append('E  ')
			s6.append('E  ')
			
			hex20.new_hex(v8,e12,s6)
			quad_m = quad_m1
			quad_mm = quad_mm1
			
			# 2nd layer
			#quad_m,quad_mm is known.
			#quad_m1, quad_mm1 
			quad_m1 = []
			quad_mm1 = []
			for iv in range(0,4):
				xyz1 = quad_m[iv]
				xyz2 = bot_quads.xyz[qbot_id][iv]
				new_vert = line_split(xyz1,xyz2,0.333)
				quad_m1.append(new_vert)
			
			v8 = []
			for iv in range(0,4):
				v8.append(quad_m1[iv])
			for iv in range(0,4):
				v8.append(quad_m[iv])
			
			e12 = []
			for ie in range(0,4):
				xyz1 = quad_mm[ie]
				xyz2 = bot_quads.mid_xyz[qbot_id][ie]
				new_vert = line_split(xyz1,xyz2,0.333)
				quad_mm1.append(new_vert)
				e12.append(new_vert)
			for ie in range(0,4):
				e12.append(quad_mm[ie])
			for ie in range(0,4):
				mxyz = line_split(quad_m[ie],quad_m1[ie],0.5)
				e12.append(mxyz)
		
			# no need to update s6, but s6[5] change to E
			s6 = []
			for iface in range(0,4):
				s6.append(bot_quads.edge_tag[qbot_id][iface])
			s6.append('E  ')
			s6.append('E  ')
			
			hex20.new_hex(v8,e12,s6)
			quad_m = quad_m1
			quad_mm = quad_mm1
			
			# 3nd layer
			#quad_m,quad_mm is known.
			#quad_m1, quad_mm1 
			quad_m1 = []
			quad_mm1 = []
			for iv in range(0,4):
				xyz1 = quad_m[iv]
				xyz2 = bot_quads.xyz[qbot_id][iv]
				new_vert = line_split(xyz1,xyz2,0.5)
				quad_m1.append(new_vert)
			
			v8 = []
			for iv in range(0,4):
				v8.append(quad_m1[iv])
			for iv in range(0,4):
				v8.append(quad_m[iv])
			
			e12 = []
			for ie in range(0,4):
				xyz1 = quad_mm[ie]
				xyz2 = bot_quads.mid_xyz[qbot_id][ie]
				new_vert = line_split(xyz1,xyz2,0.5)
				quad_mm1.append(new_vert)
				e12.append(new_vert)
			for ie in range(0,4):
				e12.append(quad_mm[ie])
			for ie in range(0,4):
				mxyz = line_split(quad_m[ie],quad_m1[ie],0.5)
				e12.append(mxyz)
		
			# no need to update s6, but s6[5] change to E
			s6 = []
			for iface in range(0,4):
				s6.append(bot_quads.edge_tag[qbot_id][iface])
			s6.append('E  ')
			s6.append('E  ')
			
			hex20.new_hex(v8,e12,s6)
			quad_m = quad_m1
			quad_mm = quad_mm1
	
	
	# construct bottom layer
	v8 = []
	for iv in range(0,4):
		v8.append(bot_quads.xyz[qbot_id][iv])
	for iv in range(0,4):
		v8.append(quad_m[iv])
	
	e12 = []
	for ie in range(0,4):
		e12.append(bot_quads.mid_xyz[qbot_id][ie])
	for ie in range(0,4):
		e12.append(quad_mm[ie])
	for ie in range(0,4):
		mx = (quad_m[ie][0]+bot_quads.xyz[qbot_id][ie][0])/2.0
		my = (quad_m[ie][1]+bot_quads.xyz[qbot_id][ie][1])/2.0
		mz = (quad_m[ie][2]+bot_quads.xyz[qbot_id][ie][2])/2.0
		e12.append([mx,my,mz])
	
	s6 = []
	for iface in range(0,4):
		s6.append(bot_quads.edge_tag[qbot_id][iface])
	s6.append('E  ')
	s6.append('E  ')
	
	hex20.new_hex(v8,e12,s6)

#hex20.fix_non_right_hand_elements()


#===========================================================================
# need new method to extrude up and down elements
# 
#
#
# only project to 50%
target_dn_plane_z = -0.08*scaling
target_up_plaen_z = 1.96*scaling

nlayers_up = 30
nlayers_dn = 20

delta_sw = 0.01

# to help fast extact dn and up corner quads
#
#
# 
top_side_quads = voro_quads('top_side_quads')
bot_side_quads = voro_quads('bot_side_quads')

for iquad in range (0,bot_quads.nquads):
	tag = bot_quads.tag[iquad]
	
	# for top layers
	if tag == 'TOP':
		
		for ilayer in range (0,nlayers_up):
			s6 = []
			for iface in range(0,6):
				s6.append('E  ')
			
			v8 = []
			for iv in range(0,4):
				zbot = bot_quads.xyz[iquad][iv][2]
				delta_up = (target_up_plaen_z-zbot)/nlayers_up
				
				offset = [0,0,delta_up*ilayer]
				vert_new =offset_vert(bot_quads.xyz[iquad][iv],offset)
				v8.append(vert_new)
			for iv in range(0,4):
				zbot = bot_quads.xyz[iquad][iv][2]
				delta_up = (target_up_plaen_z-zbot)/nlayers_up
				
				offset = [0,0,delta_up*(ilayer+1)]
				vert_new =offset_vert(bot_quads.xyz[iquad][iv],offset)
				v8.append(vert_new)
			
			e12 = v8_to_e12(v8)	
			if ilayer == (nlayers_up-1): s6[5] = 'TOP'
			hex20.new_hex(v8,e12,s6)
			
			ihex = hex20.nhex-1
			for iface in range(0,4):
				vc,fn,vquad,equad,edge_tag = hex20.return_face_info(ihex,iface)
				cxyz = [0,0,vc[2]]
				
				ifsidewall = False
				# judge vc an fn, determine if it is on sidewall
				rr =  math.sqrt(vc[0]**2.0 + vc[1]**2.0)
				vec1 = vector_minus(vc,cxyz)
				angle = angles_two_vectors(vec1,fn)
				
				if (rr>(cyl_radius-0.02)) and ((angle<15)or(angle>165)): ifsidewall = True
				
				# sidewall
				if (ifsidewall):
					# add to top_side_quads
					top_side_quads.new_quad(vquad,0,0)
					top_side_quads.edge_tag[top_side_quads.nquads-1] = edge_tag
					top_side_quads.add_mid(top_side_quads.nquads-1,equad)
				
			
	# for bot layers
	if tag == 'BOT':
	
		for ilayer in range (0,nlayers_dn):
			s6 = []
			for iface in range(0,6):
				s6.append('E  ')
			
			v8 = []
			for iv in range(0,4):
				zbot = bot_quads.xyz[iquad][iv][2]
				delta_dn = (zbot-target_dn_plane_z)/nlayers_dn
			
				offset = [0,0,-delta_dn*ilayer]
				vert_new =offset_vert(bot_quads.xyz[iquad][iv],offset)
				v8.append(vert_new)
			for iv in range(0,4):
				zbot = bot_quads.xyz[iquad][iv][2]
				delta_dn = (zbot-target_dn_plane_z)/nlayers_dn
				
				offset = [0,0,-delta_dn*(ilayer+1)]
				vert_new =offset_vert(bot_quads.xyz[iquad][iv],offset)
				v8.append(vert_new)
			
			e12 = v8_to_e12(v8)	
			if ilayer == (nlayers_dn-1): s6[5] = 'BOT'
			hex20.new_hex(v8,e12,s6)
			
			ihex = hex20.nhex-1
			for iface in range(0,4):
				vc,fn,vquad,equad,edge_tag = hex20.return_face_info(ihex,iface)
				cxyz = [0,0,vc[2]]
				
				ifsidewall = False
				# judge vc an fn, determine if it is on sidewall
				rr =  math.sqrt(vc[0]**2.0 + vc[1]**2.0)
				vec1 = vector_minus(vc,cxyz)
				angle = angles_two_vectors(vec1,fn)
				
				if (rr>(cyl_radius-0.02)) and ((angle<15)or(angle>165)): ifsidewall = True
				
				# sidewall
				if (ifsidewall):
					# add to bot_side_quads
					bot_side_quads.new_quad(vquad,0,0)
					bot_side_quads.edge_tag[bot_side_quads.nquads-1] = edge_tag
					bot_side_quads.add_mid(bot_side_quads.nquads-1,equad)
	
	# for side wall
	if tag == 'SW ':
		
		s6 = []
		for iface in range(0,4):
			s6.append(bot_quads.edge_tag[iquad][iface])
		s6.append('E  ')
		s6.append('SW ')

		cxyz = [0,0,0]
		v8 = []
		for iv in range(0,4):
			v8.append(bot_quads.xyz[iquad][iv])
		for iv in range(0,4):
			vbot = bot_quads.xyz[iquad][iv]
			cxyz[2] = vbot[2]
			vert_new = line_split_for_boundary_layer(vbot,cxyz,-delta_sw)
			v8.append(vert_new)
		
		e12 = v8_to_e12(v8)
		
		hex20.new_hex(v8,e12,s6)
		
vtkFileName = 'top_side_quads.vtk'
dump_quad_2nd_vtk(top_side_quads,vtkFileName)
vtkFileName = 'bot_side_quads.vtk'
dump_quad_2nd_vtk(bot_side_quads,vtkFileName)	
	
# based on top_side_quads and bot_side_quads
# generate sw boundadry layers

for iquad in range (0,top_side_quads.nquads):
	s6 = []
	for iface in range(0,4):
		s6.append(top_side_quads.edge_tag[iquad][iface])
	s6.append('E  ')
	s6.append('SW ')

	cxyz = [0,0,0]
	v8 = []
	for iv in range(0,4):
		v8.append(top_side_quads.xyz[iquad][iv])
	for iv in range(0,4):
		vbot = top_side_quads.xyz[iquad][iv]
		cxyz[2] = vbot[2]
		vert_new = line_split_for_boundary_layer(vbot,cxyz,-delta_sw)
		v8.append(vert_new)
		
	e12 = v8_to_e12(v8)
	
	hex20.new_hex(v8,e12,s6)
	
for iquad in range (0,bot_side_quads.nquads):
	s6 = []
	for iface in range(0,4):
		s6.append(bot_side_quads.edge_tag[iquad][iface])
	s6.append('E  ')
	s6.append('SW ')

	cxyz = [0,0,0]
	v8 = []
	for iv in range(0,4):
		v8.append(bot_side_quads.xyz[iquad][iv])
	for iv in range(0,4):
		vbot = bot_side_quads.xyz[iquad][iv]
		cxyz[2] = vbot[2]
		vert_new = line_split_for_boundary_layer(vbot,cxyz,-delta_sw)
		v8.append(vert_new)
		
	e12 = v8_to_e12(v8)
	
	hex20.new_hex(v8,e12,s6)

hex20.fix_non_right_hand_elements()
dump_hex_2nd_vtk_nnrh(hex20,'hex20_nnrh.vtk')

#hex20.check_max_aspect_ratio()
#dump_hex_2nd_vtk_nnrh(hex20,'hex20_har.vtk')

#hex20.check_non_right_hand_elements_use_nek_method()
#dump_hex_2nd_vtk_nnrh(hex20,'hex20_nnrh.vtk')

# parallel version for check mesh

#np = 32
#ftag = 'har'
#check_max_aspect_ratio_parallel(hex20,np,ftag)
#dump_hex_2nd_vtk_nnrh(hex20,'hex20_har.vtk')

#ftag = 'nnrh'
#check_non_right_hand_elements_use_nek_method_parallel(hex20,np,ftag)
#dump_hex_2nd_vtk_nnrh(hex20,'hex20_nnrh.vtk')


#dump_hex_vtk(hex20,'hex8.vtk')
#dump_hex_2nd_vtk(hex20,'hex20.vtk')
#dump_hex_2nd_vtk_nnrh(hex20,'hex20_nnrh.vtk')

#=================================================================================
# part 4
# generate rea file based on hex20 information.
# no edge curvature right now

# explicit fix neg-jacobian element by element number provided by nek.
# fixing it by linearize this element.
# some mesh does not need this step.....
#fix_negative_jacobian_elements_from_nek(hex20)

endtime = time.time()
print 'time used:' + str(endtime-starttime)

baseRea = 'base.rea'
newReaFile = 'kp3b.rea'
alpha = 1.0
#hex20.scaling(alpha)
write_rea(baseRea,newReaFile,hex20)

#==============================================================================
endtime = time.time()
print 'time used:' + str(endtime-starttime)
