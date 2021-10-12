# python module
import numpy as np
import math
from scipy.spatial import Voronoi, voronoi_plot_2d
import multiprocessing
import time
#=================================================================================
# user defined files
from data_class import *
from merge_vertices import *
from splitting_subroutines import *
from polygon_divide import *
from write_rea import *
from dump_vtk import *
from fix_elements import *
from polygons_to_hexs import *

starttime = time.time()
#=================================================================================
# all input here:

if_top_bot_ghost_pebble_from_files = True

pebble_file_path = 'pebbles.dat'
top_gpebble_file_path = 'top_ghost.dat'
bot_gpebble_file_path = 'bot_ghost.dat'

nprocs = 16 # number of parallel processors

#=================================================================================
# only touch below when you know what you are doing !!! 
ghost_pebble_tol_side = 1.5
ghost_pebble_tol_top_bot = 2.0

r_chamfer = 0.1 # chamfer radius relative to pebble radius

nm = 5          # number of vertice merge iteration
tol = [0.05,0.08,0.1,0.08,0.05] # tolerance of each vertice merge iteration
#=================================================================================

geo = []
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

geo.append(cyl_radius)
geo.append(cyl_bot)
geo.append(cyl_top)
geo.append(r_chamfer)
geo.append(scaling)

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

# construct top and bot ghost pebbles

if if_top_bot_ghost_pebble_from_files:	
	# 1. top ghost pebbles
	top_gpebble_file = open(top_gpebble_file_path) 
	number_of_gpebbles = int(top_gpebble_file.readline().split()[0])
	for ipebble in range(0,number_of_gpebbles):
		line = top_gpebble_file.readline()
		p_x = float(line.split()[0])*scaling
		p_y = float(line.split()[1])*scaling
		p_z = float(line.split()[2])*scaling+0.5
		ifghost = True
		tag = 'TOP'
		random_pebbles.new_pebble(p_x,p_y,p_z,ifghost,tag)    # add new pebble
	top_gpebble_file.close()

	# 2. bot ghost pebbles
	bot_gpebble_file = open(bot_gpebble_file_path) 
	number_of_gpebbles = int(bot_gpebble_file.readline().split()[0])
	for ipebble in range(0,number_of_gpebbles):
		line = bot_gpebble_file.readline()
		p_x = float(line.split()[0])*scaling
		p_y = float(line.split()[1])*scaling
		p_z = float(line.split()[2])*scaling-0.5
		ifghost = True
		tag = 'BOT'
		random_pebbles.new_pebble(p_x,p_y,p_z,ifghost,tag)    # add new pebble
	bot_gpebble_file.close()
else:
	# use mirror approach to generate top and bot ghost pebbles
	number_of_pebbles = int(random_pebbles.number_of_pebbles())
	for ipebble in range(0,number_of_pebbles):
		p_x = random_pebbles.pebble_xyz(ipebble)[0]
		p_y = random_pebbles.pebble_xyz(ipebble)[1]
		p_z = random_pebbles.pebble_xyz(ipebble)[2]

		d_to_top = cyl_top - p_z
		if d_to_top <= pebble_diameter*ghost_pebble_tol_top_bot: 
			# construct ghost pebble
			g_x = p_x
			g_y = p_y
			g_z = cyl_top + d_to_top
			ifghost = True
			tag = 'TOP'
			random_pebbles.new_pebble(g_x,g_y,g_z,ifghost,tag)
		
		d_to_bot = p_z - cyl_bot
		if d_to_bot <= pebble_diameter*ghost_pebble_tol_top_bot: 
		# construct ghost pebble
			g_x = p_x
			g_y = p_y
			g_z = cyl_bot - d_to_bot
			ifghost = True
			tag = 'BOT'
			random_pebbles.new_pebble(g_x,g_y,g_z,ifghost,tag)

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


minAngle = 20.0

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

#
if_bend_polygon_corner = False

if if_bend_polygon_corner:
	# bend corner vertice to reduce non-right-hand elements
	beta = 0.2  # distance of move
	delta = 15.0 # angles for search
	for iface in range(0,vor_face.nf):
		vor_face.bend_corner_vertice(iface,vor_vert,random_pebbles,beta,delta)
	
	beta = 0.1  # distance of move
	delta = 10.0 # angles for search
	for iface in range(0,vor_face.nf):
		vor_face.bend_corner_vertice(iface,vor_vert,random_pebbles,beta,delta)

	
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
	if not if_top_bot_ghost_pebble_from_files:
		tag = 'BOT'
		vor_face.project_to_plane(iface,vor_vert,random_pebbles,cyl_bot,tag)
		tag = 'TOP'
		vor_face.project_to_plane(iface,vor_vert,random_pebbles,cyl_top,tag)
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
longEdge = 0.6*pebble_diameter
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

#=============================================================================
#============================================================================

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
	
	
vtkFileName = 'polygon_7.vtk'
print ' dumping '+ vtkFileName
dump_poly_vtk(vor_vert,vor_face,vtkFileName)

#=============================================================================
# new divide for octagon,heptagon,hexagon, pentagon
for iface in range(0,vor_face.nf):
	polygon_divide_for_octagon_new(iface,vor_face,vor_vert)
	polygon_divide_for_heptagon_new(iface,vor_face,vor_vert)
	polygon_divide_for_hexagon_new(iface,vor_face,vor_vert)
	#polygon_divide_for_pentagon_new(iface,vor_face,vor_vert)
	polygon_divide_for_pentagon_tri_split(iface,vor_face,vor_vert)
	polygon_divide_for_pentagon_tri_split2(iface,vor_face,vor_vert)
	polygon_divide_for_pentagon_bi_split(iface,vor_face,vor_vert)

for iface in range(0,vor_face.nf):
	vor_face.find_edge_info(iface)
	vor_face.calculate_edge_length(iface,vor_vert)
	vor_face.clean_up_zero_length_edge(iface,vor_vert)

for iface in range(0,vor_face.nf):
	vor_face.find_edge_info(iface)
	vor_face.calculate_edge_length(iface,vor_vert)
	vor_face.add_f_center(iface,vor_vert,random_pebbles,pebble_diameter)
	vor_face.calculate_vertice_angles(iface,vor_vert)
	
#============================================================================

vtkFileName = 'polygon_8.vtk'
print ' dumping '+ vtkFileName
dump_poly_vtk(vor_vert,vor_face,vtkFileName)
#============================================================================

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

vtkFileName = 'polygon_9.vtk'
print ' dumping '+ vtkFileName
dump_poly_vtk(vor_vert,vor_face,vtkFileName)

#================================================================================
#=================================================================================
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

endtime = time.time()
print 'time used:' + str(endtime-starttime)

#================================================================================
# make the valid polygons evenly distributed in processors
n_valid_face = 0 # number of valid face
valid_faces = [] # array of valid faces
eqtot = 0
for iface in range(0,vor_face.nf):
	ifghost = vor_face.ifghost[iface]
	ifcollapse = vor_face.ifcollapse[iface]
	ifchamfer = vor_face.ifchamfer[iface]
	if (not ifghost) and (not ifcollapse):
		n_valid_face = n_valid_face + 1
		valid_faces.append(iface)
		
		nedges = len(vor_face.f_to_v[iface])
		
		if ifchamfer:
			eqtot = eqtot + nedges*10
		else:
			if nedges ==5:
				eqtot = eqtot + 15
			else:
				eqtot = eqtot + nedges
		
print ' total valid polygons: '+ str(n_valid_face)
print ' total estimated quads: '+ str(eqtot)

#=================================================================================
# now we have the all polygons.
# paralel divie to quads, and the quads to hexs

# rearange polygons, so wordload is more even

nqrange = int(math.ceil(eqtot/nprocs))+2
n_valid_face = 0 # number of valid face
valid_faces = [] # array of valid faces

start = [0]
iquad = 0
iquad2 = nqrange
for iface in range(0,vor_face.nf):
	ifghost = vor_face.ifghost[iface]
	ifcollapse = vor_face.ifcollapse[iface]
	ifchamfer = vor_face.ifchamfer[iface]
	if (not ifghost) and (not ifcollapse):
		n_valid_face = n_valid_face + 1
		valid_faces.append(iface)
		
		
		nedges = len(vor_face.f_to_v[iface])
		
		if ifchamfer:
			iquad = iquad + nedges*10
		else:
			if nedges ==5:
				iquad = iquad + 15
			else:
				iquad = iquad + nedges

		if iquad > iquad2:
			start.append(len(valid_faces)-1)
			iquad2 = iquad2 + nqrange

#print 'str(len(start)); '+ str(len(start))
#print start

end = []

for ip in range(0,nprocs-1):
	end.append(start[ip+1])	
end.append(n_valid_face)


jobs = []
	
for ip in range(0,nprocs):
	if if_top_bot_ghost_pebble_from_files:
		process = multiprocessing.Process(target=polygons_to_hexs,args=(vor_face,vor_vert,random_pebbles,pebble_diameter,geo,start[ip],end[ip],valid_faces,ip))
	else:
		process = multiprocessing.Process(target=polygons_to_hexs2,args=(vor_face,vor_vert,random_pebbles,pebble_diameter,geo,start[ip],end[ip],valid_faces,ip))
	jobs.append(process)
	
for j in jobs:
	j.start()

for j in jobs:
	j.join()

endtime = time.time()
print 'time used:' + str(endtime-starttime)
