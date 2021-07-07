import numpy as np
import math

from data_class import *
from splitting_subroutines import *


def set_vertices(v1,vert,new):
	# set vertices to new location
	# including mergeed vertics
	xyz1 = vert.xyz[v1]
	vert.xyz[v1] = new
	
	for iv_mgd in vert.mgd_v[v1]:
		vert.xyz[iv_mgd] = new

	return

def find_nbr_vertices(vert,face):
	print 'find_nbr_vertices'
	for iv in range(0,vert.nv):
		# 1st layer
		for iface1 in vert.v_to_f[iv]:
			for ivf1 in face.f_to_v[iface1]:
				vnbr = ivf1
				if (vnbr not in vert.nbr_v[iv]) and (vnbr!=iv):
					vert.nbr_v[iv].append(vnbr)
					vert.nbr_v[vnbr].append(iv)
					
	for iv in range(0,vert.nv):
		nbrv = vert.nbr_v[iv]
		nbrv = remove_duplicates(nbrv)
		vert.nbr_v[iv] = nbrv
	
def merge_vertices(vert,tol):
	print 'merge_vertices'
	for iv in range(0,vert.nv):
		xyz = vert.xyz[iv]
		nbrlist = vert.nbr_v[iv]

		for inbr in nbrlist:
			if (inbr!=iv) :
				nbr_xyz = vert.xyz[inbr]
				d = distance(xyz,nbr_xyz)
				if d < tol:
					if inbr not in vert.mgd_v[iv]:
						mid = line_split(xyz,nbr_xyz,0.5)
					
						# inherit merged vertice list
						vert.mgd_v[iv].extend(vert.mgd_v[inbr])
						#vert.mgd_v[iv].append(inbr)
						vert.mgd_v[iv]= remove_duplicates(vert.mgd_v[iv])

						# copy merged vertice list to all merged vertices
						for mgd in vert.mgd_v[iv]:
							if (mgd!=iv) :
								vert.mgd_v[mgd] = vert.mgd_v[iv]
					
						set_vertices(iv,vert,mid)
					
						# inherit nbr vertice list
						vert.nbr_v[iv].extend(vert.nbr_v[inbr])
						#vert.nbr_v[iv].append(inbr)
						vert.nbr_v[iv]= remove_duplicates(vert.nbr_v[iv])
						#vert.nbr_v[inbr] = vert.nbr_v[iv]
						for mgd in vert.mgd_v[iv]:
							if (mgd!=iv) :
								vert.nbr_v[mgd] = vert.nbr_v[iv]
					
					
def merge_vertices_zplane(vert,tol,plz,dlz):
	# merge_verticece on z plane defined by plz+-dlz
	print 'merge_vertices_zplane'
	for iv in range(0,vert.nv):
		xyz = vert.xyz[iv]
		zz = xyz[2]
		if ((plz-dlz)<zz) and (zz<(plz+dlz)):
			nbrlist = vert.nbr_v[iv]
			
			for inbr in nbrlist:
				nbr_xyz = vert.xyz[inbr]
				d = distance(xyz,nbr_xyz)
				if d < tol:
					if inbr not in vert.mgd_v[iv]:
						mid = line_split(xyz,nbr_xyz,0.5)
						
						# inherit merged vertice list
						vert.mgd_v[iv].extend(vert.mgd_v[inbr])
						vert.mgd_v[iv].append(inbr)
						vert.mgd_v[iv]= remove_duplicates(vert.mgd_v[iv])

						# copy merged vertice list to all merged vertices
						for mgd in vert.mgd_v[iv]:
							vert.mgd_v[mgd] = vert.mgd_v[iv]
					
						set_vertices(iv,vert,mid)
					
						# inherit merged vertice list
						vert.nbr_v[iv].extend(vert.nbr_v[inbr])
						vert.nbr_v[iv].append(inbr)
						vert.nbr_v[iv]= remove_duplicates(vert.nbr_v[iv])
						vert.nbr_v[inbr] = vert.nbr_v[iv]
	
def merge_vertices_zplane(vert,tol,plz,dlz):
	# merge_verticece on z plane defined by plz+-dlz
	print 'merge_vertices_zplane'
	for iv in range(0,vert.nv):
		xyz = vert.xyz[iv]
		zz = xyz[2]
		if ((plz-dlz)<zz) and (zz<(plz+dlz)):
			nbrlist = vert.nbr_v[iv]
			
			for inbr in nbrlist:
				nbr_xyz = vert.xyz[inbr]
				d = distance(xyz,nbr_xyz)
				if d < tol:
					if inbr not in vert.mgd_v[iv]:
						mid = line_split(xyz,nbr_xyz,0.5)
						
						# inherit merged vertice list
						vert.mgd_v[iv].extend(vert.mgd_v[inbr])
						vert.mgd_v[iv].append(inbr)
						vert.mgd_v[iv]= remove_duplicates(vert.mgd_v[iv])

						# copy merged vertice list to all merged vertices
						for mgd in vert.mgd_v[iv]:
							vert.mgd_v[mgd] = vert.mgd_v[iv]
					
						set_vertices(iv,vert,mid)
					
						# inherit merged vertice list
						vert.nbr_v[iv].extend(vert.nbr_v[inbr])
						vert.nbr_v[iv].append(inbr)
						vert.nbr_v[iv]= remove_duplicates(vert.nbr_v[iv])
						vert.nbr_v[inbr] = vert.nbr_v[iv]	