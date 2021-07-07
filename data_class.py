import numpy as np
import math
# define classes for data

from splitting_subroutines import *
from merge_vertices import *

class pebbles:						# class for pebble list
	def __init__(self,name):
		self.name = name			# name of instance
		self.npebbles = 0			# total number of pebbles
		self.rpebbles = 0			# real pebble number
		self.gpebbles = 0			# ghost pebble number
		self.xyz = []				# pebbles coordinates
		self.p_to_f = []			# pebble_to_face pointer
		self.p_ghost = []			# list for real or ghost pebbles
		self.p_to_r = []			# pebble to region
		self.tag = []
		self.top_quads = []
		self.bot_quads = []
		self.voro_quads = []
	
	def new_pebble(self,x,y,z,ifghost,tag):			# add one new pebble
		self.xyz.append([x,y,z])
		self.p_to_f.append([])
		self.p_to_r.append(-1)
		self.npebbles = self.npebbles+1
		self.tag.append(tag) # tag: '' for real,'side' for side ghost pebble,'top' for top ghost pebble,'bot' for bot ghost pebble
		if ifghost:
			self.gpebbles = self.gpebbles +1
			self.p_ghost.append(ifghost)
		else:
			self.rpebbles = self.rpebbles +1
			self.p_ghost.append(ifghost)
		
		self.top_quads.append([])
		self.bot_quads.append([])
		self.voro_quads.append([])

	def pebble_xyz(self,ipebble):
		if ipebble > (self.npebbles-1):
			print 'ERROR: pebble ',str(ipebble),' has not been added yet !'
		return self.xyz[ipebble]
	
	def all_pebbles_xyz(self):
		return self.xyz

	def pebble_if_ghost(self,ipebble):
		if ipebble > (self.npebbles-1):
			print 'ERROR: pebble ',str(ipebble),' has not been added yet !'
		return self.p_ghost[ipebble]
	
	def number_of_pebbles(self):
		return self.npebbles
	
	def number_of_real_pebbles(self):
		return self.rpebbles
		
	def number_of_ghost_pebbles(self):
		return self.gpebbles
	
	def add_face_to_pebble(self,ipebble,iface):  # add face(ifae) to pebble(ipebble) 
		if ipebble > self.npebbles:
			print 'ERROR: pebble ',str(ipebble),' has not been added yet !'
		self.p_to_f[ipebble].append(iface)
	
	def add_region_to_pebble(self,ipebble,iregion):
		if ipebble > self.npebbles:
			print 'ERROR: pebble ',str(ipebble),' has not been added yet !'
		self.p_to_r[ipebble] = iregion
	
class voro_faces:
	def __init__(self,name):
		self.name = name
		self.nf = 0
		self.nf_real = 0
		self.nf_ghost = 0
		self.f_to_v = []
		self.f_to_r = []
		self.f_to_e = []
		self.f_to_p = []
		self.f_center = []
		self.f_pp_center = []
		self.f_normal = []
		self.quad = []
		self.ifghost = []
		self.ifcollapse = []
		self.ifchamfer = []
		self.angles = []
		self.lengths = []
		
		self.f_to_mid = [] 	# special list for chamfer polygon
							# mid edge point of each edge
	
	def new_face(self,vlist,plist):					# add new face
		self.f_to_v.append(vlist)
		self.f_to_e.append([])
		self.f_to_r.append([])
		self.f_to_p.append(plist)
		self.f_center.append([])
		self.f_pp_center.append([])
		self.f_normal.append([])
		self.angles.append([]) 
		self.lengths.append([])
		self.quad.append([])
		self.ifcollapse.append(False)
		self.ifchamfer.append(False)
		self.nf = self.nf +1
		self.ifghost.append(False)
		
		self.f_to_mid.append([])
		
		'''
		try:	
			vlist.index(-1)
		except ValueError: 
			self.ifghost.append(False) 		# not ghost face
			self.nf_real = self.nf_real + 1
		else:
			self.ifghost.append(True) # if find -1 in vlist, then this face is ghost face
			self.nf_ghost = self.nf_ghost + 1
		'''
	
	def add_pebble_to_face(self,iface,ipebble):				# add new pebble to face
		if iface > (self.nf-1):
			print 'ERROR: face ',str(iface),' has not been added yet !'
		self.f_to_p[iface].append(ipebble)
		if(len(f_to_p[iface])>2):
			print 'ERROR: more than 2 pebbles for a voro_face'
		# add warning, no more than 2 pebbles per face

	def face_to_point(self,iface):
		if iface > (self.nf-1):
			print 'ERROR: face ',str(iface),' has not been added yet !'
		return f_to_p[iface]

	def add_f_pp_center(self,iface,p1xyz,p2xyz,vert):
		# add pebble-pebble midpoint to this face.
		# however pebble-pebble midpoint does not necessary inside this face
		if iface > (self.nf-1):
			print 'ERROR: face ',str(iface),' has not been added yet !'

		fcx = (p1xyz[0]+p2xyz[0])/2.0
		fcy = (p1xyz[1]+p2xyz[1])/2.0
		fcz = (p1xyz[2]+p2xyz[2])/2.0
		self.f_pp_center[iface].append(fcx)
		self.f_pp_center[iface].append(fcy)
		self.f_pp_center[iface].append(fcz)
		fnx = (p1xyz[0]-p2xyz[0])
		fny = (p1xyz[1]-p2xyz[1])
		fnz = (p1xyz[2]-p2xyz[2])
		fn = [fnx,fny,fnz]
		fn_norm = np.linalg.norm(fn)
		self.f_normal[iface].append(fnx/fn_norm)
		self.f_normal[iface].append(fny/fn_norm)
		self.f_normal[iface].append(fnz/fn_norm)
	
	def redistribute(self,iface,vert,pebbles,diameter,r_chamfer):
		# redistribute face vertice if too close to pebbles
		if (not self.ifghost[iface]) and (not self.ifcollapse[iface]): 
			p1 = self.f_to_p[iface][0]
			p2 = self.f_to_p[iface][1]
			p1xyz = pebbles.xyz[p1]
			p2xyz = pebbles.xyz[p2]
			ppxyz = [p1xyz,p2xyz]
			
			p1tag = pebbles.tag[p1]
			p2tag =	pebbles.tag[p2]
			ptag = [p1tag,p2tag]
			
			if p1tag == '' and p2tag == '': # only deal with internal polygon
			
				# move vertice away from pebble 
				nv = len(self.f_to_v[iface])
				for iv in range(0,nv):
					v = self.f_to_v[iface][iv]
					vxyz = vert.xyz[v]
					new,moved = move_away_from_pebbles(vxyz,ppxyz,diameter,1.02,ptag) #1.02
					
					if moved:
						set_vertices(v,vert,new)
				
				#move vertice away from chamfer center as well
				nv = len(self.f_to_v[iface])
				for iv in range(0,nv):
					v = self.f_to_v[iface][iv]
					vxyz = vert.xyz[v]
					new,moved = move_away_from_chamfer(vxyz,ppxyz,diameter,r_chamfer)
					
					if moved:
						set_vertices(v,vert,new)
				
				
				# move vertice away if edge is too close to pebbles
				nedge = len(self.f_to_e[iface])
				for iedge in range(0,nedge):
					iv1 = self.f_to_e[iface][iedge][0]
					iv2 = self.f_to_e[iface][iedge][1]
					
					v1xyz = vert.xyz[iv1]
					v2xyz = vert.xyz[iv2]
					
					for ip in range(0,2):
						
						# d is distance between pebble center and edge
						# v is the normal vector from pebble center to edge
						dist,nvec,inside = edge_pebble_distance(v1xyz,v2xyz,ppxyz[ip])
						
						if dist <= (diameter/2.0)*1.02 and inside: #1.03,1.02
							shift = (diameter/2.0)*1.02 - dist 
							v1xyz[0] =  v1xyz[0] + shift*nvec[0]
							v1xyz[1] =  v1xyz[1] + shift*nvec[1]
							v1xyz[2] =  v1xyz[2] + shift*nvec[2]
							
							v2xyz[0] =  v2xyz[0] + shift*nvec[0]
							v2xyz[1] =  v2xyz[1] + shift*nvec[1]
							v2xyz[2] =  v2xyz[2] + shift*nvec[2]
							
							set_vertices(iv1,vert,v1xyz)
							set_vertices(iv2,vert,v2xyz)
				
				# move vertice away if edge is too close to chamfers
				#cxyz = [0,0,0]
				#for i in range(0,3):
				#	cxyz[i] = (ppxyz[0][i]+ppxyz[1][i])/2.0
				
				#nedge = len(self.f_to_e[iface])
				#for iedge in range(0,nedge):
				#	iv1 = self.f_to_e[iface][iedge][0]
				#	iv2 = self.f_to_e[iface][iedge][1]
				#	
				#	v1xyz = vert.xyz[iv1]
				#	v2xyz = vert.xyz[iv2]
				#	
				#	# d is distance between pebble center and edge
				#	# v is the normal vector from pebble center to edge
				#	dist,nvec,inside = edge_pebble_distance(v1xyz,v2xyz,cxyz)
				#		
				#	if dist <= (diameter/2.0)*r_chamfer*1.1 :
				#		shift = (diameter/2.0)*r_chamfer*1.1 - dist 
				#		v1xyz[0] =  v1xyz[0] + shift*nvec[0]
				#		v1xyz[1] =  v1xyz[1] + shift*nvec[1]
				#		v1xyz[2] =  v1xyz[2] + shift*nvec[2]
				#			
				#		v2xyz[0] =  v2xyz[0] + shift*nvec[0]
				#		v2xyz[1] =  v2xyz[1] + shift*nvec[1]
				#		v2xyz[2] =  v2xyz[2] + shift*nvec[2]
				#			
				#		set_vertices(iv1,vert,v1xyz)
				#		set_vertices(iv2,vert,v2xyz)
				

		return
		
	def add_f_center(self,iface,vert,pebbles,diameter):
		# add face center
		self.f_center[iface] = []
		if (not self.ifghost[iface]) and (not self.ifcollapse[iface]): 
			nedges = len(self.f_to_v[iface])
			totLength = 0.0
			length = []
			mid = []
			for i in range(0,nedges):
				i1 = i
				i2 = i+1
				if(i2>nedges-1):
					i2 = 0
			
				v1 = self.f_to_v[iface][i1]
				v2 = self.f_to_v[iface][i2]
			
				midv = line_split(vert.xyz[v1],vert.xyz[v2],0.5)
				mid.append(midv)
				d = distance(vert.xyz[v1],vert.xyz[v2])
				length.append(d)
				totLength = totLength + d
		
			fcx = 0.0
			fcy = 0.0
			fcz = 0.0
		
			for i in range(0,nedges):
				#print length[i]
				fcx = fcx + mid[i][0]*length[i]/totLength
				fcy = fcy + mid[i][1]*length[i]/totLength
				fcz = fcz + mid[i][2]*length[i]/totLength
		
			self.f_center[iface] = [fcx,fcy,fcz]

			center = self.f_center[iface]
			
			# however, if f_center is inside the one of the two pebbles
			p1 = self.f_to_p[iface][0]
			p2 = self.f_to_p[iface][1]
			p1xyz = pebbles.xyz[p1]
			p2xyz = pebbles.xyz[p2]
			pxyz = [p1xyz,p2xyz]
			
			p1tag = pebbles.tag[p1]
			p2tag =	pebbles.tag[p2]
			ptag = [p1tag,p2tag]
			if p1tag == '' and p2tag == '': # only deal with internal polygon
				center,moved = move_away_from_pebbles(center,pxyz,diameter,1.02,ptag)
				self.f_center[iface] = center
		return
		
	def find_edge_info(self,iface):

		self.f_to_e[iface] = []
		if not self.ifghost[iface] and not self.ifcollapse[iface]: # only do this for non-ghost face
				
			nverts_in_face = len(self.f_to_v[iface])	
			iv1 = 0
			iv2 = 1
			self.f_to_e[iface] = []
			# for all edges except last one
			for iedge in range(0,nverts_in_face-1):
			
				v1 = self.f_to_v[iface][iv1]
 				v2 = self.f_to_v[iface][iv2]
 														
				self.f_to_e[iface].append([v1,v2])					
				iv1 = iv1+1
				iv2 = iv2+1
				
			# for last edge
			iv1 = nverts_in_face-1
			iv2 = 0
			v1 = self.f_to_v[iface][iv1]
 			v2 = self.f_to_v[iface][iv2]
 			self.f_to_e[iface].append([v1,v2])
	
	def clean_up_zero_length_edge(self,iface,vert):
		if not self.ifghost[iface] and not self.ifcollapse[iface]:
			nedge = len(self.f_to_v[iface])
			nedge_reduced = nedge
			
			new_ftv = []
			for iedge in range(0,nedge):

				v1 = self.f_to_e[iface][iedge][0]
				v2 = self.f_to_e[iface][iedge][1]
					
				if self.lengths[iface][iedge]< 1e-5:
					# remove this edge
					nedge_reduced = nedge_reduced -1
				else:
					# construct new
					new_ftv.append(v1)
				
				if nedge_reduced < 3: self.ifcollapse[iface] = True
			
			self.f_to_v[iface] = new_ftv

		return

	def fix_concave_angle(self,iface,vert):
		nvert = len(self.f_to_v[iface])
		if (not self.ifghost[iface]) and (not self.ifcollapse[iface]):
			for ivert in range(0,nvert):
				angle = self.angles[iface][ivert]
				if angle > 180:
					v1  = ivert
					v2 = v1+1
					if v2 > (nvert-1): v2 = 0
					v3 = v1-1
					if v3 < 0: v3 = nvert-1
	
					v1 = self.f_to_v[iface][v1]
					v2 = self.f_to_v[iface][v2]
					v3 = self.f_to_v[iface][v3]
				
					xyz2 = vert.xyz[v2]
					xyz3 = vert.xyz[v3]
				
					mid = line_split(xyz2,xyz3,0.5)
			
					# move v1 to mid
					set_vertices(v1,vert,mid)
		return 
	
	def add_mid_vertice(self,iface,vert,longEdge):
		# add mid edge vertice of edge if it is too long
 		#
		# nvert = len(self.f_to_v[iface])
		nedge = len(self.f_to_v[iface])
		if (not self.ifghost[iface]) and (not self.ifcollapse[iface]):
			
			new_ftv = []
			for iedge in range(0,nedge):
				edgeLength = self.lengths[iface][iedge]

				v1 = self.f_to_e[iface][iedge][0]
				v2 = self.f_to_e[iface][iedge][1]
				
				if edgeLength > longEdge:
					
					xyz1 = vert.xyz[v1]
					xyz2 = vert.xyz[v2]
				
					mid = line_split(xyz1,xyz2,0.5)
			
					# create new vertice
					vert.new_vertice(mid[0],mid[1],mid[2])
					
					newvert = vert.nv - 1
					# adjust f_to_v list
					new_ftv.append(v1)
					new_ftv.append(newvert)
					new_ftv.append(v2)
				else:
					new_ftv.append(v1)
					new_ftv.append(v2)
					
			new_ftv = remove_duplicates(new_ftv)
			self.f_to_v[iface] = new_ftv
		return
		
	def fix_sharp_triangle(self,iface,vert,minAngle):
		# if single angle is very small.
		# collapse the other two vertices
		#
		# if two angles are both small, expend the other vertices
		nedge = len(self.f_to_v[iface])
		if (not self.ifghost[iface]) and (not self.ifcollapse[iface]) and (not self.ifchamfer[iface]) and nedge ==3:
				
			if min(self.angles[iface]) < minAngle:
				# if there is only one angle is very small
				# 
				v1  = self.angles[iface].index(min(self.angles[iface]))
				v2 = v1+1
				if v2 > 2: v2 = 0
				v3 = v1-1
				if v3 < 0: v3 = 2
				
				# merge v2,v3
				v1 = self.f_to_v[iface][v1]
				v2 = self.f_to_v[iface][v2]
				v3 = self.f_to_v[iface][v3]
				
				xyz2 = vert.xyz[v2]
				xyz3 = vert.xyz[v3]
				
				mid = line_split(xyz2,xyz3,0.5)
					
				vert.mgd_v[v2].extend(vert.mgd_v[v3])
				vert.mgd_v[v2].append(v3)
				vert.mgd_v[v2]= remove_duplicates(vert.mgd_v[v2])

					
				for mgd in vert.mgd_v[v2]:
					vert.mgd_v[mgd] = vert.mgd_v[v2]
					
				set_vertices(v2,vert,mid)
				#set_vertices(v3,vert,mid)
				
				self.ifcollapse[iface] = True
		return
		
	def remove_tiny_face(self,iface,vert,minArea):
		# if single angle is very small.
		# collapse the other two vertices
		#
		# if two angles are both small, expend the other vertices
		nedge = len(self.f_to_v[iface])
		if (not self.ifghost[iface]) and (not self.ifcollapse[iface]) and (not self.ifchamfer[iface]):
				
			# 1st step, get area of face, assume non-concate shape
			area = 0.0
			for iv1 in range(0,nedge):
				iv2= iv1+1
				if iv2 > (nedge -1): iv2 = 0
				
				v1 = self.f_to_v[iface][iv1]
				v2 = self.f_to_v[iface][iv2]
				
				xyz1 = vert.xyz[v1]
				xyz2 = vert.xyz[v2]
				
				sub_area = triangle_area(self.f_center[iface],xyz1,xyz2)
				area = area + sub_area
				
			if area  < minArea:
				iv1 = 0
				v1 = self.f_to_v[iface][iv1]
				for iv2 in range(1,nedge):

					v2 = self.f_to_v[iface][iv2]
				
					vert.mgd_v[v1].extend(vert.mgd_v[v2])
					vert.mgd_v[v1].append(v2)
					vert.mgd_v[v1]= remove_duplicates(vert.mgd_v[v1])

					
				for mgd in vert.mgd_v[v1]:
					vert.mgd_v[mgd] = vert.mgd_v[v1]
					
				set_vertices(v1,vert,self.f_center[iface])
				
				self.ifcollapse[iface] = True
		return

	
	def adjust_tri_longest_edge_for_top_bot(self,iface,vert,cyl_bot,cyl_top,pebble_diameter):
		# adjust the longest edge for triangle at top and bot 
		
		nedge = len(self.f_to_v[iface])
		if (not self.ifghost[iface]) and (not self.ifcollapse[iface]) and nedge ==3:
			if (self.f_center[iface][2] > cyl_top-0.5*pebble_diameter) or (self.f_center[iface][2] < cyl_bot+0.5*pebble_diameter):
				nedge = len(self.f_to_v[iface])
				avg_length = 0.0
				for iedge in range(0,nedge):			
					avg_length = avg_length + self.lengths[iface][iedge]/nedge
			
				shrinkEdge = -1
				longestEdge = max(self.lengths[iface])
				if longestEdge > 1.25*avg_length:
					#
					shrinkEdge = self.lengths[iface].index(longestEdge)
				
				if shrinkEdge >= 0:
					# if this edge is too long, shrink it to its 80%
					v1 = self.f_to_e[iface][shrinkEdge][0]
					v2 = self.f_to_e[iface][shrinkEdge][1]
					close_two_vertices(v1,v2,vert,0.85)
				
	
	def calculate_edge_length(self,iface,vert):
		self.lengths[iface] = []
		if (not self.ifghost[iface]) and (not self.ifcollapse[iface]):
			nedge = len(self.f_to_v[iface])
			for iedge in range(0,nedge):			
				v1 = self.f_to_e[iface][iedge][0]
				v2 = self.f_to_e[iface][iedge][1]
				d = distance(vert.xyz[v1],vert.xyz[v2])
				self.lengths[iface].append(d)
	
	def calculate_vertice_angles(self,iface,vert):
		#calculate angle of each vertices
		self.angles[iface] = []
		nvert = len(self.f_to_v[iface])
		if (not self.ifghost[iface]) and (not self.ifcollapse[iface]):
			for ivert in range(0,nvert):
				iv1  = ivert
				iv2 = iv1+1
				if iv2 > nvert-1: iv2 = 0
				iv3 = iv1-1
				if iv3 < 0: iv3 = nvert-1
				
				v1 = self.f_to_v[iface][iv1]
				v2 = self.f_to_v[iface][iv2]
				v3 = self.f_to_v[iface][iv3]
				
				v1 = self.f_to_v[iface][iv1]
				v2 = self.f_to_v[iface][iv2]
				v3 = self.f_to_v[iface][iv3]
				
				v1xyz = vert.xyz[v1]
				v2xyz = vert.xyz[v2]
				v3xyz = vert.xyz[v3]
								
				fc = self.f_center[iface]
				
				angle1,nv = angles_three_vertices(v1xyz,v2xyz,fc)
				angle2,nv = angles_three_vertices(v1xyz,v3xyz,fc)
				
				angle = angle1 + angle2
				
				self.angles[iface].append(angle)
								
		return
		
	def project_to_plane(self,iface,vert,pebbles,planez,tag):
		nvert = len(self.f_to_v[iface])
		if (not self.ifghost[iface]) and (not self.ifcollapse[iface]):
			p1 = self.f_to_p[iface][0]
			p2 = self.f_to_p[iface][1]
			p1tag = pebbles.tag[p1]
			p2tag =	pebbles.tag[p2]
			if p1tag == tag or p2tag == tag: 
				for iv in range(0,nvert):
					v =  self.f_to_v[iface][iv]
					xyz = vert.xyz[v]
					new = [xyz[0],xyz[1],planez]
					set_vertices(v,vert,new)
					
	def project_to_sidewall(self,iface,vert,pebbles,cyl_radius):
		nvert = len(self.f_to_v[iface])
		if (not self.ifghost[iface]) and (not self.ifcollapse[iface]):
			p1 = self.f_to_p[iface][0]
			p2 = self.f_to_p[iface][1]
			p1tag = pebbles.tag[p1]
			p2tag =	pebbles.tag[p2]
			if p1tag == 'SW ' or p2tag == 'SW ': 
				for iv in range(0,nvert):
					v =  self.f_to_v[iface][iv]
					xyz = vert.xyz[v]
					
					xx = xyz[0]
					yy = xyz[1]
					zz = xyz[2]
					rr = (xx**2.0+yy**2.0)**0.5
					
					new_xx = cyl_radius*xx/rr
					new_yy = cyl_radius*yy/rr
					
					new = [new_xx,new_yy,zz]
					set_vertices(v,vert,new)
					
	def detect_chamfer(self,iface,p1xyz,p2xyz,diameter,alpha):
		if (not self.ifghost[iface]) and (not self.ifcollapse[iface]):
			pdist = distance(p1xyz,p2xyz)
			if pdist < (alpha*diameter):
				# if pebble distance is too close. 
				# flag this polygon with chamfer
				self.ifchamfer[iface] = True
				
	def bend_corner_vertice(self,iface,vert,pebbles,beta,delta):
		# bend cornver vetice, so it project to pebble, it wont no the a line with with edges of polygon
		#delta = 15.0 # degree
		#beta = 0.2  # ratio
		
		p1 = self.f_to_p[iface][0]
		p2 = self.f_to_p[iface][1]
		p1tag = pebbles.tag[p1]
		p2tag =	pebbles.tag[p2]
		
		p1xyz = pebbles.xyz[p1]
		p2xyz = pebbles.xyz[p2]
		pxyz = [p1xyz,p2xyz]
		
		nvert = len(self.f_to_v[iface])
		if (not self.ifghost[iface]) and (not self.ifcollapse[iface]) and (p1tag != 'SW ') and (p2tag != 'SW '):
			for ivert in range(0,nvert):
				iv1  = ivert   # this vertice
				iv2 = iv1+1    # nearby vertice
				if iv2 > nvert-1: iv2 = 0
				iv3 = iv1-1    # nearby vertice
				if iv3 < 0: iv3 = nvert-1
				
				v1 = self.f_to_v[iface][iv1]
				v2 = self.f_to_v[iface][iv2]
				v3 = self.f_to_v[iface][iv3]
				
				v1xyz = vert.xyz[v1]
				v2xyz = vert.xyz[v2]
				v3xyz = vert.xyz[v3]
				
				d12 = distance(v1xyz,v2xyz)
				d13 = distance(v1xyz,v3xyz)
				dmin = min([d12,d13])
				
				if d12 < 1e-5:
					iv2 = iv2 + 1
				if d13 < 1e-5:
					iv3 = iv3 - 1
					
				v1 = self.f_to_v[iface][iv1]
				v2 = self.f_to_v[iface][iv2]
				v3 = self.f_to_v[iface][iv3]
				
				v1xyz = vert.xyz[v1]
				v2xyz = vert.xyz[v2]
				v3xyz = vert.xyz[v3]
				
				d12 = distance(v1xyz,v2xyz)
				d13 = distance(v1xyz,v3xyz)
				dmin = min([d12,d13])
					
					
				vec12 = vector_minus(v2xyz,v1xyz)
				vec13 = vector_minus(v3xyz,v1xyz)
				
				nv23 = np.cross(vec12,vec13)
				#nv23_norm = np.linalg.norm(nv23)
				nv23_norm = (nv23[0]**2.0+nv23[1]**2.0+nv23[2]**2.0)**0.5
				nv = [nv23[0]/nv23_norm,nv23[1]/nv23_norm,nv23[2]/nv23_norm]
				
				
				v1p1 = vector_minus(p1xyz,v1xyz)
				v1p1_norm = np.linalg.norm(v1p1)
				vp1 = [v1p1[0]/v1p1_norm,v1p1[1]/v1p1_norm,v1p1[2]/v1p1_norm]
				
				v1p2 = vector_minus(p2xyz,v1xyz)
				v1p2_norm = np.linalg.norm(v1p2)
				vp2 = [v1p2[0]/v1p2_norm,v1p2[1]/v1p2_norm,v1p2[2]/v1p2_norm]
				
				
				ap1 = angles_two_vectors(vp1,nv)
				ap2 = angles_two_vectors(vp2,nv)
				
				vangle = self.angles[iface][iv1]
				vangle_max = 120.0
				
				
				if vangle < vangle_max and ap1 < ap2 and ap2 <= 90.0+delta:  # bend negative nv
					v1new = [0,0,0]
					v1new[0] = v1xyz[0] - nv[0]*dmin*beta
					v1new[1] = v1xyz[1] - nv[1]*dmin*beta
					v1new[2] = v1xyz[2] - nv[2]*dmin*beta
					set_vertices(v1,vert,v1new)
				
				if vangle < vangle_max and ap1 < ap2 and ap1 >= 90.0-delta:  # bend postive nv
					v1new = [0,0,0]
					v1new[0] = v1xyz[0] + nv[0]*dmin*beta
					v1new[1] = v1xyz[1] + nv[1]*dmin*beta
					v1new[2] = v1xyz[2] + nv[2]*dmin*beta
					set_vertices(v1,vert,v1new)
				
				if vangle < vangle_max and ap2 < ap1 and ap1 <= 90.0+delta:  # bend negative nv
					v1new = [0,0,0]
					v1new[0] = v1xyz[0] - nv[0]*dmin*beta
					v1new[1] = v1xyz[1] - nv[1]*dmin*beta
					v1new[2] = v1xyz[2] - nv[2]*dmin*beta
					set_vertices(v1,vert,v1new)
				
				if vangle < vangle_max and ap2 < ap1 and ap2 >= 90.0-delta:  # bend postive nv
					v1new = [0,0,0]
					v1new[0] = v1xyz[0] + nv[0]*dmin*beta
					v1new[1] = v1xyz[1] + nv[1]*dmin*beta
					v1new[2] = v1xyz[2] + nv[2]*dmin*beta
					set_vertices(v1,vert,v1new)
				
class voro_regions():
	def __init__(self,name):
		self.name = name
		self.nregions = 0
		self.r_to_v = []
		self.r_to_p = []

		self.r_to_v_for_non_ghost = []
		self.nregion_non_ghost = 0

	def new_region(self,vertices_list):
		self.r_to_v.append(vertices_list)
		self.r_to_p.append(-1)
		self.nregions = self.nregions +1

	def add_vert_to_region(self,iregion,ivert):
		if iregion > (self.nregions-1):
			print 'ERROR: region ',str(iregion),' has not been added yet !'
		return self.r_to_v[iregion].append(ivert)
	
	def add_pebble_to_region(self,iregion,ipebble):
		if iregion > (self.nregions-1):
			print 'ERROR: region ',str(iregion),' has not been added yet !'
		self.r_to_p[iregion] = ipebble

	def vertices_in_region(self,iregion):
		if iregion > (self.nregions-1):
			print 'ERROR: region ',str(iregion),' has not been added yet !'
		return self.r_to_v[iregion]
		
	def all_regions(self):
		return self.r_to_v
	
	def all_non_ghost_regions(self):
		return self.r_to_v_for_non_ghost
	
	def get_non_ghost_regions(self):
		for iregion in range(0,self.nregions):
			ifghost = False
			for ivert in range(0,len(self.r_to_v[iregion])):
				if self.r_to_v[iregion][ivert] < 0:
					ifghost = True
			
			if not ifghost:
				self.r_to_v_for_non_ghost.append(self.r_to_v[iregion])
				self.nregion_non_ghost = self.nregion_non_ghost + 1
		

class voro_edges():
	def __init__(self,name):
		self.name = name
		self.nedges = 0
		self.e_to_v = []

	def new_edge(self,ivert1,ivert2):
		self.e_to_v.append([ivert1,ivert2])
		self.nedges = self.nedges +1
	
	def edge_vertices(self,iedge):
		if iedge > (self.nedges-1):
			print 'ERROR: edge ',str(iedge),' has not been added yet !'
		return self.e_to_v[iedge]

	def all_edges(self):
		return self.e_to_v

class voro_vertices():
	def __init__(self,name):
		self.name = name
		self.nv = 0
		self.xyz = []
		self.e_v_list = [] # equivalent vertice list
		self.merged = []
		self.v_to_f = []
		self.nbr_v = []
		self.mgd_v = []
		
	def new_vertice(self,x,y,z):
		self.xyz.append([x,y,z])
		self.nv = self.nv +1
		self.e_v_list.append([])
		self.merged.append(False)
		self.v_to_f.append([])
		self.nbr_v.append([self.nv-1])
		self.mgd_v.append([self.nv-1])
	
	def vertice_xyz(self,ivert):
		if ivert > (self.nv-1):
			print 'ERROR: vertice ',str(ivert),' has not been added yet !'
		return self.xyz[ivert]

	def all_vertices_xyz(self):
		return self.xyz
		
	def adjust_vertice(self,ivert,x,y,z):
		if ivert > (self.nv-1):
			print 'ERROR: vertice ',str(ivert),' has not been added yet !'
		self.xyz[ivert][0] = x
		self.xyz[ivert][1] = y
		self.xyz[ivert][2] = z
	

class voro_quads():
	def __init__(self,name):
		self.name = name
		self.nquads = 0
		self.quad_to_p = []
		self.xyz = []
		self.mid_xyz = []
		self.tag=[] # tag for boundary condition
		self.edge_tag = [] # tag for each edge
		self.linearFlag = []
		
		# the rest is for quads smoothing.
		self.connect_quads_and_vertices = []
		self.newxyz = []
		self.quad_to_quad = []
		
		
	def new_quad(self,qxyz,p1,p2):
		#qxyz =[4][3]
		self.xyz.append(qxyz)
		self.quad_to_p.append([p1,p2])
		self.nquads = self.nquads +1
		self.mid_xyz.append([])
		self.tag.append('')
		self.edge_tag.append(['E  ','E  ','E  ','E  '])
		self.linearFlag.append(False)
		
		self.connect_quads_and_vertices.append([[],[],[],[]])
		self.newxyz.append(qxyz)
		self.quad_to_quad.append([])
		
	def set_tag(self,iquad,tag):
		self.tag[iquad] = tag
	
	def add_mid(self,iquad,equad):
		for i in range(0,4):
			self.mid_xyz[iquad].append(equad[i])
	
	def add_linear_mid(self,iquad):
		for i in range(0,4):
			i1 = i
			i2 = i+1
			if i2 > 3:
				i2 = 0
			midx = (self.xyz[iquad][i1][0] + self.xyz[iquad][i2][0])/2.0
			midy = (self.xyz[iquad][i1][1] + self.xyz[iquad][i2][1])/2.0
			midz = (self.xyz[iquad][i1][2] + self.xyz[iquad][i2][2])/2.0
			mid = [midx,midy,midz]
			
			#midnew,moved = move_away_from_pebbles_global(mid,all_pebbles,pebble_diameter,1.01,ptag)
			
			self.mid_xyz[iquad].append(mid)
			
	def linear_mid(self,iquad):
		for i in range(0,4):
			i1 = i
			i2 = i+1
			if i2 > 3:
				i2 = 0
			midx = (self.xyz[iquad][i1][0] + self.xyz[iquad][i2][0])/2.0
			midy = (self.xyz[iquad][i1][1] + self.xyz[iquad][i2][1])/2.0
			midz = (self.xyz[iquad][i1][2] + self.xyz[iquad][i2][2])/2.0
			mid = [midx,midy,midz]
			
			#midnew,moved = move_away_from_pebbles_global(mid,all_pebbles,pebble_diameter,1.01,ptag)
			
			self.mid_xyz[iquad][i] = mid
	
	def add_curve_mid(self,iquad,pxyz,pdiameter):
		for i in range(0,4):
			i1 = i
			i2 = i+1
			if i2 > 3:
				i2 = 0

			curveMid = line_mid_project_to_sphere(self.xyz[iquad][i1],self.xyz[iquad][i2],pxyz,0.5)	
			self.mid_xyz[iquad].append(curveMid)

class voro_hex8s():
	def __init__(self,name):
		self.name = name
		self.nhex8s= 0
		self.quad_bot = []  # far-pebble side quad
		self.quad_top = []  # pebble side quad
		self.hex_to_p = []

	def new_hex8(self,q1,q2,p):

		self.quad_bot.append(q1)
		self.quad_top.append(q2)
		self.hex_to_p.append(p)
		self.nhex8s = self.nhex8s +1

class nek_hex20s():
	def __init__(self,name):
		self.name = name
		self.nhex = 0
		self.v8 = []   		# 8 vertices
		self.e12 = []		# 12 edges 2nd order xyz
		self.e12type = []	# 2nd order type, mid-point, sphere
		self.s6 =[]			# side-set for each side
		self.nnrh = []
		self.lar = []

	def new_hex(self,v8,e12,s6):
		#v8 = [8][3]
		self.v8.append(v8)
		self.e12.append(e12)
		self.e12type.append([])
		self.s6.append(s6)
		self.nnrh.append(False)
		self.lar.append(False)
		self.nhex = self.nhex + 1
		
	def scaling(self,alpha):
		for ih in range(0,self.nhex):
			for i in range(0,8):
				for j in range(0,3):
					new = self.v8[ih][i][j] *alpha
					self.v8[ih][i][j] = new
		
			for i in range(0,12):
				for j in range(0,3):
					new =  self.e12[ih][i][j] *alpha
					self.e12[ih][i][j] = new
	def check_non_right_hand_elements_use_nek_method(self):
		print 'check_non_right_hand_elements_use_nek_method'
		nrh = 0
		for ih in range(0,self.nhex):
			self.nnrh[ih] = False
			# nek code for non-right-hand check
			#VOLUM0(A,B,C,O) = (A-O)X(B-O).(C-O)
			#V1= VOLUM0(XYZ(1,2,IE),XYZ(1,3,IE),XYZ(1,5,IE),XYZ(1,1,IE))
			#V2= VOLUM0(XYZ(1,4,IE),XYZ(1,1,IE),XYZ(1,6,IE),XYZ(1,2,IE))
			#V3= VOLUM0(XYZ(1,1,IE),XYZ(1,4,IE),XYZ(1,7,IE),XYZ(1,3,IE))
			#V4= VOLUM0(XYZ(1,3,IE),XYZ(1,2,IE),XYZ(1,8,IE),XYZ(1,4,IE))
			#V5=-VOLUM0(XYZ(1,6,IE),XYZ(1,7,IE),XYZ(1,1,IE),XYZ(1,5,IE))
			#V6=-VOLUM0(XYZ(1,8,IE),XYZ(1,5,IE),XYZ(1,2,IE),XYZ(1,6,IE))
			#V7=-VOLUM0(XYZ(1,5,IE),XYZ(1,8,IE),XYZ(1,3,IE),XYZ(1,7,IE))
			#V8=-VOLUM0(XYZ(1,7,IE),XYZ(1,6,IE),XYZ(1,4,IE),XYZ(1,8,IE))
			#
			#IF (V1.LE.0.0.OR.V2.LE.0.0.OR.V3.LE.0.0.OR.V4.LE.0.0.OR.V5.LE.0.0.OR.V6.LE.0.0.OR.V7.LE.0.0.OR.V8.LE.0.0    ) THEN
			
			# swap here
			v8here = []
			v8here.append(self.v8[ih][0])
			v8here.append(self.v8[ih][1])
			v8here.append(self.v8[ih][3])
			v8here.append(self.v8[ih][2])
			v8here.append(self.v8[ih][4])
			v8here.append(self.v8[ih][5])
			v8here.append(self.v8[ih][7])
			v8here.append(self.v8[ih][6])
			
			a1 = VOLUM0(v8here[1],v8here[2],v8here[4],v8here[0])
			a2 = VOLUM0(v8here[3],v8here[0],v8here[5],v8here[1])
			a3 = VOLUM0(v8here[0],v8here[3],v8here[6],v8here[2])
			a4 = VOLUM0(v8here[2],v8here[1],v8here[7],v8here[3])
			a5 = -VOLUM0(v8here[5],v8here[6],v8here[0],v8here[4])
			a6 = -VOLUM0(v8here[7],v8here[4],v8here[1],v8here[5])
			a7 = -VOLUM0(v8here[4],v8here[7],v8here[2],v8here[6])
			a8 = -VOLUM0(v8here[6],v8here[5],v8here[3],v8here[7])
			if (a1<=0.0) or (a2<=0.0) or (a3<=0.0) or (a4<=0.0) or (a5<=0.0) or (a6<=0.0) or (a7<=0.0) or (a8<=0.0):
				#print a1,a2,a3,a4,a5,a6,a6,a8
				nrh = nrh + 1
				self.nnrh[ih] = True
				# nnrh: nek non right hand
			
		print '============================================'
		print str(nrh)+' non right hand elements detected using nek method'
		print '============================================'
		
	def check_non_right_hand_elements_use_nek_method_range(self,start,end,ftag,ip):
		print 'check_non_right_hand_elements_use_nek_method'
		nrh = 0
		nrhh = []
		for ih in range(start,end):
			# nek code for non-right-hand check
			#VOLUM0(A,B,C,O) = (A-O)X(B-O).(C-O)
			#V1= VOLUM0(XYZ(1,2,IE),XYZ(1,3,IE),XYZ(1,5,IE),XYZ(1,1,IE))
			#V2= VOLUM0(XYZ(1,4,IE),XYZ(1,1,IE),XYZ(1,6,IE),XYZ(1,2,IE))
			#V3= VOLUM0(XYZ(1,1,IE),XYZ(1,4,IE),XYZ(1,7,IE),XYZ(1,3,IE))
			#V4= VOLUM0(XYZ(1,3,IE),XYZ(1,2,IE),XYZ(1,8,IE),XYZ(1,4,IE))
			#V5=-VOLUM0(XYZ(1,6,IE),XYZ(1,7,IE),XYZ(1,1,IE),XYZ(1,5,IE))
			#V6=-VOLUM0(XYZ(1,8,IE),XYZ(1,5,IE),XYZ(1,2,IE),XYZ(1,6,IE))
			#V7=-VOLUM0(XYZ(1,5,IE),XYZ(1,8,IE),XYZ(1,3,IE),XYZ(1,7,IE))
			#V8=-VOLUM0(XYZ(1,7,IE),XYZ(1,6,IE),XYZ(1,4,IE),XYZ(1,8,IE))
			#
			#IF (V1.LE.0.0.OR.V2.LE.0.0.OR.V3.LE.0.0.OR.V4.LE.0.0.OR.V5.LE.0.0.OR.V6.LE.0.0.OR.V7.LE.0.0.OR.V8.LE.0.0    ) THEN
			
			# swap here
			v8here = []
			v8here.append(self.v8[ih][0])
			v8here.append(self.v8[ih][1])
			v8here.append(self.v8[ih][3])
			v8here.append(self.v8[ih][2])
			v8here.append(self.v8[ih][4])
			v8here.append(self.v8[ih][5])
			v8here.append(self.v8[ih][7])
			v8here.append(self.v8[ih][6])
			
			a1 = VOLUM0(v8here[1],v8here[2],v8here[4],v8here[0])
			a2 = VOLUM0(v8here[3],v8here[0],v8here[5],v8here[1])
			a3 = VOLUM0(v8here[0],v8here[3],v8here[6],v8here[2])
			a4 = VOLUM0(v8here[2],v8here[1],v8here[7],v8here[3])
			a5 = -VOLUM0(v8here[5],v8here[6],v8here[0],v8here[4])
			a6 = -VOLUM0(v8here[7],v8here[4],v8here[1],v8here[5])
			a7 = -VOLUM0(v8here[4],v8here[7],v8here[2],v8here[6])
			a8 = -VOLUM0(v8here[6],v8here[5],v8here[3],v8here[7])
			if (a1<=0.0) or (a2<=0.0) or (a3<=0.0) or (a4<=0.0) or (a5<=0.0) or (a6<=0.0) or (a7<=0.0) or (a8<=0.0):
				#print a1,a2,a3,a4,a5,a6,a6,a8
				nrh = nrh + 1
				self.nnrh[ih] = True
				# nnrh: nek non right hand
				nrhh.append(ih)
			
		print '============================================'
		print str(nrh)+' non right hand elements detected using nek method'
		print '============================================'
		
		filename = 'dummy_'+ftag+'/dummy_'+str(ip)
		print 'writing files to '+filename
		dummyfile = open(filename,'w')
		nlines = len(nrhh)
		line = str(nlines)+'\n'
		dummyfile.write(line)
		for iline in range(0,nlines):
			line = str(nrhh[iline]) +'\n'
			dummyfile.write(line)
		dummyfile.close()
		print 'done writing files to '+filename
		
	def check_non_right_hand_elements_use_nek_method_one_element(self,ih):
		# nek code for non-right-hand check
		#VOLUM0(A,B,C,O) = (A-O)X(B-O).(C-O)
		#V1= VOLUM0(XYZ(1,2,IE),XYZ(1,3,IE),XYZ(1,5,IE),XYZ(1,1,IE))
		#V2= VOLUM0(XYZ(1,4,IE),XYZ(1,1,IE),XYZ(1,6,IE),XYZ(1,2,IE))
		#V3= VOLUM0(XYZ(1,1,IE),XYZ(1,4,IE),XYZ(1,7,IE),XYZ(1,3,IE))
		#V4= VOLUM0(XYZ(1,3,IE),XYZ(1,2,IE),XYZ(1,8,IE),XYZ(1,4,IE))
		#V5=-VOLUM0(XYZ(1,6,IE),XYZ(1,7,IE),XYZ(1,1,IE),XYZ(1,5,IE))
		#V6=-VOLUM0(XYZ(1,8,IE),XYZ(1,5,IE),XYZ(1,2,IE),XYZ(1,6,IE))
		#V7=-VOLUM0(XYZ(1,5,IE),XYZ(1,8,IE),XYZ(1,3,IE),XYZ(1,7,IE))
		#V8=-VOLUM0(XYZ(1,7,IE),XYZ(1,6,IE),XYZ(1,4,IE),XYZ(1,8,IE))
		#
		#IF (V1.LE.0.0.OR.V2.LE.0.0.OR.V3.LE.0.0.OR.V4.LE.0.0.OR.V5.LE.0.0.OR.V6.LE.0.0.OR.V7.LE.0.0.OR.V8.LE.0.0    ) THEN
		
		# swap here
		v8here = []
		v8here.append(self.v8[ih][0])
		v8here.append(self.v8[ih][1])
		v8here.append(self.v8[ih][3])
		v8here.append(self.v8[ih][2])
		v8here.append(self.v8[ih][4])
		v8here.append(self.v8[ih][5])
		v8here.append(self.v8[ih][7])
		v8here.append(self.v8[ih][6])
		
		a1 = VOLUM0(v8here[1],v8here[2],v8here[4],v8here[0])
		a2 = VOLUM0(v8here[3],v8here[0],v8here[5],v8here[1])
		a3 = VOLUM0(v8here[0],v8here[3],v8here[6],v8here[2])
		a4 = VOLUM0(v8here[2],v8here[1],v8here[7],v8here[3])
		a5 = -VOLUM0(v8here[5],v8here[6],v8here[0],v8here[4])
		a6 = -VOLUM0(v8here[7],v8here[4],v8here[1],v8here[5])
		a7 = -VOLUM0(v8here[4],v8here[7],v8here[2],v8here[6])
		a8 = -VOLUM0(v8here[6],v8here[5],v8here[3],v8here[7])
		if (a1<=0.0) or (a2<=0.0) or (a3<=0.0) or (a4<=0.0) or (a5<=0.0) or (a6<=0.0) or (a7<=0.0) or (a8<=0.0):
			#print a1,a2,a3,a4,a5,a6,a6,a8
			self.nnrh[ih] = True
			# nnrh: nek non right hand
	
		
	def check_max_aspect_ratio(self):
		print 'check_max_aspect_ratio'
		max_ar = 0
		n_large_ar = 0
		for ih in range(0,self.nhex):
			length = []
			v8 = self.v8[ih]
			
			# edge1
			dist = distance(v8[0],v8[1])
			length.append(dist)
			# edge2
			dist = distance(v8[1],v8[2])
			length.append(dist)
			# edge3
			dist = distance(v8[2],v8[3])
			length.append(dist)
			# edge4
			dist = distance(v8[0],v8[3])
			length.append(dist)
			# edge5
			dist = distance(v8[4],v8[5])
			length.append(dist)
			# edge6
			dist = distance(v8[5],v8[6])
			length.append(dist)
			# edge7
			dist = distance(v8[6],v8[7])
			length.append(dist)
			# edge8
			dist = distance(v8[4],v8[7])
			length.append(dist)
			# edge9
			dist = distance(v8[0],v8[4])
			length.append(dist)
			# edge10
			dist = distance(v8[1],v8[5])
			length.append(dist)
			# edge11
			dist = distance(v8[2],v8[6])
			length.append(dist)
			# edge12
			dist = distance(v8[3],v8[7])
			length.append(dist)
			
			ar = max(length)/min(length)
			if ar>max_ar:
				max_ar = ar
				max_ar_x = v8[0][0]
				max_ar_y = v8[0][1]
				max_ar_z = v8[0][2]
				
			if ar>250:
				n_large_ar = n_large_ar +1
				print 'ar >250 element near ',v8[0][0],',',v8[0][1],',', v8[0][2]
				self.nnrh[ih] = True
				
		print 'max apect ratio: ', max_ar
		print 'near ',max_ar_x,',',max_ar_y,',',max_ar_z
		print n_large_ar,' elements with apect ratio > 250 '
		
	def check_max_aspect_ratio_range(self,start,end,ftag,ip):
		print 'check_max_aspect_ratio'
		max_ar = 0
		n_large_ar = 0
		
		nrhh = []
		for ih in range(start,end):
			length = []
			v8 = self.v8[ih]
			
			# edge1
			dist = distance(v8[0],v8[1])
			length.append(dist)
			# edge2
			dist = distance(v8[1],v8[2])
			length.append(dist)
			# edge3
			dist = distance(v8[2],v8[3])
			length.append(dist)
			# edge4
			dist = distance(v8[0],v8[3])
			length.append(dist)
			# edge5
			dist = distance(v8[4],v8[5])
			length.append(dist)
			# edge6
			dist = distance(v8[5],v8[6])
			length.append(dist)
			# edge7
			dist = distance(v8[6],v8[7])
			length.append(dist)
			# edge8
			dist = distance(v8[4],v8[7])
			length.append(dist)
			# edge9
			dist = distance(v8[0],v8[4])
			length.append(dist)
			# edge10
			dist = distance(v8[1],v8[5])
			length.append(dist)
			# edge11
			dist = distance(v8[2],v8[6])
			length.append(dist)
			# edge12
			dist = distance(v8[3],v8[7])
			length.append(dist)
			
			ar = max(length)/min(length)
			if ar>max_ar:
				max_ar = ar
				max_ar_x = v8[0][0]
				max_ar_y = v8[0][1]
				max_ar_z = v8[0][2]
				
			if ar>200:
				n_large_ar = n_large_ar +1
				#print 'ar >200 element near ',v8[0][0],',',v8[0][1],',', v8[0][2]
				nrhh.append([ih,ar])
				#nrhh.append(ih)

				
		print 'max apect ratio: ', max_ar
		print 'near ',max_ar_x,',',max_ar_y,',',max_ar_z
		print n_large_ar,' elements with apect ratio > 200 '
		
		
		filename = 'dummy_'+ftag+'/dummy_'+str(ip)
		print 'writing files to '+filename
		dummyfile = open(filename,'w')
		nlines = len(nrhh)
		line = str(nlines)+'\n'
		dummyfile.write(line)
		for iline in range(0,nlines):
			line = str(nrhh[iline][0])+' '+ str(nrhh[iline][1])+'\n'
			#line = str(nrhh[iline][0])+'\n'
			dummyfile.write(line)
		dummyfile.close()
		print 'done writing files to '+filename
		
	def check_max_aspect_ratio_one_element(self,ih):

		length = []
		v8 = self.v8[ih]
			
		# edge1
		dist = distance(v8[0],v8[1])
		length.append(dist)
		# edge2
		dist = distance(v8[1],v8[2])
		length.append(dist)
		# edge3
		dist = distance(v8[2],v8[3])
		length.append(dist)
		# edge4
		dist = distance(v8[0],v8[3])
		length.append(dist)
		# edge5
		dist = distance(v8[4],v8[5])
		length.append(dist)
		# edge6
		dist = distance(v8[5],v8[6])
		length.append(dist)
		# edge7
		dist = distance(v8[6],v8[7])
		length.append(dist)
		# edge8
		dist = distance(v8[4],v8[7])
		length.append(dist)
		# edge9
		dist = distance(v8[0],v8[4])
		length.append(dist)
		# edge10
		dist = distance(v8[1],v8[5])
		length.append(dist)
		# edge11
		dist = distance(v8[2],v8[6])
		length.append(dist)
		# edge12
		dist = distance(v8[3],v8[7])
		length.append(dist)
		
		ar = max(length)/min(length)
			
		if ar>250:
			n_large_ar = n_large_ar +1
			print 'ar >250 element near ',v8[0][0],',',v8[0][1],',', v8[0][2]
				
	def fix_non_right_hand_elements(self,ip):
		# fix non right hand element
		print 'fixing non right hand elements'
		nrh = 0
		nrh2 = 0
		
		max_ar = 0
		n_large_ar = 0
		for ih in range(0,self.nhex):
			v12 =[self.v8[ih][1][0]-self.v8[ih][0][0],self.v8[ih][1][1]-self.v8[ih][0][1],self.v8[ih][1][2]-self.v8[ih][0][2]] 
			v14 =[self.v8[ih][3][0]-self.v8[ih][0][0],self.v8[ih][3][1]-self.v8[ih][0][1],self.v8[ih][3][2]-self.v8[ih][0][2]] 
			v15 =[self.v8[ih][4][0]-self.v8[ih][0][0],self.v8[ih][4][1]-self.v8[ih][0][1],self.v8[ih][4][2]-self.v8[ih][0][2]] 
			
			vec1 = cross_vector(v12,v14)
			vdot = dot_vector(vec1,v15)
			
			if vdot < 0 :
				# non_right_hand element, fix it
				nrh = nrh + 1
				v8new = []
				for iv in range(0,4):
					v8new.append(self.v8[ih][iv+4])
				for iv in range(0,4):
					v8new.append(self.v8[ih][iv])
				self.v8[ih] = v8new
				
				e12new = []
				for ie in range(0,4):
					e12new.append(self.e12[ih][ie+4])
				for ie in range(0,4):
					e12new.append(self.e12[ih][ie])
				for ie in range(0,4):
					e12new.append(self.e12[ih][ie+8])		
				self.e12[ih] = e12new
				
				s6new_5 = self.s6[ih][4]
				s6new_6 = self.s6[ih][5]
				
				self.s6[ih][4] = s6new_6
				self.s6[ih][5] = s6new_5
			
			
			# check non-right-elements based on nek method
			
			self.nnrh[ih] = False
			# nek code for non-right-hand check
			#VOLUM0(A,B,C,O) = (A-O)X(B-O).(C-O)
			#V1= VOLUM0(XYZ(1,2,IE),XYZ(1,3,IE),XYZ(1,5,IE),XYZ(1,1,IE))
			#V2= VOLUM0(XYZ(1,4,IE),XYZ(1,1,IE),XYZ(1,6,IE),XYZ(1,2,IE))
			#V3= VOLUM0(XYZ(1,1,IE),XYZ(1,4,IE),XYZ(1,7,IE),XYZ(1,3,IE))
			#V4= VOLUM0(XYZ(1,3,IE),XYZ(1,2,IE),XYZ(1,8,IE),XYZ(1,4,IE))
			#V5=-VOLUM0(XYZ(1,6,IE),XYZ(1,7,IE),XYZ(1,1,IE),XYZ(1,5,IE))
			#V6=-VOLUM0(XYZ(1,8,IE),XYZ(1,5,IE),XYZ(1,2,IE),XYZ(1,6,IE))
			#V7=-VOLUM0(XYZ(1,5,IE),XYZ(1,8,IE),XYZ(1,3,IE),XYZ(1,7,IE))
			#V8=-VOLUM0(XYZ(1,7,IE),XYZ(1,6,IE),XYZ(1,4,IE),XYZ(1,8,IE))
			#
			#IF (V1.LE.0.0.OR.V2.LE.0.0.OR.V3.LE.0.0.OR.V4.LE.0.0.OR.V5.LE.0.0.OR.V6.LE.0.0.OR.V7.LE.0.0.OR.V8.LE.0.0    ) THEN
			
			# swap here
			v8here = []
			v8here.append(self.v8[ih][0])
			v8here.append(self.v8[ih][1])
			v8here.append(self.v8[ih][3])
			v8here.append(self.v8[ih][2])
			v8here.append(self.v8[ih][4])
			v8here.append(self.v8[ih][5])
			v8here.append(self.v8[ih][7])
			v8here.append(self.v8[ih][6])
			
			a1 = VOLUM0(v8here[1],v8here[2],v8here[4],v8here[0])
			a2 = VOLUM0(v8here[3],v8here[0],v8here[5],v8here[1])
			a3 = VOLUM0(v8here[0],v8here[3],v8here[6],v8here[2])
			a4 = VOLUM0(v8here[2],v8here[1],v8here[7],v8here[3])
			a5 = -VOLUM0(v8here[5],v8here[6],v8here[0],v8here[4])
			a6 = -VOLUM0(v8here[7],v8here[4],v8here[1],v8here[5])
			a7 = -VOLUM0(v8here[4],v8here[7],v8here[2],v8here[6])
			a8 = -VOLUM0(v8here[6],v8here[5],v8here[3],v8here[7])
			if (a1<=0.0) or (a2<=0.0) or (a3<=0.0) or (a4<=0.0) or (a5<=0.0) or (a6<=0.0) or (a7<=0.0) or (a8<=0.0):
				#print a1,a2,a3,a4,a5,a6,a6,a8
				nrh2 = nrh2 + 1
				self.nnrh[ih] = True
				# nnrh: nek non right hand
			
			
			# calculate aspect ratio, and collect max aspect ratio values
			length = []
			v8 = self.v8[ih]
			
			# edge1
			dist = distance(v8[0],v8[1])
			length.append(dist)
			# edge2
			dist = distance(v8[1],v8[2])
			length.append(dist)
			# edge3
			dist = distance(v8[2],v8[3])
			length.append(dist)
			# edge4
			dist = distance(v8[0],v8[3])
			length.append(dist)
			# edge5
			dist = distance(v8[4],v8[5])
			length.append(dist)
			# edge6
			dist = distance(v8[5],v8[6])
			length.append(dist)
			# edge7
			dist = distance(v8[6],v8[7])
			length.append(dist)
			# edge8
			dist = distance(v8[4],v8[7])
			length.append(dist)
			# edge9
			dist = distance(v8[0],v8[4])
			length.append(dist)
			# edge10
			dist = distance(v8[1],v8[5])
			length.append(dist)
			# edge11
			dist = distance(v8[2],v8[6])
			length.append(dist)
			# edge12
			dist = distance(v8[3],v8[7])
			length.append(dist)
			
			ar = max(length)/min(length)
			if ar>max_ar:
				max_ar = ar
				max_ar_x = v8[0][0]
				max_ar_y = v8[0][1]
				max_ar_z = v8[0][2]
			
			self.lar[ih] = False
			if ar>250:
				n_large_ar = n_large_ar +1
				print 'ar >250 element near ',v8[0][0],',',v8[0][1],',', v8[0][2],' in proc' +str(ip) + '\n'
				self.lar[ih] = True
		
		#print '============================================'
		#print str(nrh)+' non right hand elements are fixed'
		#print '============================================'
	
		#print '============================================'
		if nrh2>0: print str(nrh2)+' non right hand elements detected in proc' +str(ip) + '\n'
		#print '============================================'
		
		#print '============================================'
		if n_large_ar > 0:
			print str(n_large_ar) + ' elements with apect ratio > 250 in proc' +str(ip)  + '\n'
			print 'max aspect ratio: '+ str(max_ar) + ' in proc' +str(ip)  + '\n'
		#print '============================================'
	
	def fix_non_right_hand_elements_range(self,start,end):
		# fix non right hand element
		print 'fixing non right hand elements'
		nrh = 0
		for ih in range(start,end):
			v12 =[self.v8[ih][1][0]-self.v8[ih][0][0],self.v8[ih][1][1]-self.v8[ih][0][1],self.v8[ih][1][2]-self.v8[ih][0][2]] 
			v14 =[self.v8[ih][3][0]-self.v8[ih][0][0],self.v8[ih][3][1]-self.v8[ih][0][1],self.v8[ih][3][2]-self.v8[ih][0][2]] 
			v15 =[self.v8[ih][4][0]-self.v8[ih][0][0],self.v8[ih][4][1]-self.v8[ih][0][1],self.v8[ih][4][2]-self.v8[ih][0][2]] 
			
			vec1 = np.cross(v12,v14)
			vdot = np.dot(vec1,v15)
			
			if vdot < 0 :
				# non_right_hand element, fix it
				nrh = nrh + 1
				v8new = []
				for iv in range(0,4):
					v8new.append(self.v8[ih][iv+4])
				for iv in range(0,4):
					v8new.append(self.v8[ih][iv])
				self.v8[ih] = v8new
				
				e12new = []
				for ie in range(0,4):
					e12new.append(self.e12[ih][ie+4])
				for ie in range(0,4):
					e12new.append(self.e12[ih][ie])
				for ie in range(0,4):
					e12new.append(self.e12[ih][ie+8])		
				self.e12[ih] = e12new
				
				s6new_5 = self.s6[ih][4]
				s6new_6 = self.s6[ih][5]
				
				self.s6[ih][4] = s6new_6
				self.s6[ih][5] = s6new_5
		print '============================================'
		print str(nrh)+' non right hand elements are fixed'
		print '============================================'

	def return_face_info(self,e,f):
		# return face center xyz, face normal,face vertices
		#
		#
		v8 = self.v8[e]
		e12 = self.e12[e]
		s6 = self.s6[e]
		
		if f == 4:
			v1 = v8[0]
			v2 = v8[1]
			v3 = v8[2]
			v4 = v8[3]			
			vc = four_to_one(v1,v2,v3,v4)
			fn = get_v1_normal(v1,v4,v2)
			vquad = [v1,v2,v3,v4]
			
			
			e1 = e12[0]
			e2 = e12[1]
			e3 = e12[2]
			e4 = e12[3]
			equad = [e1,e2,e3,e4]
			
			s1 = s6[0]
			s2 = s6[1]
			s3 = s6[2]
			s4 = s6[3]
			ftag = [s1,s2,s3,s4]
		if f == 5:
			v1 = v8[4]
			v2 = v8[5]
			v3 = v8[6]
			v4 = v8[7]
			vc = four_to_one(v1,v2,v3,v4)
			fn = get_v1_normal(v1,v2,v4)
			vquad = [v1,v2,v3,v4]
			
			e1 = e12[4]
			e2 = e12[5]
			e3 = e12[6]
			e4 = e12[7]
			equad = [e1,e2,e3,e4]
			
			s1 = s6[0]
			s2 = s6[1]
			s3 = s6[2]
			s4 = s6[3]
			ftag = [s1,s2,s3,s4]
			
		if f == 0:
			v1 = v8[0]
			v2 = v8[1]
			v3 = v8[5]
			v4 = v8[4]
			vc = four_to_one(v1,v2,v3,v4)
			fn = get_v1_normal(v1,v2,v4)
			vquad = [v1,v2,v3,v4]
			
			e1 = e12[0]
			e2 = e12[9]
			e3 = e12[4]
			e4 = e12[8]
			equad = [e1,e2,e3,e4]
			
			s1 = s6[4]
			s2 = s6[1]
			s3 = s6[5]
			s4 = s6[3]
			ftag = [s1,s2,s3,s4]
			
		if f == 1:
			v1 = v8[1]
			v2 = v8[2]
			v3 = v8[6]
			v4 = v8[5]
			vc = four_to_one(v1,v2,v3,v4)
			fn = get_v1_normal(v1,v2,v4)
			vquad = [v1,v2,v3,v4]
			
		
			e1 = e12[1]
			e2 = e12[10]
			e3 = e12[5]
			e4 = e12[9]
			equad = [e1,e2,e3,e4]
			
			s1 = s6[4]
			s2 = s6[2]
			s3 = s6[5]
			s4 = s6[0]
			ftag = [s1,s2,s3,s4]

			
		if f == 2:
			v1 = v8[2]
			v2 = v8[3]
			v3 = v8[7]
			v4 = v8[6]
			vc = four_to_one(v1,v2,v3,v4)
			fn = get_v1_normal(v1,v2,v4)
			vquad = [v1,v2,v3,v4]
			
					
			e1 = e12[2]
			e2 = e12[11]
			e3 = e12[6]
			e4 = e12[10]
			equad = [e1,e2,e3,e4]
			
			
			s1 = s6[4]
			s2 = s6[3]
			s3 = s6[5]
			s4 = s6[1]
			ftag = [s1,s2,s3,s4]
			
		if f == 3:
			v1 = v8[3]
			v2 = v8[0]
			v3 = v8[4]
			v4 = v8[7]
			vc = four_to_one(v1,v2,v3,v4)
			fn = get_v1_normal(v1,v2,v4)
			vquad = [v1,v2,v3,v4]
						
			e1 = e12[3]
			e2 = e12[8]
			e3 = e12[7]
			e4 = e12[11]
			equad = [e1,e2,e3,e4]
			
						
			s1 = s6[4]
			s2 = s6[0]
			s3 = s6[5]
			s4 = s6[2]
			ftag = [s1,s2,s3,s4]
			
		return vc,fn,vquad,equad,ftag
		
		
		