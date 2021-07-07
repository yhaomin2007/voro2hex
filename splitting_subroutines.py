import numpy as np
import math

def tri_to_quad1(tri,pxyz,fc,diameter,ptag,all_pebbles):
	# try[3][3] triangle 3 points xyz
	# pxyz[2][3] two pebbles xyz
	#   pebble diameter

	m12 = line_split(tri[0],tri[1],0.5)
	m23 = line_split(tri[1],tri[2],0.5)
	m13 = line_split(tri[0],tri[2],0.5)

	#m12,moved = move_away_from_pebbles_global(m12,all_pebbles,diameter,1.03,ptag)
	#m23,moved = move_away_from_pebbles_global(m23,all_pebbles,diameter,1.03,ptag)
	#m13,moved = move_away_from_pebbles_global(m13,all_pebbles,diameter,1.03,ptag)
	
	quads = []
	quads.append([tri[0],m12,fc,m13])
	quads.append([tri[1],m23,fc,m12])
	quads.append([tri[2],m13,fc,m23])
	
	# quads[3][4][3]
	# 3 quads, 4 points per quad, xyz for each point
	return quads
	
def tri_to_quad2(tri,pxyz,fc,diameter,ifInternal,ptag,all_pebbles):
	# try[3][3] triangle 3 points xyz
	# pxyz[2][3] two pebbles xyz
	#   pebble diameter

	m12 = line_split(tri[0],tri[1],0.5)
	m23 = line_split(tri[1],tri[2],0.5)
	m13 = line_split(tri[0],tri[2],0.5)
	
	if ifInternal:
		m12,moved = move_away_from_pebbles(m12,pxyz,diameter,1.03,ptag)
		#m23,moved = move_away_from_pebbles_global(m23,all_pebbles,diameter,1.03,ptag)
		m13,moved = move_away_from_pebbles(m13,pxyz,diameter,1.03,ptag)
	
	quads = []
	quads.append([tri[0],m12,fc,m13])
	quads.append([tri[1],m23,fc,m12])
	quads.append([tri[2],m13,fc,m23])
	
	# quads[3][4][3]
	# 3 quads, 4 points per quad, xyz for each point
	return quads

def rec_to_quad(rec,pxyz,fc,diameter,ptag,all_pebbles):

	m12 =  line_split(rec[0],rec[1],0.5)
	m23 =  line_split(rec[1],rec[2],0.5)
	m34 =  line_split(rec[2],rec[3],0.5)
	m14 =  line_split(rec[3],rec[0],0.5)
	
	#m12,moved = move_away_from_pebbles_global(m12,all_pebbles,diameter,1.03,ptag)
	#m23,moved = move_away_from_pebbles_global(m23,all_pebbles,diameter,1.03,ptag)
	#m34,moved = move_away_from_pebbles_global(m34,all_pebbles,diameter,1.03,ptag)
	#m14,moved = move_away_from_pebbles_global(m14,all_pebbles,diameter,1.03,ptag)
	
	quads = []
	
	quads.append([rec[0],m12,fc,m14])
	quads.append([rec[1],m23,fc,m12])
	quads.append([rec[2],m34,fc,m23])
	quads.append([rec[3],m14,fc,m34])
	
	return quads
	

def rec_to_quad_for_chamfer(rec,vcm,pxyz,diameter,ifInternal,ptag,all_pebbles):

	m12 =  line_split(rec[0],rec[1],0.5)
	#m23 =  line_split(rec[1],rec[2],0.5)
	m23 =  vcm
	m34 =  line_split(rec[2],rec[3],0.5)
	m14 =  line_split(rec[3],rec[0],0.5)
	
	if ifInternal:
		m12,moved = move_away_from_pebbles(m12,pxyz,diameter,1.03,ptag)
		m23,moved = move_away_from_pebbles(m23,pxyz,diameter,1.005,ptag)
		m34,moved = move_away_from_pebbles(m34,pxyz,diameter,1.03,ptag)
		#m14,moved = move_away_from_pebbles_global(m14,all_pebbles,diameter,1.03,ptag)
	
	fc = line_split(m23,m14,0.5)
	fc,moved = move_away_from_pebbles(fc,pxyz,diameter,1.025,ptag)
	
	quads = []
	
	quads.append([rec[0],m12,fc,m14])
	quads.append([rec[1],m23,fc,m12])
	quads.append([rec[2],m34,fc,m23])
	quads.append([rec[3],m14,fc,m34])
	
	return quads
	
def rec_to_quad_for_chamfer2(rec,vcm,pxyz,diameter,ifInternal,ptag,all_pebbles):
	# split more in edge to chamfer directions

	m12 =  line_split(rec[0],rec[1],0.5)
	#m23 =  line_split(rec[1],rec[2],0.5)
	m23 =  vcm
	m34 =  line_split(rec[2],rec[3],0.5)
	m14 =  line_split(rec[3],rec[0],0.5)
	
	if ifInternal:
		m12,moved = move_away_from_pebbles(m12,pxyz,diameter,1.01,ptag)
		#m23,moved = move_away_from_pebbles(m23,pxyz,diameter,1.005,ptag)
		m34,moved = move_away_from_pebbles(m34,pxyz,diameter,1.01,ptag)
		#m14,moved = move_away_from_pebbles_global(m14,all_pebbles,diameter,1.03,ptag)
	
	fc = line_split(m23,m14,0.5)
	fc,moved = move_away_from_pebbles(fc,pxyz,diameter,1.01,ptag)
	
	p0 = line_split(rec[0],m12,0.5)
	p1 = line_split(fc,m14,0.5)
	p2 = line_split(rec[3],m34,0.5)
	
	p3 = line_split(rec[1],m12,0.5)
	p4 = line_split(fc,m23,0.5)
	p5 = line_split(rec[2],m34,0.5)
	
	p0,moved = move_away_from_pebbles(p0,pxyz,diameter,1.01,ptag)
	p1,moved = move_away_from_pebbles(p1,pxyz,diameter,1.01,ptag)
	p2,moved = move_away_from_pebbles(p2,pxyz,diameter,1.01,ptag)
	p3,moved = move_away_from_pebbles(p3,pxyz,diameter,1.01,ptag)
	p4,moved = move_away_from_pebbles(p4,pxyz,diameter,1.01,ptag)
	p5,moved = move_away_from_pebbles(p5,pxyz,diameter,1.01,ptag)
		
	quads = []
	
	quads.append([rec[0],p0,p1,m14])
	quads.append([p0,m12,fc,p1])
	quads.append([rec[1],m23,p4,p3])
	quads.append([p3,p4,fc,m12])
	quads.append([rec[2],p5,p4,m23])
	quads.append([p5,m34,fc,p4])
	quads.append([rec[3],m14,p1,p2])
	quads.append([p2,p1,fc,m34])
	
	return quads

def rec_to_quad_for_chamfer3(rec,vcm,pxyz,diameter,ifInternal,ptag,all_pebbles):
	# split more in edge to chamfer directions

	m12 =  line_split(rec[0],rec[1],0.5)
	#m23 =  line_split(rec[1],rec[2],0.5)
	m23 =  vcm
	m34 =  line_split(rec[2],rec[3],0.5)
	m14 =  line_split(rec[3],rec[0],0.5)
	
	if ifInternal:
		m12,moved = move_away_from_pebbles(m12,pxyz,diameter,1.01,ptag)
		#m23,moved = move_away_from_pebbles(m23,pxyz,diameter,1.005,ptag)
		m34,moved = move_away_from_pebbles(m34,pxyz,diameter,1.01,ptag)
		#m14,moved = move_away_from_pebbles_global(m14,all_pebbles,diameter,1.03,ptag)
	
	fc = line_split(m23,m14,0.5)
	fc,moved = move_away_from_pebbles(fc,pxyz,diameter,1.01,ptag)
	
	p0 = line_split(rec[0],m12,0.5)
	p1 = line_split(fc,m14,0.5)
	p2 = line_split(rec[3],m34,0.5)
	
	p3 = line_split(rec[1],m12,0.5)
	p4 = line_split(fc,m23,0.5)
	p5 = line_split(rec[2],m34,0.5)
	
	p0,moved = move_away_from_pebbles(p0,pxyz,diameter,1.01,ptag)
	p1,moved = move_away_from_pebbles(p1,pxyz,diameter,1.01,ptag)
	p2,moved = move_away_from_pebbles(p2,pxyz,diameter,1.01,ptag)
	p3,moved = move_away_from_pebbles(p3,pxyz,diameter,1.01,ptag)
	p4,moved = move_away_from_pebbles(p4,pxyz,diameter,1.01,ptag)
	p5,moved = move_away_from_pebbles(p5,pxyz,diameter,1.01,ptag)
	
	
	p6 = line_split(rec[1],p3,0.3)
	p7 = line_split(m23,p4,0.3)
	p8 = line_split(rec[2],p5,0.3)
	
	p6,moved = move_away_from_pebbles(p6,pxyz,diameter,1.005,ptag)
	p7,moved = move_away_from_pebbles(p7,pxyz,diameter,1.005,ptag)
	p8,moved = move_away_from_pebbles(p8,pxyz,diameter,1.005,ptag)
	
		
	quads = []
	
	quads.append([rec[0],p0,p1,m14])
	quads.append([p0,m12,fc,p1])
	quads.append([rec[1],m23,p7,p6])
	quads.append([p6,p7,p4,p3])
	quads.append([p3,p4,fc,m12])
	
	quads.append([rec[2],p8,p7,m23])
	quads.append([p8,p5,p4,p7])
	quads.append([p5,m34,fc,p4])
	quads.append([rec[3],m14,p1,p2])
	quads.append([p2,p1,fc,m34])
	
	return quads
	
def quad_to_4quads(quad,tag,pxyz,diameter,ptag):
	# 
	four_quads = []
	
	m12 =  line_split(quad[0],quad[1],0.5)
	m23 =  line_split(quad[1],quad[2],0.5)
	m34 =  line_split(quad[2],quad[3],0.5)
	m14 =  line_split(quad[3],quad[0],0.5)
	fc = line_split(m12,m34,0.5)
		
	m12,moved = move_away_from_pebbles(m12,pxyz,diameter,1.02,ptag)
	m23,moved = move_away_from_pebbles(m23,pxyz,diameter,1.02,ptag)
	m34,moved = move_away_from_pebbles(m34,pxyz,diameter,1.02,ptag)
	m14,moved = move_away_from_pebbles(m14,pxyz,diameter,1.02,ptag)
	fc,moved = move_away_from_pebbles(fc,pxyz,diameter,1.02,ptag)
	
	four_quads.append([quad[0],m12,fc,m14])
	four_quads.append([quad[1],m23,fc,m12])
	four_quads.append([quad[2],m34,fc,m23])
	four_quads.append([quad[3],m14,fc,m34])
	
	four_tags = []
	four_tags.append([tag[0],'E  ','E  ',tag[3]])
	four_tags.append([tag[1],'E  ','E  ',tag[0]])
	four_tags.append([tag[2],'E  ','E  ',tag[1]])
	four_tags.append([tag[3],'E  ','E  ',tag[2]])
	
	return four_quads,four_tags
	
def quad_to_4quads_for_chamfer(quad,tag,pxyz,diameter,r_chamfer,ptag):
	# 
	four_quads = []
	
	r_c = r_chamfer*(diameter/2.0)
	
	m12 =  line_split(quad[0],quad[1],0.5)
	m23 =  line_split(quad[1],quad[2],0.5)
	m34 =  line_split(quad[2],quad[3],0.5)
	m14 =  line_split(quad[3],quad[0],0.5)
	
	#pmiddle = line_split(pxyz[0],pxyz[1],0.5)
	#
	#cindex = tag.index('C  ')
	#if cindex == 0:
	#	d = distance(m34,pmiddle)
	#	if r_c > 0.95*d: r_c = 0.95*d
	#	m12 =  line_split(pmiddle,m34,r_c/d)
	#	fc = line_split(m12,m34,0.5)
	#elif cindex == 1:
	#	d = distance(m14,pmiddle)
	#	if r_c > 0.95*d: r_c = 0.95*d
	#	m23 =  line_split(pmiddle,m14,r_c/d)
	#	fc = line_split(m23,m14,0.5)
	#elif cindex == 2:
	#	d = distance(m12,pmiddle)
	#	if r_c > 0.95*d: r_c = 0.95*d
	#	m34 =  line_split(pmiddle,m12,r_c/d)
	#	fc = line_split(m12,m34,0.5)
	#elif cindex == 3:
	#	d = distance(m23,pmiddle)
	#	if r_c > 0.95*d: r_c = 0.95*d
	#	m14 =  line_split(pmiddle,m23,r_c/d)
	#	fc = line_split(m23,m14,0.5)
	#
	#m12,moved = move_away_from_pebbles(m12,pxyz,diameter,1.02,ptag)
	#m23,moved = move_away_from_pebbles(m23,pxyz,diameter,1.02,ptag)
	#m34,moved = move_away_from_pebbles(m34,pxyz,diameter,1.02,ptag)
	#m14,moved = move_away_from_pebbles(m14,pxyz,diameter,1.02,ptag)
	#fc,moved = move_away_from_pebbles(fc,pxyz,diameter,1.02,ptag)
	
	four_quads.append([quad[0],m12,fc,m14])
	four_quads.append([quad[1],m23,fc,m12])
	four_quads.append([quad[2],m34,fc,m23])
	four_quads.append([quad[3],m14,fc,m34])
	
	four_tags = []
	four_tags.append([tag[0],'E  ','E  ',tag[3]])
	four_tags.append([tag[1],'E  ','E  ',tag[0]])
	four_tags.append([tag[2],'E  ','E  ',tag[1]])
	four_tags.append([tag[3],'E  ','E  ',tag[2]])
	
	return four_quads,four_tags
	
def edge_pebble_distance(v1xyz,v2xyz,pxyz):
	# calculate distance between pxyz and edge(v1-v2)
	
	v12 = vector_minus(v2xyz,v1xyz)
	vp1 = vector_minus(v1xyz,pxyz)
	
	norm_v12 = np.linalg.norm(v12)
	
	alpha = - np.dot(v12,vp1)/norm_v12**2.0
	
	if alpha < 0 or alpha > 1: 
		inside = False
	else:
		inside = True
		
	pd = [0,0,0]
	pd[0] = v1xyz[0]+alpha*v12[0]
	pd[1] = v1xyz[1]+alpha*v12[1]
	pd[2] = v1xyz[2]+alpha*v12[2]
	
	dist = distance(pd,pxyz)
	
	nvec = vector_minus(pd,pxyz)
	nvec[0] = nvec[0]/dist
	nvec[1] = nvec[1]/dist
	nvec[2] = nvec[2]/dist
	
	return dist,nvec,inside
	
def move_away_from_pebbles(xyz,ppxyz,diameter,r1,ptag):

	#return xyz,False
	## move away from pebbles
	## this is a very brutal force approach
	
	dp1 = distance(xyz,ppxyz[0])
	dp2 = distance(xyz,ppxyz[1])
	dp_avg = 0.5*(dp1+dp2)
	dp = [dp1,dp2]
		
	if min(dp) > r1*(diameter/2.0):
		return xyz, False
	else:
		if ptag[0]== '' and ptag[1] == '':
		# if internal face/quad
			if dp1 <= dp2:
				pxyz = ppxyz[0]
			else:
				pxyz = ppxyz[1]

			nvc = [xyz[0]-pxyz[0],xyz[1]-pxyz[1],xyz[2]-pxyz[2]]
			nvc_length = np.linalg.norm(nvc)
			nvc = [nvc[0]/nvc_length,nvc[1]/nvc_length,nvc[2]/nvc_length]
			new_dp = (diameter/2.0)*r1
			#new_dp =  dp_avg
			new = [pxyz[0]+nvc[0]*new_dp,pxyz[1]+nvc[1]*new_dp,pxyz[2]+nvc[2]*new_dp]

			return new, True
		
		else:
	
			return xyz, False

def move_away_from_pebbles_global(xyz,all_pebbles,diameter,r1,ptag):

	#return xyz,False
	## move away from all pebbles
	## this is a very brutal force approach,
	# search all pebbles
	
	if ptag[0]== '' and ptag[1] == '':
		minDist = 100.0
		number_of_pebbles = int(all_pebbles.rpebbles)
		#print number_of_pebbles
		for ipebble in range(0,number_of_pebbles):
		
			pxyz = all_pebbles.xyz[ipebble]
			dist = distance(xyz,pxyz)
			if dist < minDist:
				minDist = dist
				min_pxyz = pxyz
		
		if minDist > r1*(diameter/2.0):
			return xyz, False
		else:
			pxyz = min_pxyz
			nvc = [xyz[0]-pxyz[0],xyz[1]-pxyz[1],xyz[2]-pxyz[2]]
			nvc_length = np.linalg.norm(nvc)
			nvc = [nvc[0]/nvc_length,nvc[1]/nvc_length,nvc[2]/nvc_length]
			new_dp = (diameter/2.0)*r1
			#new_dp =  dp_avg
			new = [pxyz[0]+nvc[0]*new_dp,pxyz[1]+nvc[1]*new_dp,pxyz[2]+nvc[2]*new_dp]

			return new, True
	else:
		return xyz, False

def move_away_from_chamfer(xyz,ppxyz,diameter,r_chamfer):
	# move away from chamfer
	cxyz= [0,0,0]
	for i in range(0,3):
		cxyz[i] = (ppxyz[0][i]+ppxyz[1][i])/2.0
	
	cradius = r_chamfer*diameter/2.0
	
	dist = distance(xyz,cxyz)
	
	if dist < 1.1*cradius:
		
		nvc = [xyz[0]-cxyz[0],xyz[1]-cxyz[1],xyz[2]-cxyz[2]]
		nvc_length = np.linalg.norm(nvc)
		nvc = [nvc[0]/nvc_length,nvc[1]/nvc_length,nvc[2]/nvc_length]
		new_radius = 1.1*cradius
		new = [cxyz[0]+nvc[0]*new_radius,cxyz[1]+nvc[1]*new_radius,cxyz[2]+nvc[2]*new_radius]
		
		return new, True
	else:
		return xyz, False
		

def project_to_pebble_pebble_midplane(xyz,ppxyz,diameter):
	
		
	dr,dvec = distance_to_pebble_pebble_midline(xyz,ppxyz,diameter)
	

	d = dr/(diameter/2.0)
	R = [xyz[0]+dvec[0],xyz[1]+dvec[1],xyz[2]+dvec[2]]

	return R

		
def distance_to_pebble_pebble_midline(xyz,ppxyz,diameter):
	
	pmid = line_split(ppxyz[0],ppxyz[1],0.5)
	
	p12 = vector_minus(ppxyz[1],ppxyz[0])
	d12 = np.linalg.norm(p12)
	nv12 = [p12[0]/d12,p12[1]/d12,p12[2]/d12]
	
	theta = vector_minus(pmid,xyz)
	alpha = -np.dot(theta,nv12) 
	pp =   [pmid[0]+alpha*nv12[0],pmid[1]+alpha*nv12[1],pmid[2]+alpha*nv12[2]]
	dr = distance(pp,xyz)
	
	dvec = [-alpha*nv12[0],-alpha*nv12[1],-alpha*nv12[2]]
	# dr is distance
	# dvec is projection vector to pebble-pebble midplane
	return dr,dvec
		
def line_mid_project_to_sphere(v1,v2,p,alpha):
	# v1[3],v2[3].p[3],
	# v3 = v1*alpha + v2*(1-alpha)
	#
	# when alpha = 1, v3 = v1
	# when alpha = 0, v3 = v2
 	beta = 1 - alpha

	pv1 = [v1[0]-p[0],v1[1]-p[1],v1[2]-p[2]]
	pv2 = [v2[0]-p[0],v2[1]-p[1],v2[2]-p[2]]
	r1 = np.linalg.norm(pv1)
	r2 = np.linalg.norm(pv2)
	pv1n = [pv1[0]/r1,pv1[1]/r1,pv1[2]/r1]
	pv2n = [pv2[0]/r2,pv2[1]/r2,pv2[2]/r2]

	pv3 = [pv1n[0]*alpha+pv2n[0]*beta,pv1n[1]*alpha+pv2n[1]*beta,pv1n[2]*alpha+pv2n[2]*beta]
	r3 = np.linalg.norm(pv3)
	pv3n = [pv3[0]/r3,pv3[1]/r3,pv3[2]/r3]

	r3 = r1*alpha + r2*beta
	mid = [pv3n[0]*r3+p[0],pv3n[1]*r3+p[1],pv3n[2]*r3+p[2]]

	return mid

def project_to_pebble_ratio(quad,pxyz,pdiameter,ratio,tol):
	quad_on_p = []
	for ipoint in range(0,4):
		d = ((quad[ipoint][0]-pxyz[0])**2.0+(quad[ipoint][1]-pxyz[1])**2.0+(quad[ipoint][2]-pxyz[2])**2.0)**0.5
		r = (pdiameter/2.0)*(1-ratio) + d*(ratio)
		
		if d <= (pdiameter/2.0+tol): r = d - tol
		
		x_proj = pxyz[0]+ (r/d)*(quad[ipoint][0]-pxyz[0])
		y_proj = pxyz[1]+ (r/d)*(quad[ipoint][1]-pxyz[1])
		z_proj = pxyz[2]+ (r/d)*(quad[ipoint][2]-pxyz[2])
		
		quad_on_p.append([x_proj,y_proj,z_proj])

	# quad_on_p[3]
	return quad_on_p

def project_to_pebble(quad,pxyz,pdiameter,tol):
	# quad[4][3] quad xyz 
	# pxyz[3]	 pebble xyz
	
	quad_on_p = []
	for ipoint in range(0,4):
		d = ((quad[ipoint][0]-pxyz[0])**2.0+(quad[ipoint][1]-pxyz[1])**2.0+(quad[ipoint][2]-pxyz[2])**2.0)**0.5
		r = pdiameter/2.0
		
		if d <= (r+tol): r = d - tol
		#if d < r: 
		#	print 'ERROR: negative projection'
		#	print quad[ipoint]
		
		x_proj = pxyz[0]+ (r/d)*(quad[ipoint][0]-pxyz[0])
		y_proj = pxyz[1]+ (r/d)*(quad[ipoint][1]-pxyz[1])
		z_proj = pxyz[2]+ (r/d)*(quad[ipoint][2]-pxyz[2])
		
		quad_on_p.append([x_proj,y_proj,z_proj])

	# quad_on_p[3]
	return quad_on_p

def adjust_top_quad_for_chamfer(quad_bot,quad_top,four_edge_tags,mid,pxyz):
	
	m2p = vector_minus(pxyz,mid)
	norm = np.linalg.norm(m2p)
	
	dist = 0.005
	m2p = [dist*m2p[0]/norm,dist*m2p[1]/norm,dist*m2p[2]/norm]
	

	edge_id = four_edge_tags.index('C  ')
	
	v1 = edge_id
	v2 = v1 + 1
	if v2 > 3: v2 = 0

	adjusted_top_quad = quad_top
	adjusted_top_quad[v1] = offset_vert(quad_bot[v1],m2p)
	adjusted_top_quad[v2] = offset_vert(quad_bot[v2],m2p)
	
	return adjusted_top_quad
	
def ifQuadsShareTwoNodes(quad1,quad2):
	
	ifShare = False
	
	shareNode = 0
	
	for iv1 in range(0,4):
		xyz1 = quad1[iv1]
		for iv2 in range(0,4):
			xyz2 = quad2[iv2]
			dist = distance(xyz1,xyz2)
			if dist < 1e-5: shareNode = shareNode + 1
	
	if shareNode == 2: ifShare= True
	return ifShare
	
def line_split(xyz1,xyz2,alpha):
	# xyz1[3]
	# xyz2[3]
	# bl 
	#
	# return xyz3
	#
	# 1 <-alpha-> 3 <- (1-alpha)->2
	
	x3 = xyz1[0]+ alpha*(xyz2[0]-xyz1[0])
	y3 = xyz1[1]+ alpha*(xyz2[1]-xyz1[1])
	z3 = xyz1[2]+ alpha*(xyz2[2]-xyz1[2])
	
	xyz3 = [x3,y3,z3]
	
	return xyz3
	
def distance(xyz1,xyz2):
	d = ((xyz1[0]-xyz2[0])**2.0+(xyz1[1]-xyz2[1])**2.0+(xyz1[2]-xyz2[2])**2.0)**0.5
	
	return d
		
def line_split_for_boundary_layer(xyz1,xyz2,bl):
	# xyz1[3]
	# xyz2[3]
	# bl 
	#
	# return xyz3
	#
	# 1 <-alpha-> 3 <- (1-alpha)->2
	length = ((xyz1[0]-xyz2[0])**2.0+(xyz1[1]-xyz2[1])**2.0+(xyz1[2]-xyz2[2])**2.0)**0.5
	if bl > 0.5*length: bl = 0.5*length
	
	x3 = xyz1[0]+ (bl/length)*(xyz2[0]-xyz1[0])
	y3 = xyz1[1]+ (bl/length)*(xyz2[1]-xyz1[1])
	z3 = xyz1[2]+ (bl/length)*(xyz2[2]-xyz1[2])
	
	xyz3 = [x3,y3,z3]
	
	return xyz3
	
def offset_vert(xyz1,offset):
	# xyz1[3]
	# return xyz2 = xyz1 + offset 
	x2 = xyz1[0]+ offset[0]
	y2 = xyz1[1]+ offset[1]
	z2 = xyz1[2]+ offset[2]
	xyz2 = [x2,y2,z2]
	return xyz2

def proj_vertice_to_cyl(xyz,cly_radius):
	xx = xyz[0]
	yy = xyz[1]
	zz = xyz[2]
	rr =(xx**2.0+yy**2.0)**0.5
	xx_new = cly_radius*xx/rr
	yy_new = cly_radius*yy/rr
	new = [xx_new,yy_new,zz]
	
	return new

def v8_to_e12(v8):
	# calcualte e12 based v8
	e12 = []
	for i in range(0,4):
		i1 = i
		i2 = i+1
		if i2 > 3:
			i2 = 0
		midx = (v8[i1][0] + v8[i2][0])/2.0
		midy = (v8[i1][1] + v8[i2][1])/2.0
		midz = (v8[i1][2] + v8[i2][2])/2.0
		e12.append([midx,midy,midz])
	
	for i in range(4,8):
		i1 = i
		i2 = i+1
		if i2 > 7:
			i2 = 4
		midx = (v8[i1][0] + v8[i2][0])/2.0
		midy = (v8[i1][1] + v8[i2][1])/2.0
		midz = (v8[i1][2] + v8[i2][2])/2.0
		e12.append([midx,midy,midz])
		
	for i in range(0,4):
		i1 = i
		i2 = i1+4
		midx = (v8[i1][0] + v8[i2][0])/2.0
		midy = (v8[i1][1] + v8[i2][1])/2.0
		midz = (v8[i1][2] + v8[i2][2])/2.0
		e12.append([midx,midy,midz])
	
	return e12

def four_to_one(v1,v2,v3,v4):
	
	vcx = (v1[0]+v2[0]+v3[0]+v4[0])*0.25
	vcy = (v1[1]+v2[1]+v3[1]+v4[1])*0.25
	vcz = (v1[2]+v2[2]+v3[2]+v4[2])*0.25
	vc = [vcx,vcy,vcz]	
	return vc
	
def three_to_one(v1,v2,v3):
	
	vcx = (v1[0]+v2[0]+v3[0])*1.0/3.0
	vcy = (v1[1]+v2[1]+v3[1])*1.0/3.0
	vcz = (v1[2]+v2[2]+v3[2])*1.0/3.0
	vc = [vcx,vcy,vcz]	
	return vc

def get_v1_normal(v1,v2,v4):
	vec12 = [v2[0]-v1[0],v2[1]-v1[1],v2[2]-v1[2]]
	vec14 = [v4[0]-v1[0],v4[1]-v1[1],v4[2]-v1[2]]
	vcross =  np.cross(vec12,vec14)
	norm = np.linalg.norm(vcross)
	v1n = [vcross[0]/norm,vcross[1]/norm,vcross[2]/norm]
	return v1n
	
def merge_two_vertices(v1,v2,vert):
	# merge_two_vertices to middle
	xyz1 = vert.xyz[v1]
	xyz2 = vert.xyz[v2]
	mid = line_split(xyz1,xyz2,0.5)
	vert.xyz[v1] = mid
	vert.xyz[v2] = mid

	for iv in range(0,vert.nv):
		if iv != v1 and iv != v2:
			dist1 = distance(xyz1,vert.xyz[iv])
			dist2 = distance(xyz2,vert.xyz[iv])
			if(dist1 < 1e-5) or (dist2 < 1e-5):
				vert.xyz[iv] = mid
	return

def close_two_vertices(v1,v2,vert,ratio):
	# make two vertices closer
	xyz1 = vert.xyz[v1]
	xyz2 = vert.xyz[v2]
	mid = line_split(xyz1,xyz2,0.5)
	
	mid1 = line_split(mid,xyz1,ratio)
	mid2 = line_split(mid,xyz2,ratio)

	vert.xyz[v1] = mid1
	vert.xyz[v2] = mid2
	
	for iv in range(0,vert.nv):
		if iv != v1:
			dist = distance(xyz1,vert.xyz[iv])
			if dist < 1e-5: vert.xyz[iv] = mid1
	
	for iv in range(0,vert.nv):
		if iv != v2:
			dist = distance(xyz2,vert.xyz[iv])
			if dist < 1e-5: vert.xyz[iv] = mid2
			
	return

	
def angles_three_vertices(xyz1,xyz2,xyz3):
	# get angle v2-v1-v3

	vec12 = [xyz2[0]-xyz1[0],xyz2[1]-xyz1[1],xyz2[2]-xyz1[2]]
	vec13 = [xyz3[0]-xyz1[0],xyz3[1]-xyz1[1],xyz3[2]-xyz1[2]]
	
	angle = np.dot(vec12,vec13)/(np.linalg.norm(vec12)*np.linalg.norm(vec13))

	if angle > 1.0: angle = 1.0
	if angle < -1.0: angle = -1.0
	
	angle = math.acos(angle)*180.0/math.pi
	
	nv = np.cross(vec12,vec13)
	norm = np.linalg.norm(nv)

	nv = [nv[0]/norm,nv[1]/norm,nv[2]/norm]
	
	return angle,nv

def angles_two_vectors(vec1,vec2):
	
	angle = np.dot(vec1,vec2)/(np.linalg.norm(vec1)*np.linalg.norm(vec2))

	if angle > 1.0: angle = 1.0
	if angle < -1.0: angle = -1.0
	
	angle = math.acos(angle)*180.0/math.pi
	
	return angle

def triangle_area(xyz1,xyz2,xyz3):
	# get angle v2-v1-v3

	vec12 = [xyz2[0]-xyz1[0],xyz2[1]-xyz1[1],xyz2[2]-xyz1[2]]
	vec13 = [xyz3[0]-xyz1[0],xyz3[1]-xyz1[1],xyz3[2]-xyz1[2]]
	
	nv = np.cross(vec12,vec13)
	area = 0.5* np.linalg.norm(nv)
	
	return area

def vector_minus(A,B):
	C = [A[0]-B[0],A[1]-B[1],A[2]-B[2]]
	return C
	
def VOLUM0(A,B,C,D): 
	# same subroutine in nek
	vecDA = vector_minus(A,D)
	vecDB = vector_minus(B,D)
	devDC = vector_minus(C,D)
	
	#vec1 = np.cross(vecDA,vecDB)
	#vdot = np.dot(vec1,devDC)
	
	vec1 = cross_vector(vecDA,vecDB)
	vdot = dot_vector(vec1,devDC)
	
	return vdot
	
def cross_vector(v1,v2):
	res = [0,0,0]
	
	res[0] = v1[1]*v2[2] -  v1[2]*v2[1]
	res[1] = v1[2]*v2[0] -  v1[0]*v2[2]
	res[2] = v1[0]*v2[1] -  v1[1]*v2[0]
	return res
	
def dot_vector(v1,v2):
	d = v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2]
	return d

def remove_duplicates(alist):
	# remove duplicates in the list, but remain the original order
	res = []
	for i in alist:
		if i not in res:
			res.append(i)
	return res
	

def remove_duplicates_cqv(cqv):
	# remove duplicates in the list, but remain the original order
	res = []
	for c1 in cqv:
		if len(res)>0:
			ifadd = True
			for c2 in res:
				if if_a2_equal(c1,c2):
					ifadd = False
			if ifadd:
				res.append(c1)
		else:
			res.append(c1)
	return res

def if_a2_equal(a,b):
	if(a[0]==b[0]) and (a[1]==b[1]):
		return True
	else:
		return False
