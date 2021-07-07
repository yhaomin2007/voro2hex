from data_class import *
from splitting_subroutines import * 
import math

def polygon_divide_for_chamfer_polygon(iface,face,vert,pebbles,pdiameter):
	nvert = len(face.f_to_v[iface])
	if (not face.ifghost[iface]) and (not face.ifcollapse[iface]) and (face.ifchamfer[iface]):
		lengths = face.lengths[iface]
		minlength = min(lengths)
		minlengthIndex = lengths.index(minlength)
		plist = face.f_to_p[iface]
		p1 = plist[0]
		p1xyz = pebbles.xyz[p1]
		p2 = plist[1]
		p2xyz = pebbles.xyz[p2]
		pmiddle = line_split(p1xyz,p2xyz,0.5)
		
		ftov = face.f_to_v[iface]
		
		newftov = []
		for iv in range(0,nvert):
			newftov.append(ftov[(minlengthIndex+iv)%nvert])
			
		if minlength < 0.2*pdiameter:
			# shrink into smaller polygon,
			face.ifcollapse[iface] = True
			
			# generate 
			
	

def polygon_divide_for_octagon_new(iface,face,vert):
	nvert = len(face.f_to_v[iface])
	if (not face.ifghost[iface]) and (not face.ifcollapse[iface]) and (not face.ifchamfer[iface]) and (nvert == 8):
		# pack all vertices
		# pack all edges
		# find smallest edge length and start with it.
		ftov = face.f_to_v[iface]
		plist = face.f_to_p[iface]
		lengths = face.lengths[iface]
		minlength = min(lengths)
		minlengthIndex = lengths.index(minlength)
		
		startv = minlengthIndex
		
		newftov = []
		for iv in range(0,nvert):
			newftov.append(ftov[(minlengthIndex+iv)%nvert])
		
		# generate mid points (4)
	
		mv1 = line_split(vert.xyz[newftov[0]],vert.xyz[newftov[1]],0.5)
		mv2 = line_split(vert.xyz[newftov[2]],vert.xyz[newftov[3]],0.5)
		mv3 = line_split(vert.xyz[newftov[4]],vert.xyz[newftov[5]],0.5)
		mv4 = line_split(vert.xyz[newftov[6]],vert.xyz[newftov[7]],0.5)
		
		fc = face.f_center[iface]
		m1 = line_split(mv1,fc,0.4)
		m2 = line_split(mv2,fc,0.4)
		m3 = line_split(mv3,fc,0.4)
		m4 = line_split(mv4,fc,0.4)
		
		vert.new_vertice(m1[0],m1[1],m1[2])
		m1v = vert.nv -1
		vert.new_vertice(m2[0],m2[1],m2[2])
		m2v = vert.nv -1
		vert.new_vertice(m3[0],m3[1],m3[2])
		m3v = vert.nv -1
		vert.new_vertice(m4[0],m4[1],m4[2])
		m4v = vert.nv -1
		
		
		# generate new polygon (quads and triangles)
		face.ifcollapse[iface] = True
		
		# triangles
		new_vlist = [m1v,newftov[0],newftov[1]]
		face.new_face(new_vlist,plist)
		
		new_vlist = [m2v,newftov[2],newftov[3]]
		face.new_face(new_vlist,plist)
		
		new_vlist = [m3v,newftov[4],newftov[5]]
		face.new_face(new_vlist,plist)
		
		new_vlist = [m4v,newftov[6],newftov[7]]
		face.new_face(new_vlist,plist)
		
		# quads
		new_vlist = [m1v,newftov[1],newftov[2],m2v]
		face.new_face(new_vlist,plist)
		
		new_vlist = [m2v,newftov[3],newftov[4],m3v]
		face.new_face(new_vlist,plist)
		
		new_vlist = [m3v,newftov[5],newftov[6],m4v]
		face.new_face(new_vlist,plist)
		
		new_vlist = [m4v,newftov[7],newftov[0],m1v]
		face.new_face(new_vlist,plist)
		
		#center quad
		new_vlist = [m1v,m2v,m3v,m4v]
		face.new_face(new_vlist,plist)
	
	
def polygon_divide_for_heptagon_new(iface,face,vert):
	nvert = len(face.f_to_v[iface])
	if (not face.ifghost[iface]) and (not face.ifcollapse[iface]) and (not face.ifchamfer[iface]) and (nvert == 7):
		# pack all vertices
		# pack all edges
		# find smallest edge length and start with it.
		ftov = face.f_to_v[iface]
		plist = face.f_to_p[iface]
		lengths = face.lengths[iface]
		minlength = min(lengths)
		minlengthIndex = lengths.index(minlength)
		
		startv = minlengthIndex
		
		newftov = []
		for iv in range(0,nvert):
			newftov.append(ftov[(minlengthIndex+iv)%nvert])

		# generate mid points (4)
		mv1 = line_split(vert.xyz[newftov[0]],vert.xyz[newftov[1]],0.5)
		mv2 = line_split(vert.xyz[newftov[2]],vert.xyz[newftov[3]],0.5)
		mv3 = vert.xyz[newftov[4]]
		mv4 = line_split(vert.xyz[newftov[5]],vert.xyz[newftov[6]],0.5)
		
		fc = face.f_center[iface]
		m1 = line_split(mv1,fc,0.4)
		m2 = line_split(mv2,fc,0.4)
		m3 = line_split(mv3,fc,0.4)
		m4 = line_split(mv4,fc,0.4)
		
		vert.new_vertice(m1[0],m1[1],m1[2])
		m1v = vert.nv -1
		vert.new_vertice(m2[0],m2[1],m2[2])
		m2v = vert.nv -1
		vert.new_vertice(m3[0],m3[1],m3[2])
		m3v = vert.nv -1
		vert.new_vertice(m4[0],m4[1],m4[2])
		m4v = vert.nv -1
		
		
		# generate new polygon (quads and triangles)
		face.ifcollapse[iface] = True
		
		# triangles
		new_vlist = [m1v,newftov[0],newftov[1]]
		face.new_face(new_vlist,plist)
		
		new_vlist = [m2v,newftov[2],newftov[3]]
		face.new_face(new_vlist,plist)
		
		#new_vlist = [m3v,newftov[4],newftov[5]]
		#face.new_face(new_vlist,plist)
		
		new_vlist = [m4v,newftov[5],newftov[6]]
		face.new_face(new_vlist,plist)
		
		# quads
		new_vlist = [m1v,newftov[1],newftov[2],m2v]
		face.new_face(new_vlist,plist)
		
		new_vlist = [m2v,newftov[3],newftov[4],m3v]
		face.new_face(new_vlist,plist)
		
		new_vlist = [m3v,newftov[4],newftov[5],m4v]
		face.new_face(new_vlist,plist)
		
		new_vlist = [m4v,newftov[6],newftov[0],m1v]
		face.new_face(new_vlist,plist)
		
		#center quad
		new_vlist = [m1v,m2v,m3v,m4v]
		face.new_face(new_vlist,plist)
	
def polygon_divide_for_hexagon_new(iface,face,vert):

	nvert = len(face.f_to_v[iface])
	if (not face.ifghost[iface]) and (not face.ifcollapse[iface]) and (not face.ifchamfer[iface]) and (nvert == 6):
		# pack all vertices
		# pack all edges
		# find smallest edge length and start with it.
		ftov = face.f_to_v[iface]
		plist = face.f_to_p[iface]
		lengths = face.lengths[iface]
		minlength = min(lengths)
		minlengthIndex = lengths.index(minlength)
		
		startv = minlengthIndex
		
		newftov = []
		for iv in range(0,nvert):
			newftov.append(ftov[(minlengthIndex+iv)%nvert])
		
		# generate mid points (4)
	
		mv1 = line_split(vert.xyz[newftov[0]],vert.xyz[newftov[1]],0.5)
		mv2 = vert.xyz[newftov[2]]
		mv3 = line_split(vert.xyz[newftov[3]],vert.xyz[newftov[4]],0.5)
		mv4 = vert.xyz[newftov[5]]
		
		fc = face.f_center[iface]
		m1 = line_split(mv1,fc,0.4)
		m2 = line_split(mv2,fc,0.4)
		m3 = line_split(mv3,fc,0.4)
		m4 = line_split(mv4,fc,0.4)
		
		vert.new_vertice(m1[0],m1[1],m1[2])
		m1v = vert.nv -1
		vert.new_vertice(m2[0],m2[1],m2[2])
		m2v = vert.nv -1
		vert.new_vertice(m3[0],m3[1],m3[2])
		m3v = vert.nv -1
		vert.new_vertice(m4[0],m4[1],m4[2])
		m4v = vert.nv -1
		
		
		# generate new polygon (quads and triangles)
		face.ifcollapse[iface] = True
		
		# triangles
		new_vlist = [m1v,newftov[0],newftov[1]]
		face.new_face(new_vlist,plist)
		
		#new_vlist = [m2v,newftov[2],newftov[3]]
		#face.new_face(new_vlist,plist)
		
		new_vlist = [m3v,newftov[3],newftov[4]]
		face.new_face(new_vlist,plist)
		
		#new_vlist = [m4v,newftov[5],newftov[6]]
		#face.new_face(new_vlist,plist)
		
		# quads
		new_vlist = [m1v,newftov[1],newftov[2],m2v]
		face.new_face(new_vlist,plist)
		
		new_vlist = [m2v,newftov[2],newftov[3],m3v]
		face.new_face(new_vlist,plist)
		
		new_vlist = [m3v,newftov[4],newftov[5],m4v]
		face.new_face(new_vlist,plist)
		
		new_vlist = [m4v,newftov[5],newftov[0],m1v]
		face.new_face(new_vlist,plist)
		
		#center quad
		new_vlist = [m1v,m2v,m3v,m4v]
		face.new_face(new_vlist,plist)


def polygon_divide_for_pentagon_new(iface,face,vert):
	nvert = len(face.f_to_v[iface])
	if (not face.ifghost[iface]) and (not face.ifcollapse[iface]) and (not face.ifchamfer[iface]) and (nvert == 5):
		# pack all vertices
		# pack all edges
		# find smallest edge length and start with it.
		ftov = face.f_to_v[iface]
		plist = face.f_to_p[iface]
		lengths = face.lengths[iface]
		minlength = min(lengths)
		minlengthIndex = lengths.index(minlength)
		
		startv = minlengthIndex
		
		newftov = []
		for iv in range(0,nvert):
			newftov.append(ftov[(minlengthIndex+iv)%nvert])
		
		# generate mid points (4)
	
		mv1 = line_split(vert.xyz[newftov[0]],vert.xyz[newftov[1]],0.5)
		mv2 =  vert.xyz[newftov[2]]
		mv3 =  vert.xyz[newftov[3]]
		mv4 =  vert.xyz[newftov[4]]
		
		fc = face.f_center[iface]
		m1 = line_split(mv1,fc,0.4)
		m2 = line_split(mv2,fc,0.4)
		m3 = line_split(mv3,fc,0.4)
		m4 = line_split(mv4,fc,0.4)
		
		vert.new_vertice(m1[0],m1[1],m1[2])
		m1v = vert.nv -1
		vert.new_vertice(m2[0],m2[1],m2[2])
		m2v = vert.nv -1
		vert.new_vertice(m3[0],m3[1],m3[2])
		m3v = vert.nv -1
		vert.new_vertice(m4[0],m4[1],m4[2])
		m4v = vert.nv -1
		
		# generate new polygon (quads and triangles)
		face.ifcollapse[iface] = True
		
		# triangles
		new_vlist = [m1v,newftov[0],newftov[1]]
		face.new_face(new_vlist,plist)
		
		# quads
		new_vlist = [m1v,newftov[1],newftov[2],m2v]
		face.new_face(new_vlist,plist)
		
		new_vlist = [m2v,newftov[2],newftov[3],m3v]
		face.new_face(new_vlist,plist)
		
		new_vlist = [m3v,newftov[3],newftov[4],m4v]
		face.new_face(new_vlist,plist)
		
		new_vlist = [m4v,newftov[4],newftov[0],m1v]
		face.new_face(new_vlist,plist)
		
		#center quad
		new_vlist = [m1v,m2v,m3v,m4v]
		face.new_face(new_vlist,plist)


def polygon_divide_for_heptagon(iface,face,vert):

	nvert = len(face.f_to_v[iface])
	if (not face.ifghost[iface]) and (not face.ifcollapse[iface]) and (not face.ifchamfer[iface]) and (nvert == 7):
		maxAngle = max(face.angles[iface])
		maxAngleIndex = face.angles[iface].index(maxAngle)
		v1 = maxAngleIndex
		
		face.ifcollapse[iface] = True
		vlist = face.f_to_v[iface]
		plist = face.f_to_p[iface]
		fc = face.f_center[iface]
		vert.new_vertice(fc[0],fc[1],fc[2])
		
		fcv = vert.nv -1
		
		for isplit in range(0,3):
			new_vlist = [fcv]
			for iv in range(v1,v1+3):
				vr = iv%(nvert)
				vr = vlist[vr]
				new_vlist.append(vr)
			v1 = v1+2
			face.new_face(new_vlist,plist)
			
		new_vlist = [fcv]
		for iv in range(v1,v1+2):
			vr = iv%(nvert)
			vr = vlist[vr]
			new_vlist.append(vr)
		face.new_face(new_vlist,plist)

def polygon_divide_for_pentagon_tri_split(iface,face,vert):
	#divide polygon for pentagon.
	# max angle is larger than 120, but second max angles is smaller then 100
	
	nvert = len(face.f_to_v[iface])
	if (not face.ifghost[iface]) and (not face.ifcollapse[iface]) and (not face.ifchamfer[iface]) and (nvert == 5):
		angles = face.angles[iface]
		maxAngle = max(angles)
		maxAngleIndex = angles.index(maxAngle)
		
		angles.remove(maxAngle)
		maxAngle2 = max(angles)
		
		angles = face.angles[iface]
		
		if maxAngle > 130 and maxAngle2 < 110:
			
			v1 = maxAngleIndex
		
			face.ifcollapse[iface] = True
			vlist = face.f_to_v[iface]
			plist = face.f_to_p[iface]
			
			fc = face.f_center[iface]
			vert.new_vertice(fc[0],fc[1],fc[2])
			fcv = vert.nv -1			

			new_vlist = [fcv]
			for iv in range(v1,v1+3):
				vr = iv%(nvert)
				vr = vlist[vr]
				new_vlist.append(vr)
			face.new_face(new_vlist,plist)
			
			new_vlist = [fcv]
			for iv in range(v1+2,v1+4):
				vr = iv%(nvert)
				vr = vlist[vr]
				new_vlist.append(vr)
			face.new_face(new_vlist,plist)
			
			new_vlist = [fcv]
			for iv in range(v1+3,v1+6):
				vr = iv%(nvert)
				vr = vlist[vr]
				new_vlist.append(vr)
			new_vlist.append(vlist[v1])
			face.new_face(new_vlist,plist)

def polygon_divide_for_pentagon_tri_split2(iface,face,vert):
	#divide polygon for pentagon.
	
	nvert = len(face.f_to_v[iface])
	if (not face.ifghost[iface]) and (not face.ifcollapse[iface]) and (not face.ifchamfer[iface]) and (nvert == 5):
		angles = face.angles[iface]
		maxAngle = max(angles)
		maxAngleIndex = angles.index(maxAngle)
		
		angles.remove(maxAngle)
		maxAngle2 = max(angles)
		
		angles = face.angles[iface]
		
		maxAngleIndex2 = angles.index(maxAngle2)
		
		v1 = maxAngleIndex
		v2 = maxAngleIndex2
		ifnbr = False
		if v2 == (v1+1)%nvert: 
			ifnbr= True
		if v2 == (v1-1)%nvert: 
			ifnbr= True
			v3 = v1
			v1 = v2
			v2 = v3
		
		if maxAngle > 120 and maxAngle2 > 120 and ifnbr:
			face.ifcollapse[iface] = True
			vlist = face.f_to_v[iface]
			plist = face.f_to_p[iface]
			
			fc = face.f_center[iface]
			vert.new_vertice(fc[0],fc[1],fc[2])
			fcv = vert.nv -1	
			
			new_vlist = [fcv]
			for iv in range(v1,v1+2):
				vr = iv%(nvert)
				vr = vlist[vr]
				new_vlist.append(vr)
			face.new_face(new_vlist,plist)
			
			new_vlist = [fcv]
			for iv in range(v1+1,v1+4):
				vr = iv%(nvert)
				vr = vlist[vr]
				new_vlist.append(vr)
			face.new_face(new_vlist,plist)
			
			new_vlist = [fcv]
			for iv in range(v1+3,v1+6):
				vr = iv%(nvert)
				vr = vlist[vr]
				new_vlist.append(vr)
			new_vlist.append(vlist[v1])
			face.new_face(new_vlist,plist)
			
def polygon_divide_for_pentagon_bi_split(iface,face,vert):
	#divide polygon for pentagon.
	
	nvert = len(face.f_to_v[iface])
	if (not face.ifghost[iface]) and (not face.ifcollapse[iface]) and (not face.ifchamfer[iface]) and (nvert == 5):
		angles = face.angles[iface]
		maxAngle = max(angles)
		maxAngleIndex = angles.index(maxAngle)
		
		angles.remove(maxAngle)
		maxAngle2 = max(angles)
		
		angles = face.angles[iface]
		
		maxAngleIndex2 = angles.index(maxAngle2)
		
		v1 = maxAngleIndex
		v2 = maxAngleIndex2
		ifnnbr = True
		if v2 == (v1+1)%nvert: ifnnbr= False
		if v2 == (v1-1)%nvert: ifnnbr= False
		
		if maxAngle > 120 and maxAngle2 > 120 and ifnnbr:

			face.ifcollapse[iface] = True
			vlist = face.f_to_v[iface]
			plist = face.f_to_p[iface]
			
			fc = face.f_center[iface]
			vert.new_vertice(fc[0],fc[1],fc[2])
			fcv = vert.nv -1		
		
		
			if v2 > v1: 
				di = v2-v1
			else:
				v3 = v1
				v1 = v2
				v2 = v3
				di = v2-v1
			
			new_vlist = []
			for iv in range(v1,v1+di+1):
				vr = iv%(nvert)
				vr = vlist[vr]
				new_vlist.append(vr)
			face.new_face(new_vlist,plist)
			
			new_vlist = []
			for iv in range(v1+di,v1+6):
				vr = iv%(nvert)
				vr = vlist[vr]
				new_vlist.append(vr)
			face.new_face(new_vlist,plist)
	
			
def polygon_divide_for_hexagon(iface,face,vert):
	#divide polygon for hexagon
	
	nvert = len(face.f_to_v[iface])
	if (not face.ifghost[iface]) and (not face.ifcollapse[iface]) and (not face.ifchamfer[iface]) and (nvert == 6):
		maxAngle = max(face.angles[iface])
		maxAngleIndex = face.angles[iface].index(maxAngle)
		v1 = maxAngleIndex
		
		face.ifcollapse[iface] = True
		vlist = face.f_to_v[iface]
		plist = face.f_to_p[iface]
			
		new_vlist = []
		for iv in range(v1,v1+3+1):
			vr = iv%(nvert)
			vr = vlist[vr]
			new_vlist.append(vr)
		face.new_face(new_vlist,plist)
		
		new_vlist = []
		for iv in range(v1+3,v1+3+3+1):
			vr = iv%(nvert)
			vr = vlist[vr]
			new_vlist.append(vr)
		face.new_face(new_vlist,plist)
	

def polygon_divide_new(iface,face,vert,maxAngle):
	#divide polygon into smaller polygons
	#by linking mid face point and all vertices with large angle
	
	nvert = len(face.f_to_v[iface])
	if (not face.ifghost[iface]) and (not face.ifcollapse[iface]) and (not face.ifchamfer[iface]) and (nvert > 4):
		
		nsplit = 0
		split_vertices = []
		iv = 0
		for angle in face.angles[iface]:
			if angle > maxAngle:
				nsplit = nsplit + 1
				split_vertices.append(iv)
			iv = iv + 1
			
		if nsplit == 0: 
			return
		elif nsplit == 1:
			# find the vertice in the oppsite side to split this polygon
			face.ifcollapse[iface] = True
			vlist = face.f_to_v[iface]
			plist = face.f_to_p[iface]
			
			
			v1 = split_vertices[0]
			v2 = (v1 - 1)%(nvert)
			v3 = (v1 + 1)%(nvert)
			angle_v1 = face.angles[iface][v1]
			
			xyz1 = vert.xyz[face.f_to_v[iface][v1]]
			xyz2 = vert.xyz[face.f_to_v[iface][v2]]
			
			
			min_include_angle = 0.5*angle_v1
			# find v_opposite
			for iv in range(0,nvert-3):
				ivo = (iv+v3+1)%(nvert)
				vro =  vlist[ivo]
				xyzo = vert.xyz[vro]
				
				angle_o,nv = angles_three_vertices(xyz1,xyz2,xyzo)
				include_angle = abs(angle_o - 0.5*angle_v1)
				if (include_angle<min_include_angle):
					min_include_angle = include_angle
					v_opposite = ivo
			
			v2 = v_opposite
			
			if v2 > v1:
				None
			else:
				# swap v1,v2
				v3 = v2
				v2 = v1
				v1 = v3
			
			n12 = v2-v1
			n21 = nvert - n12
			# split between v1 and v_opposite
						
			new_vlist = []
			for iv in range(v1,v1+n12+1):
				vr = iv%(nvert)
				vr = vlist[vr]
				new_vlist.append(vr)
			face.new_face(new_vlist,plist)
			
			new_vlist = []
			for iv in range(v2,v2+n21+1):
				vr = iv%(nvert)
				vr = vlist[vr]
				new_vlist.append(vr)
			face.new_face(new_vlist,plist)
			
			
		elif nsplit == 2:
			# split through these two vertice
			face.ifcollapse[iface] = True
			vlist = face.f_to_v[iface]
			plist = face.f_to_p[iface]
			
			v1 = split_vertices[0]
			v2 = split_vertices[1]
			
			new_vlist = []
			for iv in range(v1,v2+1):
				vr = iv%(nvert)
				vr = vlist[vr]
				new_vlist.append(vr)
			face.new_face(new_vlist,plist)
			
			new_vlist = []
			for iv in range(v2,v1+nvert+1):
				vr = iv%(nvert)
				vr = vlist[vr]
				new_vlist.append(vr)
			face.new_face(new_vlist,plist)
			
		else:
		
			face.ifcollapse[iface] = True
			vlist = face.f_to_v[iface]
			plist = face.f_to_p[iface]
			fc = face.f_center[iface]
			vert.new_vertice(fc[0],fc[1],fc[2])
			fcv = vert.nv -1
			
			for isplit in range(0,nsplit-1):
				v1 = split_vertices[isplit]
				v2 = split_vertices[isplit+1]
				
				new_vlist = [fcv]
				for iv in range(v1,v2+1):
					vr = iv%(nvert)
					vr = vlist[vr]
					new_vlist.append(vr)
				face.new_face(new_vlist,plist)
			
			v1 = split_vertices[nsplit-1]
			v2 = split_vertices[0]
			
			n12 = nvert - (v1 - v2)
			
			# for last splitted polygon
			new_vlist = [fcv]
			for iv in range(v1,v1+n12+1):
				vr = iv%(nvert)
				vr = vlist[vr]
				new_vlist.append(vr)
			face.new_face(new_vlist,plist)
			

def polygon_divide(iface,face,vert):
	#divide polygon into smaller polygons
	#add new vertices if needed
	# then turn off this face
	nvert = len(face.f_to_v[iface])
	if (not face.ifghost[iface]) and (not face.ifcollapse[iface]) and (not face.ifchamfer[iface]) and nvert > 4:
		maxAngle = max(face.angles[iface])
		maxAngleIndex = face.angles[iface].index(maxAngle)
		v1 = maxAngleIndex
		v1_in_vert = face.f_to_v[iface][v1]
		xyz1 = vert.xyz[v1_in_vert]
		fc = face.f_center[iface]
		
		## find v2,v3 to divide this polygon.
		nv_prev = None
		v2 = 2
		if nvert == 5: 
			v3 = 3
		else:
			v3 = v2 + 2
		
			v2_in_vert = face.f_to_v[iface][(v2+v1)%(nvert-1)]
			xyz2 = vert.xyz[v2_in_vert]
			angle_12,nv_vr = angles_three_vertices(fc,xyz2,xyz1)
		
			if angle_12 < 60.0: 
				v2 = v2 + 1
				v3 = v3 + 1
		
		v2_in_vert = face.f_to_v[iface][(v2+v1)%(nvert-1)]
		xyz2 = vert.xyz[v2_in_vert]
		v3_in_vert = face.f_to_v[iface][(v3+v1)%(nvert-1)]
		xyz3 = vert.xyz[v3_in_vert]
		angle_23,nv_vr = angles_three_vertices(fc,xyz2,xyz3)
		if angle_23 < 40.0:
			v3 = v3 + 1
		
			
		face.ifcollapse[iface] = True # collapes this face
		
		vert.new_vertice(fc[0],fc[1],fc[2])
		fcv = vert.nv -1
		
		vlist = face.f_to_v[iface]
		plist = face.f_to_p[iface]
		
		#print v1,v2,v3,nvert
		#print vlist
		
		new_vlist = [fcv]
		for iv in range(v1,v1+v2+1):
			vr = iv
			if(iv > nvert-1): vr = iv - nvert
			vr = vlist[vr]
			new_vlist.append(vr)
		face.new_face(new_vlist,plist)
		
		new_vlist = [fcv]
		for iv in range(v1+v2,v1+v3+1):
			vr = iv
			if(iv > nvert-1): vr = iv - nvert
			vr = vlist[vr]
			new_vlist.append(vr)
		face.new_face(new_vlist,plist)

		new_vlist = [fcv]
		for iv in range(v1+v3,v1+nvert+1):
			vr = iv
			if(iv > nvert-1): vr = iv - nvert
			vr = vlist[vr]
			new_vlist.append(vr)
		face.new_face(new_vlist,plist)

def concave_quad_fix1(iface,face,vert):
	# fix concave quad by morphing
	nvert = len(face.f_to_v[iface])
	if (not face.ifghost[iface]) and (not face.ifcollapse[iface]) and (not face.ifchamfer[iface]) and nvert == 4:
		maxAngle = max(face.angles[iface])
		if maxAngle > 170:
			maxAngleIndex = face.angles[iface].index(maxAngle)
			v1 = maxAngleIndex
	
			# pull v2 and v4
			# push v1 
			# repsective to v3
			v2 = v1 + 1
			if v2 > 3: v2 = v2 -4
			v3 = v1 + 2
			if v3 > 3: v3 = v3 -4
			v4 = v1 + 3
			if v4 > 3: v4 = v4 -4
			
			v1 = face.f_to_v[iface][v1]
			v2 = face.f_to_v[iface][v2]
			v3 = face.f_to_v[iface][v3]
			v4 = face.f_to_v[iface][v4]
			
			xyz1 =	vert.xyz[v1]
			xyz2 =	vert.xyz[v2]
			xyz3 =	vert.xyz[v3]
			xyz4 =	vert.xyz[v4]
		
			new1 = line_split(xyz3,xyz1,1.1)
			new2 = line_split(xyz3,xyz2,0.95)
			new4 = line_split(xyz3,xyz4,0.95)
			
			set_vertices(v1,vert,new1)
			set_vertices(v2,vert,new2)
			set_vertices(v4,vert,new4)
	
def concave_quad_fix2(iface,face,vert,angle):
	# fix concave quad by splitting
	nvert = len(face.f_to_v[iface])
	if (not face.ifghost[iface]) and (not face.ifcollapse[iface]) and (not face.ifchamfer[iface]) and nvert == 4:
		maxAngle = max(face.angles[iface])
		if maxAngle > angle:
			maxAngleIndex = face.angles[iface].index(maxAngle)
			v1 = maxAngleIndex

			v2 = v1 + 1
			if v2 > 3: v2 = v2 -4
			v3 = v1 + 2
			if v3 > 3: v3 = v3 -4
			v4 = v1 + 3
			if v4 > 3: v4 = v4 -4
			
			v1 = face.f_to_v[iface][v1]
			v2 = face.f_to_v[iface][v2]
			v3 = face.f_to_v[iface][v3]
			v4 = face.f_to_v[iface][v4]
			
			
			face.ifcollapse[iface] = True

			plist = face.f_to_p[iface]
			new_vlist = [v1,v2,v3]
			face.new_face(new_vlist,plist)
			new_vlist = [v3,v4,v1]
			face.new_face(new_vlist,plist)
	
def concave_quad_fix3(iface,face,vert,angle):
	# fix concave quad by merging
	nvert = len(face.f_to_v[iface])
	if (not face.ifghost[iface]) and (not face.ifcollapse[iface]) and (not face.ifchamfer[iface]) and nvert == 4:
		maxAngle = max(face.angles[iface])
		if maxAngle > angle:
			maxAngleIndex = face.angles[iface].index(maxAngle)
			v1 = maxAngleIndex
	
			# pull v2 and v4
			# push v1 
			# repsective to v3
			v2 = v1 + 1
			if v2 > 3: v2 = v2 -4
			v3 = v1 + 2
			if v3 > 3: v3 = v3 -4
			v4 = v1 + 3
			if v4 > 3: v4 = v4 -4
			
			v1 = face.f_to_v[iface][v1]
			v2 = face.f_to_v[iface][v2]
			v3 = face.f_to_v[iface][v3]
			v4 = face.f_to_v[iface][v4]
			
			v1xyz = vert.xyz[v1]
			v2xyz = vert.xyz[v2]
			v3xyz = vert.xyz[v3]

			d12 = distance(v1xyz,v2xyz)
			d13 = distance(v1xyz,v3xyz)
			
			face.f_to_v[iface].remove(v1)
			
			if d12 <= d13:
				# merge v1 to v2
				vert.mgd_v[v2].extend(vert.mgd_v[v1])
				vert.mgd_v[v2].append(v1)
				vert.mgd_v[v2]= remove_duplicates(vert.mgd_v[v2])
				for mgd in vert.mgd_v[v2]:
					vert.mgd_v[mgd] = vert.mgd_v[v2]
				set_vertices(v2,vert,v2xyz)
			else:
				# merge v1 to v3
				vert.mgd_v[v3].extend(vert.mgd_v[v1])
				vert.mgd_v[v3].append(v1)
				vert.mgd_v[v3]= remove_duplicates(vert.mgd_v[v3])
				for mgd in vert.mgd_v[v3]:
					vert.mgd_v[mgd] = vert.mgd_v[v3]
				set_vertices(v3,vert,v3xyz)


def fix_concave_top_quads(top_quads,pebbles):
	nquads = top_quads.nquads
	tot_fixed_quads = 0
	for iquad in range(0,nquads):
		# evaluate four angles of this quad
		angles = []
		verts = top_quads.xyz[iquad]
		
		center1 = line_split(verts[0],verts[1],0.5)
		center2 = line_split(verts[2],verts[3],0.5)
		center = line_split(center1,center2,0.5)
		
		# corner 0
		angle1,nv = angles_three_vertices(verts[0],verts[1],center)
		angle2,nv = angles_three_vertices(verts[0],verts[3],center)
		angle = angle1 + angle2
		angles.append(angle)
		
		# corner 1
		angle1,nv = angles_three_vertices(verts[1],verts[2],center)
		angle2,nv = angles_three_vertices(verts[1],verts[0],center)
		angle = angle1 + angle2
		angles.append(angle)
		
		# corner 2
		angle1,nv = angles_three_vertices(verts[2],verts[1],center)
		angle2,nv = angles_three_vertices(verts[2],verts[3],center)
		angle = angle1 + angle2
		angles.append(angle)
	
		# corner 3
		angle1,nv = angles_three_vertices(verts[3],verts[2],center)
		angle2,nv = angles_three_vertices(verts[3],verts[0],center)
		angle = angle1 + angle2
		angles.append(angle)
		
		# if the max angle is larger than 165, than move this corner
		if max(angles) > 160:
			#print 'one concave top-quad fixed'
			tot_fixed_quads = tot_fixed_quads +1 
			ipebble = top_quads.quad_to_p[iquad][0]
			pxyz = pebbles.xyz[ipebble]
			
			moveCorner = angles.index(max(angles))
			oppositeCorner = (moveCorner+2)%4
			
			oldCornerXYZ = verts[moveCorner]
			
			#oldCornerXYZ1 = line_split(verts[(moveCorner+1)%4],verts[(moveCorner-1)%4],0.5)
			
			oppositeCornerXYZ = verts[oppositeCorner]
			
			newCornerXYZ = line_split(oppositeCornerXYZ,oldCornerXYZ,1.035)
			
			# need to project this newCornerXYZ to pebble sueface
			dist = distance(pxyz,newCornerXYZ)
			proj_radius  =distance(pxyz,oldCornerXYZ)
			
			newCornerXYZ = line_split(pxyz,newCornerXYZ,proj_radius/dist)
			
			top_quads.xyz[iquad][moveCorner] = newCornerXYZ
			
			top_quads.linear_mid(iquad)
				
			same_pebble_top_quads_list = pebbles.top_quads[ipebble]
			
			for iquad2 in same_pebble_top_quads_list:
				if iquad2 != iquad:
					verts2 = top_quads.xyz[iquad2]
					for i in range(0,4):
						dist = distance(oldCornerXYZ,verts2[i])
						if dist < 1e-4:
							top_quads.xyz[iquad2][i] =  newCornerXYZ
							top_quads.linear_mid(iquad2)
		
	print 'tot_fixed_quads: ' + str(tot_fixed_quads) 
		
		
	
