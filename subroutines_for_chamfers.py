from splitting_subroutines import *

def shrink_towards_center(iface,face,vert,vert_reduction):

	quads = []
	nverts = len(face.f_to_v[iface])
	fc = face.f_center[iface]
	p1 = face.f_to_p[iface][0]
	p2 = face.f_to_p[iface][1]
	plist = [p1,p2]
	
	if nverts > 8:
		# 1st, create inner polygon
		vlist = []
		for iv in range(0,nverts):
			v1 = face.f_to_v[iface][iv]
			v1xyz = vert.xyz[v1]
			v2xyz = line_split(v1xyz,fc,0.2) # inner polygon vertice
			vert.new_vertice(v2xyz[0],v2xyz[1],v2xyz[2])
			vlist.append(vert.nv-1)

		face.new_face(vlist,plist)
		iface_new = face.nf-1

		face.find_edge_info(iface_new)
		face.calculate_edge_length(iface_new,vert)

		#2nd, collapse small edge of inner polygon
		# collapse the small vert_reduction edges
		lengths = face.lengths[iface_new]
		lengths.sort()
		threshold_length = lengths[vert_reduction]
		# now collapse edges that < threshold_length
		collapsed = []
		
		for iv in range(0,nverts):
			length = face.lengths[iface_new][iv]
			collapsed.append(False)
			if length < threshold_length:
				if iv > 0 and  collapsed[iv-1]== True :
					# jump this one, and increase threshold_length
					vert_reduction = vert_reduction + 1
					threshold_length = lengths[vert_reduction]
				else:
				
					collapsed[iv] = True
				
					iv1 = iv
					iv2 = (iv +1)%nverts
				
					iv1 = face.f_to_v[iface_new][iv1]
					iv2 = face.f_to_v[iface_new][iv2]
				
					v1xyz = vert.xyz[iv1]
					v2xyz = vert.xyz[iv2]
				
					midxyz = line_split(v1xyz,v2xyz,0.5)
				
					vert.xyz[iv1] = midxyz
					vert.xyz[iv2] = midxyz
				
					vert.mgd_v[iv1] = vert.mgd_v[iv1] + vert.mgd_v[iv2]
					vert.mgd_v[iv1]= remove_duplicates(vert.mgd_v[iv1])
					vert.mgd_v[iv2] = vert.mgd_v[iv1]
				
					for iv3 in vert.mgd_v[iv1]:
						vert.xyz[iv3] = midxyz
					
				
				'''
				for iv3 in vert.mgd_v[iv1]:
					vert.xyz[iv3] = midxyz
					vert.mgd_v[iv1].append(iv2)
				vert.mgd_v[iv1]= remove_duplicates(vert.mgd_v[iv1])
				
				for iv3 in vert.mgd_v[iv2]:
					vert.xyz[iv3] = midxyz
					vert.mgd_v[iv2].append(iv1)
				vert.mgd_v[iv2]= remove_duplicates(vert.mgd_v[iv2])
				'''
			
		#3rd, construct tri and ref. and split to quads

		for iv in range(0,nverts):
			iv1 = iv
			iv2 = (iv +1)%nverts
			
			if collapsed[iv]:
				# construct tri
				tri = []
				tri.append(vert.xyz[face.f_to_v[iface_new][iv1]])
				tri.append(vert.xyz[face.f_to_v[iface][iv1]])
				tri.append(vert.xyz[face.f_to_v[iface][iv2]])
				
				tri_mid = []
				mid =  line_split(vert.xyz[face.f_to_v[iface_new][iv1]],vert.xyz[face.f_to_v[iface][iv1]],0.5)
				tri_mid.append(mid)
				tri_mid.append(vert.xyz[face.f_to_mid[iface][iv1]])
				
				mid =  line_split(vert.xyz[face.f_to_v[iface_new][iv1]],vert.xyz[face.f_to_v[iface][iv2]],0.5)
				tri_mid.append(mid)
				
				new_quads = quadratic_tri_to_3quads(tri,tri_mid)
			
				for i in range(0,3):
					quads.append(new_quads[i])
				
			else:
				# construct rec
				rec = []
				rec.append(vert.xyz[face.f_to_v[iface_new][iv1]])
				rec.append(vert.xyz[face.f_to_v[iface][iv1]])
				rec.append(vert.xyz[face.f_to_v[iface][iv2]])
				rec.append(vert.xyz[face.f_to_v[iface_new][iv2]])
				
				rec_mid = []
				
				mid =  line_split(vert.xyz[face.f_to_v[iface_new][iv1]],vert.xyz[face.f_to_v[iface][iv1]],0.5)
				rec_mid.append(mid)
				
				rec_mid.append(vert.xyz[face.f_to_mid[iface][iv1]])
				
				mid =  line_split(vert.xyz[face.f_to_v[iface_new][iv2]],vert.xyz[face.f_to_v[iface][iv2]],0.5)
				rec_mid.append(mid)
				
				mid =  line_split(vert.xyz[face.f_to_v[iface_new][iv1]],vert.xyz[face.f_to_v[iface_new][iv2]],0.5)
				rec_mid.append(mid)
		
				new_quads = quadratic_quad_to_4quads(rec,rec_mid)
				for i in range(0,4):
					quads.append(new_quads[i])
		
	
	return quads,iface_new
	
def split_chamfer(iface,face,vert):
	#
	quads = []
	nedges = len(face.f_to_v[iface])
	if nedges <=7 :
		for iedge in range(0,nedges):
			tri = []
			# 1st point of tri is always the face center
			trix= face.f_center[iface][0]
			triy= face.f_center[iface][1]
			triz= face.f_center[iface][2]
			tri.append([trix,triy,triz])
					
			v1 = face.f_to_e[iface][iedge][0]
			trix = vert.xyz[v1][0]
			triy = vert.xyz[v1][1]
			triz = vert.xyz[v1][2]
			tri.append([trix,triy,triz])
					
			v2 = face.f_to_e[iface][iedge][1]
			trix = vert.xyz[v2][0]
			triy = vert.xyz[v2][1]
			triz = vert.xyz[v2][2]
			tri.append([trix,triy,triz])
			
			tri_mid = []
			m1fc = line_split(tri[0],tri[1],0.5)
			tri_mid.append(m1fc)
			
			vm12=face.f_to_mid[iface][iedge]
			vm12xyz = vert.xyz[vm12]
			tri_mid.append(vm12xyz)
			
			m2fc = line_split(tri[0],tri[2],0.5)
			tri_mid.append(m2fc)
			
			new_quads = quadratic_tri_to_3quads(tri,tri_mid)
			
			for i in range(0,3):
				quads.append(new_quads[i])
			
	return quads
	
def split_chamfer_octagon(iface,face,vert):
	quads = []
	nedges = len(face.f_to_v[iface])
	if nedges == 8:
		vlist = face.f_to_v[iface]
		vmidlist = face.f_to_mid[iface]

		m = []
		v = []
		for i in range(0,8):
			v.append(vert.xyz[vlist[i]])
			m.append(vert.xyz[vmidlist[i]])
			
		c = []
		c01 = line_split(v[0],v[3],0.2)
		c02 = line_split(v[1],v[6],0.2)
		c1 = line_split(c01,c02,0.5)
		c.append(c1)
		c01 = line_split(v[0],v[3],0.8)
		c02 = line_split(v[2],v[5],0.2)
		c2 = line_split(c01,c02,0.5)
		c.append(c2)
		c01 = line_split(v[7],v[4],0.8)
		c02 = line_split(v[2],v[5],0.8)
		c3 = line_split(c01,c02,0.5)
		c.append(c3)
		c01 = line_split(v[7],v[4],0.2)
		c02 = line_split(v[1],v[6],0.8)
		c4 = line_split(c01,c02,0.5)
		c.append(c4)
		
		#c = redistribute_quad(c)
		
		# from now, split octogon into quads and tris, 
		# then split into final quads
		
		# 1. 
		rec = [c[0],v[1],v[2],c[1]]
		rec_mid = []
		mid = line_split(c[0],v[1],0.5)
		rec_mid.append(mid)
		rec_mid.append(m[1])
		mid = line_split(c[1],v[2],0.5)
		rec_mid.append(mid)
		mid = line_split(c[0],c[1],0.5)
		rec_mid.append(mid)
		
		new_quads = quadratic_quad_to_4quads(rec,rec_mid)
		for i in range(0,4):
			quads.append(new_quads[i])
		
		# 2. 
		rec = [c[2],c[1],v[3],v[4]]
		rec_mid = []
		mid = line_split(c[2],c[1],0.5)
		rec_mid.append(mid)
		mid = line_split(v[3],c[1],0.5)
		rec_mid.append(mid)
		rec_mid.append(m[3])		
		mid = line_split(c[2],v[4],0.5)
		rec_mid.append(mid)
				
		new_quads = quadratic_quad_to_4quads(rec,rec_mid)
		for i in range(0,4):
			quads.append(new_quads[i])
		
		# 3. 
		rec = [v[6],c[3],c[2],v[5]]
		rec_mid= []
		mid = line_split(v[6],c[3],0.5)
		rec_mid.append(mid)
		mid = line_split(c[2],c[3],0.5)
		rec_mid.append(mid)
		mid = line_split(c[2],v[5],0.5)
		rec_mid.append(mid)
		rec_mid.append(m[5])		
		
		new_quads = quadratic_quad_to_4quads(rec,rec_mid)
		for i in range(0,4):
			quads.append(new_quads[i])
		
		# 4. 
		rec = [v[7],v[0],c[0],c[3]]
		rec_mid = []
		rec_mid.append(m[7])
		mid = line_split(v[0],c[0],0.5)
		rec_mid.append(mid)
		mid = line_split(c[0],c[3],0.5)
		rec_mid.append(mid)
		mid = line_split(c[3],v[7],0.5)
		rec_mid.append(mid)
		
		new_quads = quadratic_quad_to_4quads(rec,rec_mid)
		for i in range(0,4):
			quads.append(new_quads[i])
		
		#5
		rec = [c[3],c[0],c[1],c[2]]
		rec_mid = []
		mid = line_split(c[3],c[0],0.5)
		rec_mid.append(mid)
		mid = line_split(c[0],c[1],0.5)
		rec_mid.append(mid)
		mid = line_split(c[1],c[2],0.5)
		rec_mid.append(mid)
		mid = line_split(c[2],c[3],0.5)
		rec_mid.append(mid)
			
		new_quads = quadratic_quad_to_4quads(rec,rec_mid)
		for i in range(0,4):
			quads.append(new_quads[i])
		
		#1 
		tri = [v[0],v[1],c[0]]
		tri_mid = []
		tri_mid.append(m[0])
		mid = line_split(c[0],v[1],0.5)
		tri_mid.append(mid)
		mid = line_split(c[0],v[0],0.5)
		tri_mid.append(mid)
		
		new_quads = quadratic_tri_to_3quads(tri,tri_mid)	
		for i in range(0,3):
			quads.append(new_quads[i])
		
		#2 
		tri = [v[2],v[3],c[1]]
		tri_mid = []
		tri_mid.append(m[2])
		mid = line_split(c[1],v[3],0.5)
		tri_mid.append(mid)
		mid = line_split(c[1],v[2],0.5)
		tri_mid.append(mid)
		
		new_quads = quadratic_tri_to_3quads(tri,tri_mid)	
		for i in range(0,3):
			quads.append(new_quads[i])
		
				
		#3
		tri = [v[4],v[5],c[2]]
		tri_mid = []
		tri_mid.append(m[4])
		mid = line_split(c[2],v[5],0.5)
		tri_mid.append(mid)
		mid = line_split(c[2],v[4],0.5)
		tri_mid.append(mid)
		
		new_quads = quadratic_tri_to_3quads(tri,tri_mid)	
		for i in range(0,3):
			quads.append(new_quads[i])
		
		#4
		tri = [v[6],v[7],c[3]]
		tri_mid = []
		tri_mid.append(m[6])
		mid = line_split(c[3],v[7],0.5)
		tri_mid.append(mid)
		mid = line_split(c[3],v[6],0.5)
		tri_mid.append(mid)
		
		new_quads = quadratic_tri_to_3quads(tri,tri_mid)	
		for i in range(0,3):
			quads.append(new_quads[i])

	return quads
	
def quadratic_quad_to_4quads(quad,quad_mid):

	fc = four_to_one(quad_mid[0],quad_mid[1],quad_mid[2],quad_mid[3])

	four_quads = []
	
	four_quads.append([quad[0],quad_mid[0],fc,quad_mid[3]])
	four_quads.append([quad[1],quad_mid[1],fc,quad_mid[0]])
	four_quads.append([quad[2],quad_mid[2],fc,quad_mid[1]])
	four_quads.append([quad[3],quad_mid[3],fc,quad_mid[2]])

	return four_quads
	
def quadratic_tri_to_3quads(tri,tri_mid):

	fc = three_to_one(tri_mid[0],tri_mid[1],tri_mid[2])
	
	three_quads = []
	
	three_quads.append([tri[0],tri_mid[0],fc,tri_mid[2]])
	three_quads.append([tri[1],tri_mid[1],fc,tri_mid[0]])
	three_quads.append([tri[2],tri_mid[2],fc,tri_mid[1]])
	
	return three_quads

def polygon_add_linear_mid(iface,face,vert):
	#add linear edge mid for polygon
	nverts = len(face.f_to_v[iface])
	vmlist = []
	for iv in range(0,nverts):				
		iv1 = iv
		iv2 = (iv +1)%nverts
		iv1 = face.f_to_v[iface][iv1]
		iv2 = face.f_to_v[iface][iv2]
		v1xyz = vert.xyz[iv1]
		v2xyz = vert.xyz[iv2]
				
		midxyz = line_split(v1xyz,v2xyz,0.5)
		vert.new_vertice(midxyz[0],midxyz[1],midxyz[2])
		
		vmlist.append(vert.nv-1)
		
	face.f_to_mid[iface] = vmlist
	
	return