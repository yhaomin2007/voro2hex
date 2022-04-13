import os
#=================================================================================
# user defined files
from data_class import *
from merge_vertices import *
from splitting_subroutines import *
from polygon_divide import *
from write_rea import *
from dump_vtk import *
from fix_elements import *
from hex_jacobian_subroutines import *

def polygons_to_hexs(vor_face,vor_vert,random_pebbles,pebble_diameter,geo,start,end,valid_faces,ip):
	
	cyl_radius = geo[0]
	cyl_bot = geo[1]
	cyl_top = geo[2]
	r_chamfer = geo[3]
	scaling = geo[4]
	n_radius = geo[5]
	bl_thickness = geo[6]
	domain_dn_plane_offset = geo[7]
	domain_up_plane_offset = geo[8]
	nlayers_dn = geo[9]
	nlayers_up = geo[10]
	delta_sw = geo[11]
	
	target_dn_plane_z = cyl_bot - domain_dn_plane_offset*pebble_diameter
	target_up_plane_z = cyl_top + domain_up_plane_offset*pebble_diameter
	#=================================================================================
	# split polygon into quads
	if ip == 0: print 'starting generating quads'

	vor_quads = voro_quads('voro_quads')

	#for iface in range(start,end):
	for iface in valid_faces[start:end]:
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
				if maxAngle <= 120:
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
						fc,moved = move_away_from_pebbles(fc,pxyz,pebble_diameter,1.03,ptag)
						
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
					cyl_radius_no_bl = cyl_radius - 0.01
					pmiddle = proj_vertice_to_cyl(pmiddle,cyl_radius_no_bl)
				
				
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
				
				
	os.system('rm -rf proc'+str(ip))
	os.system('mkdir proc'+str(ip))
	
	vtkFileName = 'proc'+str(ip)+'/voro_quad.vtk'
	dump_quad_vtk(vor_quads,vtkFileName)

	bot_quads = voro_quads('bot_quads')			
	top_quads = voro_quads('top_quads')
	vor_hex8s= voro_hex8s('voro_hex8s')

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
		
		bot_quads.add_linear_mid(bot_quads.nquads-1)
		
		
		for ipebble in range(0,2):
			pebble_id = vor_quads.quad_to_p[iquad][ipebble]
			ifghost = random_pebbles.p_ghost[pebble_id]
			ptag = random_pebbles.tag[pebble_id]
			if (not ifghost) or (ptag == 'INT'):
				# only construct hex for non-ghost pebble
				pxyz = random_pebbles.xyz[pebble_id]

				prj_center = random_pebbles.prj_center[pebble_id]
				tol = 0.01 # 0.0075
				#if (ptag == 'INT'):
					# for intersection pebbles, project to pebble surface, but with a different projecting center
				#	quad_top = project_to_pebble_for_interecting_pebbles(bot_quads.xyz[bot_quads.nquads-1],pxyz,prj_center,pebble_diameter,tol)
				#else:
					
				quad_top = project_to_pebble(bot_quads.xyz[bot_quads.nquads-1],pxyz,pebble_diameter,tol)
				
				top_quads.new_quad(quad_top,pebble_id,pebble_id)
				top_quads.edge_tag[top_quads.nquads-1] = vor_quads.edge_tag[iquad]
				top_quads.add_linear_mid(top_quads.nquads-1)
						
				qbot_id = bot_quads.nquads-1
				qtop_id = top_quads.nquads-1
				vor_hex8s.new_hex8(qbot_id,qtop_id,pebble_id)
		
				if (not linearFlag):
					bot_m = bot_quads.mid_xyz[qbot_id]
					#if (ptag == 'INT'):
					#	top_quads.mid_xyz[qtop_id] = project_to_pebble_for_interecting_pebbles(bot_m,pxyz,prj_center,pebble_diameter,tol)
					#else:
						
					top_quads.mid_xyz[qtop_id] = project_to_pebble(bot_m,pxyz,pebble_diameter,tol)
							
	
	for iquad in range (0,bot_quads.nquads):
		p1 = bot_quads.quad_to_p[iquad][0]
		p2 = bot_quads.quad_to_p[iquad][1]
		p1tag = random_pebbles.tag[p1]
		p2tag =	random_pebbles.tag[p2]
		qtag = 'E  '
		
		#if p1tag != '':  qtag = p1tag
		#if p2tag != '':  qtag = p2tag
		
		if 'SW ' in p12tag: qtag = 'SW '
		if 'TOP' in p12tag: qtag = 'TOP'
		if 'BOT' in p12tag: qtag = 'BOT'
		
		bot_quads.set_tag(iquad,qtag)
	
		
	vtkFileName = 'proc'+str(ip)+'/top_quad.vtk'
	dump_quad_2nd_vtk(top_quads,vtkFileName)
				
	vtkFileName = 'proc'+str(ip)+'/bot_quad.vtk'
	dump_quad_2nd_vtk(bot_quads,vtkFileName)
	
	if ip == 0: print 'generating hex20 elements '
	
	# from voro_hex8 convert to nek hex20 elements
	hex20 = nek_hex20s('nek_hex20s')
	
	# construct n_radius layers from polygon to sphere...
	# 1. construct boundary layer element
	# 2. construct internal layer element
	n_radius_no_bl = n_radius - 1
	
	# for internal hex20 elements.
	for ihex8 in range(0,vor_hex8s.nhex8s):
		qbot_id = vor_hex8s.quad_bot[ihex8]
		qtop_id = vor_hex8s.quad_top[ihex8]

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
		
		for ilayer in range(0,n_radius_no_bl):
			#
			# quad_m, quad_mm
			# bot_quads.xyz, bot_quads.mid_xyz
			#
			# 
			alpha1 = (float(n_radius_no_bl)-float(ilayer)-1.0)/float(n_radius_no_bl)
			alpha2 = (float(n_radius_no_bl)-float(ilayer))/float(n_radius_no_bl)
			
			v8 = []
			for iv in range(0,4):
				new_vert = line_split(bot_quads.xyz[qbot_id][iv],quad_m[iv],alpha1)
				v8.append(new_vert)
			for iv in range(0,4):
				new_vert = line_split(bot_quads.xyz[qbot_id][iv],quad_m[iv],alpha2)
				v8.append(new_vert)
			
			e12 = []
			for ie in range(0,4):
				new_vert = line_split(bot_quads.mid_xyz[qbot_id][ie],quad_mm[ie],alpha1)
				e12.append(new_vert)
			for ie in range(0,4):
				new_vert = line_split(bot_quads.mid_xyz[qbot_id][ie],quad_mm[ie],alpha2)
				e12.append(new_vert)
			for ie in range(0,4):
				mx = (v8[ie][0]+v8[ie+4][0])/2.0
				my = (v8[ie][1]+v8[ie+4][1])/2.0
				mz = (v8[ie][2]+v8[ie+4][2])/2.0
				new_vert = [mx,my,mz]
				e12.append(new_vert)
			
			s6 = []
			for iface in range(0,4):
				s6.append(bot_quads.edge_tag[qbot_id][iface])
			s6.append('E  ')
			s6.append('E  ')
			hex20.new_hex(v8,e12,s6)

	#===========================================================================
	if ip == 0:	print 'generate unstream, downstream, and cylinder sidewall hexs'
	
	#target_dn_plane_z = 0.3*scaling
	#target_up_plane_z = 1.0*scaling
	
	#target_dn_plane_z = cyl_bot - 2.0*pebble_diameter
	#target_up_plane_z = cyl_top + 4.0*pebble_diameter
	#
	#nlayers_dn = 20
	#nlayers_up = 30
	#
	#delta_sw = 0.01
	
	# to help fast extact dn and up corner quads

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
					delta_up = (target_up_plane_z-zbot)/nlayers_up
					
					offset = [0,0,delta_up*ilayer]
					vert_new =offset_vert(bot_quads.xyz[iquad][iv],offset)
					v8.append(vert_new)
				for iv in range(0,4):
					zbot = bot_quads.xyz[iquad][iv][2]
					delta_up = (target_up_plane_z-zbot)/nlayers_up
					
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
			
	vtkFileName = 'proc'+str(ip)+'/top_side_quads.vtk'
	dump_quad_2nd_vtk(top_side_quads,vtkFileName)
	vtkFileName = 'proc'+str(ip)+'/bot_side_quads.vtk'
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
	
	hex20.fix_non_right_hand_elements(ip)
	
	vtkFileName = 'proc'+str(ip)+'/hex20_lar.vtk'
	dump_hex_2nd_vtk_lar(hex20,vtkFileName)
	
	vtkFileName = 'proc'+str(ip)+'/hex20_nnrh.vtk'
	dump_hex_2nd_vtk_nnrh(hex20,vtkFileName)
	
	#hex20.check_max_aspect_ratio()
	
	#dump_hex_2nd_vtk_nnrh(hex20,'hex20_har.vtk')
	
	#hex20.check_non_right_hand_elements_use_nek_method()
	
	#dump_hex_vtk(hex20,'hex8.vtk')
	
	vtkFileName = 'proc'+str(ip)+'/hex20.vtk'
	dump_hex_2nd_vtk(hex20,vtkFileName)
	
	#dump_hex_2nd_vtk_nnrh(hex20,'hex20_nnrh.vtk')
	#=================================================================================
	# part 4
	# generate rea file based on hex20 information.
	# no edge curvature right now
	
	# explicit fix neg-jacobian element by element number provided by nek.
	# fixing it by linearize this element.
	# some mesh does not need this step.....
	#fix_negative_jacobian_elements_from_nek(hex20)
	
	#if ip == 0: print 'wrting new rea files'
	#baseRea = 'base.rea'
	#newReaFile =  'proc'+str(ip)+'/pb.rea'
	#newReaFile = 'pb'+str(ip)+'.rea'
	#write_rea(baseRea,newReaFile,hex20)
	
	file =  'proc'+str(ip)+'/hfile.dat'
	write_hfiles(file,hex20)

#=================================================================================
# use old method to generate bl for top and bot extrueded mesh. 

def polygons_to_hexs2(vor_face,vor_vert,random_pebbles,pebble_diameter,geo,start,end,valid_faces,ip):
	
	cyl_radius = geo[0]
	cyl_bot = geo[1]
	cyl_top = geo[2]
	r_chamfer = geo[3]
	scaling = geo[4]
	n_radius = geo[5]
	bl_thickness = geo[6]
	domain_dn_plane_offset = geo[7]
	domain_up_plane_offset = geo[8]
	nlayers_dn = geo[9]
	nlayers_up = geo[10]
	delta_sw = geo[11]
	
	target_dn_plane_z = cyl_bot - domain_dn_plane_offset*pebble_diameter
	target_up_plane_z = cyl_top + domain_up_plane_offset*pebble_diameter
	#=================================================================================
	# split polygon into quads
	if ip == 0: print 'starting generating quads'

	vor_quads = voro_quads('voro_quads')

	#for iface in range(start,end):
	for iface in valid_faces[start:end]:
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
				if maxAngle <= 120:
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
						fc,moved = move_away_from_pebbles(fc,pxyz,pebble_diameter,1.03,ptag)
						
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
					cyl_radius_no_bl = cyl_radius - 0.01
					pmiddle = proj_vertice_to_cyl(pmiddle,cyl_radius_no_bl)
				
				
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
				
				
	os.system('rm -rf proc'+str(ip))
	os.system('mkdir proc'+str(ip))
	
	vtkFileName = 'proc'+str(ip)+'/voro_quad.vtk'
	dump_quad_vtk(vor_quads,vtkFileName)
	
	bot_quads = voro_quads('bot_quads')			
	top_quads = voro_quads('top_quads')
	vor_hex8s= voro_hex8s('voro_hex8s')

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
		
		bot_quads.add_linear_mid(bot_quads.nquads-1)
		
		
		for ipebble in range(0,2):
			pebble_id = vor_quads.quad_to_p[iquad][ipebble]
			ifghost = random_pebbles.p_ghost[pebble_id]
			ptag = random_pebbles.tag[pebble_id]
			if (not ifghost) or (ptag == 'INT'):
				# only construct hex for non-ghost pebble
				pxyz = random_pebbles.xyz[pebble_id]
				
				prj_center = random_pebbles.prj_center[pebble_id]
				tol = 0.01 # 0.0075
				#if (ptag == 'INT'):
					# for intersection pebbles, project to pebble surface, but with a different projecting center
				#	quad_top = project_to_pebble_for_interecting_pebbles(bot_quads.xyz[bot_quads.nquads-1],pxyz,prj_center,pebble_diameter,tol)
				#else:
					
				quad_top = project_to_pebble(bot_quads.xyz[bot_quads.nquads-1],pxyz,pebble_diameter,tol)
				
				top_quads.new_quad(quad_top,pebble_id,pebble_id)
				top_quads.edge_tag[top_quads.nquads-1] = vor_quads.edge_tag[iquad]
				top_quads.add_linear_mid(top_quads.nquads-1)
						
				qbot_id = bot_quads.nquads-1
				qtop_id = top_quads.nquads-1
				vor_hex8s.new_hex8(qbot_id,qtop_id,pebble_id)
		
				if (not linearFlag):
					bot_m = bot_quads.mid_xyz[qbot_id]
					#if (ptag == 'INT'):
					#	# for intersection pebbles, project to pebble surface, but with a different projecting center
					#	top_quads.mid_xyz[qtop_id] = project_to_pebble_for_interecting_pebbles(bot_m,pxyz,prj_center,pebble_diameter,tol)
					#else:
							
					top_quads.mid_xyz[qtop_id] = project_to_pebble(bot_m,pxyz,pebble_diameter,tol)
							
	
	for iquad in range (0,bot_quads.nquads):
		p1 = bot_quads.quad_to_p[iquad][0]
		p2 = bot_quads.quad_to_p[iquad][1]
		p1tag = random_pebbles.tag[p1]
		p2tag =	random_pebbles.tag[p2]
		p12tag = [p1tag,p2tag]
		qtag = 'E  '
		
		#if p1tag != '':  qtag = p1tag
		#if p2tag != '':  qtag = p2tag
		
		if 'SW ' in p12tag: qtag = 'SW '
		if 'TOP' in p12tag: qtag = 'TOP'
		if 'BOT' in p12tag: qtag = 'BOT'
		
		bot_quads.set_tag(iquad,qtag)
	
		
	vtkFileName = 'proc'+str(ip)+'/top_quad.vtk'
	dump_quad_2nd_vtk(top_quads,vtkFileName)
				
	vtkFileName = 'proc'+str(ip)+'/bot_quad.vtk'
	dump_quad_2nd_vtk(bot_quads,vtkFileName)
	
	if ip == 0: print 'generating hex20 elements '
	
	# from voro_hex8 convert to nek hex20 elements
	
	hex20 = nek_hex20s('nek_hex20s')
	
	# construct n_radius layers from polygon to sphere...
	# 1. construct boundary layer element
	# 2. construct internal layer element
	n_radius_no_bl = n_radius - 1
	
	
	# for internal hex20 elements.
	for ihex8 in range(0,vor_hex8s.nhex8s):
		qbot_id = vor_hex8s.quad_bot[ihex8]
		qtop_id = vor_hex8s.quad_top[ihex8]

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
	
		for ilayer in range(0,n_radius_no_bl):
			#
			# quad_m, quad_mm
			# bot_quads.xyz, bot_quads.mid_xyz
			#
			# 
			alpha1 = (float(n_radius_no_bl)-float(ilayer)-1.0)/float(n_radius_no_bl)
			alpha2 = (float(n_radius_no_bl)-float(ilayer))/float(n_radius_no_bl)
			
			v8 = []
			for iv in range(0,4):
				new_vert = line_split(bot_quads.xyz[qbot_id][iv],quad_m[iv],alpha1)
				v8.append(new_vert)
			for iv in range(0,4):
				new_vert = line_split(bot_quads.xyz[qbot_id][iv],quad_m[iv],alpha2)
				v8.append(new_vert)
			
			e12 = []
			for ie in range(0,4):
				new_vert = line_split(bot_quads.mid_xyz[qbot_id][ie],quad_mm[ie],alpha1)
				e12.append(new_vert)
			for ie in range(0,4):
				new_vert = line_split(bot_quads.mid_xyz[qbot_id][ie],quad_mm[ie],alpha2)
				e12.append(new_vert)
			for ie in range(0,4):
				mx = (v8[ie][0]+v8[ie+4][0])/2.0
				my = (v8[ie][1]+v8[ie+4][1])/2.0
				mz = (v8[ie][2]+v8[ie+4][2])/2.0
				new_vert = [mx,my,mz]
				e12.append(new_vert)
			
			s6 = []
			for iface in range(0,4):
				s6.append(bot_quads.edge_tag[qbot_id][iface])
			s6.append('E  ')
			s6.append('E  ')
			hex20.new_hex(v8,e12,s6)

	#===========================================================================
	if ip == 0:	print 'generate unstream, downstream, and cylinder sidewall hexs'
	
	#target_dn_plane_z = 0.3*scaling
	#target_up_plane_z = 1.0*scaling
	
	#target_dn_plane_z = cyl_bot - 2.0*pebble_diameter
	#target_up_plane_z = cyl_top + 4.0*pebble_diameter
	#
	#nlayers_dn = 20
	#nlayers_up = 30
	#
	#delta_sw = 0.01
	
	# to help fast extact dn and up corner quads

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
					delta_up = (target_up_plane_z-zbot)/nlayers_up
					
					offset = [0,0,delta_up*ilayer]
					vert_new =offset_vert(bot_quads.xyz[iquad][iv],offset)
					v8.append(vert_new)
				for iv in range(0,4):
					zbot = bot_quads.xyz[iquad][iv][2]
					delta_up = (target_up_plane_z-zbot)/nlayers_up
					
					offset = [0,0,delta_up*(ilayer+1)]
					vert_new =offset_vert(bot_quads.xyz[iquad][iv],offset)
					v8.append(vert_new)
				
				e12 = v8_to_e12(v8)	
				if ilayer == (nlayers_up-1): s6[5] = 'TOP'
				hex20.new_hex(v8,e12,s6)
				
				
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
				#vert_new = line_split_for_boundary_layer(vbot,cxyz,-delta_sw)
				# direct project to sidewall.
				r_here = distance(vbot,cxyz)
				delta_sw = cyl_radius - r_here
				vert_new = line_split_for_boundary_layer(vbot,cxyz,-delta_sw)
				
				v8.append(vert_new)
			
			e12 = v8_to_e12(v8)
			
			hex20.new_hex(v8,e12,s6)
			
			
			
			#========================================================================
			ihex = hex20.nhex-1
			for iface in range(0,4):
				vc,fn,vquad,equad,edge_tag = hex20.return_face_info(ihex,iface)

				if_top_side_quads = False				
				if_bot_side_quads = False
					
				if (vc[2] > (cyl_top-0.001)): if_top_side_quads = True
				if (vc[2] < (cyl_bot+0.001)): if_bot_side_quads = True
				
				# sidewall
				if (if_top_side_quads):
					# add to top_side_quads
					top_side_quads.new_quad(vquad,0,0)
					top_side_quads.edge_tag[top_side_quads.nquads-1] = edge_tag
					top_side_quads.add_mid(top_side_quads.nquads-1,equad)

				if (if_bot_side_quads):
					# add to bot_side_quads
					bot_side_quads.new_quad(vquad,0,0)
					bot_side_quads.edge_tag[bot_side_quads.nquads-1] = edge_tag
					bot_side_quads.add_mid(bot_side_quads.nquads-1,equad)
			
			
	vtkFileName = 'proc'+str(ip)+'/top_side_quads.vtk'
	dump_quad_2nd_vtk(top_side_quads,vtkFileName)
	vtkFileName = 'proc'+str(ip)+'/bot_side_quads.vtk'
	dump_quad_2nd_vtk(bot_side_quads,vtkFileName)	
		
	# based on top_side_quads and bot_side_quads

	for iquad in range (0,top_side_quads.nquads):

		for ilayer in range (0,nlayers_up):
	
			s6 = []
			for iface in range(0,6):
				s6.append('E  ')
			
			for iface in range(0,4):
				s6[iface] = top_side_quads.edge_tag[iquad][iface]
	
			v8 = []

			for iv in range(0,4):
							
				zbot = top_side_quads.xyz[iquad][iv][2]
				delta_up = (target_up_plane_z-zbot)/nlayers_up
				
				offset = [0,0,delta_up*ilayer]
				vert_new =offset_vert(top_side_quads.xyz[iquad][iv],offset)
				v8.append(vert_new)
				
			for iv in range(0,4):
				
				zbot = top_side_quads.xyz[iquad][iv][2]
				delta_up = (target_up_plane_z-zbot)/nlayers_up
				
				offset = [0,0,delta_up*(ilayer+1)]
				vert_new =offset_vert(top_side_quads.xyz[iquad][iv],offset)
				v8.append(vert_new)
			
			e12 = v8_to_e12(v8)
		
			if ilayer == (nlayers_up-1): s6[5] = 'TOP'
			
			hex20.new_hex(v8,e12,s6)

	for iquad in range (0,bot_side_quads.nquads):

		for ilayer in range (0,nlayers_dn):
	
			s6 = []
			for iface in range(0,6):
				s6.append('E  ')
			
			for iface in range(0,4):
				s6[iface] = bot_side_quads.edge_tag[iquad][iface]
	
			v8 = []
			for iv in range(0,4):
				
				zbot = bot_side_quads.xyz[iquad][iv][2]
				delta_dn = (zbot-target_dn_plane_z)/nlayers_dn
				
				offset = [0,0,-delta_dn*ilayer]
				vert_new =offset_vert(bot_side_quads.xyz[iquad][iv],offset)
				v8.append(vert_new)
			for iv in range(0,4):
				
				zbot = bot_side_quads.xyz[iquad][iv][2]
				delta_dn = (zbot-target_dn_plane_z)/nlayers_dn
				
				offset = [0,0,-delta_dn*(ilayer+1)]
				vert_new =offset_vert(bot_side_quads.xyz[iquad][iv],offset)
				v8.append(vert_new)
		
			e12 = v8_to_e12(v8)
			
			if ilayer == (nlayers_dn-1): s6[5] = 'BOT'
			hex20.new_hex(v8,e12,s6)
	
	#================================================
	
	hex20.fix_non_right_hand_elements(ip)
	
	vtkFileName = 'proc'+str(ip)+'/hex20_lar.vtk'
	dump_hex_2nd_vtk_lar(hex20,vtkFileName)
	
	vtkFileName = 'proc'+str(ip)+'/hex20_nnrh.vtk'
	dump_hex_2nd_vtk_nnrh(hex20,vtkFileName)
	
	#hex20.check_max_aspect_ratio()
	
	#dump_hex_2nd_vtk_nnrh(hex20,'hex20_har.vtk')
	
	#hex20.check_non_right_hand_elements_use_nek_method()
	
	vtkFileName = 'proc'+str(ip)+'/hex8.vtk'
	dump_hex_vtk(hex20,vtkFileName)
	
	#vtkFileName = 'proc'+str(ip)+'/hex20.vtk'
	#dump_hex_2nd_vtk(hex20,vtkFileName)
	
	#dump_hex_2nd_vtk_nnrh(hex20,'hex20_nnrh.vtk')
	#=================================================================================
	# part 4
	# generate rea file based on hex20 information.
	# no edge curvature right now
	
	# explicit fix neg-jacobian element by element number provided by nek.
	# fixing it by linearize this element.
	# some mesh does not need this step.....
	#fix_negative_jacobian_elements_from_nek(hex20)
	
	#if ip == 0: print 'wrting new rea files'
	#baseRea = 'base.rea'
	#newReaFile =  'proc'+str(ip)+'/pb.rea'
	#newReaFile = 'pb'+str(ip)+'.rea'
	#write_rea(baseRea,newReaFile,hex20)
	
	file =  'proc'+str(ip)+'/hfile.dat'
	write_hfiles(file,hex20)
	
def write_hfiles(file,hex20):
	
	fileHolder = open(file,'w')
	nhexs = hex20.nhex
	print 'write_hfiles: dumping '+str(nhexs) + ' elements'
	newLine = str(nhexs) + '\n'
	fileHolder.write(newLine)
	
	for ihex in range(2,nhexs):
		v8 = hex20.v8[ihex]
		e12 =  hex20.e12[ihex]
		
		if_auto_correction = True
		if (if_auto_correction):
			if_linearize = detect_jacobian_for_hex20(v8,e12)
			if (if_linearize):
				hex20.e12[ihex] = v8_to_e12(hex20.v8[ihex])
				hex20.e12[ihex-1] = v8_to_e12(hex20.v8[ihex-1])
				hex20.e12[ihex-2] = v8_to_e12(hex20.v8[ihex-2])
	
	for ihex in range(0,nhexs):
		v8 = hex20.v8[ihex]
		e12 =  hex20.e12[ihex]
		s6 = hex20.s6[ihex]
				
		ns6 = []
		for iface in range(0,6):
			ss  = s6[iface]
			if (ss == 'E  '): ns6.append(0)
			if (ss == 'PW '): ns6.append(1)
			if (ss == 'C  '): ns6.append(2)
			if (ss == 'SW '): ns6.append(3)
			if (ss == 'TOP'): ns6.append(4)
			if (ss == 'BOT'): ns6.append(5)
	
		newLine = ''
		for iv in range(0,8):
			for i in range(0,3):
				newLine = newLine + str(v8[iv][i]) + ' '
		newLine = newLine + '\n'
		fileHolder.write(newLine)
		
		newLine = ''
		for ie in range(0,12):
			for i in range(0,3):
				newLine = newLine + str(e12[ie][i]) + ' '
		newLine = newLine + '\n'
		fileHolder.write(newLine)
			
		newLine = ''
		for iface in range(0,6):
			newLine = newLine + str(ns6[iface]) + ' '
		newLine = newLine + '\n'
		fileHolder.write(newLine)
	
	fileHolder.close()
	print 'DONE: dumping '+str(nhexs) + ' elements'
