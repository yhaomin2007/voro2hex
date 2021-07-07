import math
import multiprocessing
import os

from splitting_subroutines import *

def build_connectivity_for_quads(quads,pebbles,start,end,iproc,filetag):
	print 'build_connectivity_for_quads'	
	dummy = []
	
	for iquad in range(start,end):
		ip1= quads.quad_to_p[iquad][0]
		ip2= quads.quad_to_p[iquad][1]
		same_pebble_quads_list1 = pebbles.voro_quads[ip1]
		same_pebble_quads_list2 = pebbles.voro_quads[ip2]
		
		quads_list = same_pebble_quads_list1
		quads_list.extend(same_pebble_quads_list2)
		
		for ivert in range(0,4):
			vxyz = quads.xyz[iquad][ivert]
			# store self first
			dummy.append([iquad,ivert,iquad,ivert])
			
			# find all quads,vertice that share the same location
			for iquad2 in quads_list:
				if iquad2 != iquad:
					verts2 = quads.xyz[iquad2]
					for i in range(0,4):
						dist = distance(vxyz,verts2[i])
						if dist < 1e-4:
							dummy.append([iquad,ivert,iquad2,i])
	
	# writing to a file 
	filename = 'dummy_'+filetag+'/dummy_'+str(iproc)
	print 'writing files to '+filename
	dummyfile = open(filename,'w')
	nlines = len(dummy)
	line = str(nlines)+'\n'
	dummyfile.write(line)
	for iline in range(0,nlines):
		line = str(dummy[iline][0]) + ' ' + str(dummy[iline][1]) + ' ' + str(dummy[iline][2]) + ' ' + str(dummy[iline][3]) +'\n'
		dummyfile.write(line)
	dummyfile.close()
	print 'done writing files to '+filename
	
def build_connectivity_for_quads_parallel(quads,pebbles,np,new_connectivity,filetag):
	print 'build_connectivity_for_quads_parallel'
	nquads = quads.nquads
	jobs = []
	start = []
	end = []
	
	nrange = int(math.ceil(nquads/np))

	for ip in range(0,np):
		start.append(ip*nrange)
		
	for ip in range(0,np-1):
		end.append(start[ip+1])	
	end.append(nquads)
	

	if new_connectivity:
	
		os.system('rm -rf dummy_'+filetag)
		os.system('mkdir dummy_'+filetag)
	
		for ip in range(0,np):
			process = multiprocessing.Process(target=build_connectivity_for_quads,args=(quads,pebbles,start[ip],end[ip],ip,filetag))
			jobs.append(process)
	
		for j in jobs:
			j.start()

		for j in jobs:
			j.join()

	# return connectivity
	# read all dummy files from 
	print 'combining all connectivity info from dummy files'
	
	for ip in range(0,np):
		filename = 'dummy_'+filetag+'/dummy_'+str(ip)
		dummyfile = open(filename)
		line = dummyfile.readline()
		nlines =  int(line.split()[0])
		for iline in range(0,nlines):
			line = dummyfile.readline()
			iquad1 = int(line.split()[0])
			iv1    = int(line.split()[1])
			iquad2 = int(line.split()[2])
			iv2    = int(line.split()[3])
			quads.connect_quads_and_vertices[iquad1][iv1].append([iquad2,iv2])
			
		dummyfile.close()
	print 'done combining all connectivity info'
	
	# merge all all connectivity info'
	# make sure there is no problem 
	
	for iquad1 in range(0,nquads):
		for iv1 in range(0,4):
			cqv1 = quads.connect_quads_and_vertices[iquad1][iv1]
			cqv1 = remove_duplicates(cqv1)
			quads.connect_quads_and_vertices[iquad1][iv1] = cqv1
	
	n = 3
	for i in range(0,n):
		for iquad1 in range(0,nquads):
			#print iquad1
			for iv1 in range(0,4):
				cqv1 = quads.connect_quads_and_vertices[iquad1][iv1]
				#print cqv1
				cqv3 = cqv1
				for nbr_quads_vertices in cqv1:
					#print nbr_quads_vertices
					iquad2 = nbr_quads_vertices[0]
					iv2 = nbr_quads_vertices[1]
					if (iquad2 != iquad1):
						cqv2 = quads.connect_quads_and_vertices[iquad2][iv2]
					#	cqv2 = remove_duplicates(cqv2)
						#print cqv2
						cqv3.extend(cqv2)
						cqv3 = remove_duplicates(cqv3)
			
				cqv1 = remove_duplicates(cqv3)
				#print cqv1
				quads.connect_quads_and_vertices[iquad1][iv1] = cqv1

		if i == (n-1):
			for iquad1 in range(0,nquads):
				for iv1 in range(0,4):
					cqv1 = quads.connect_quads_and_vertices[iquad1][iv1]
					for nbr_quads_vertices in cqv1:
						iquad2 = nbr_quads_vertices[0]
						iv2 = nbr_quads_vertices[1]
						if (iquad2 != iquad1):
							quads.connect_quads_and_vertices[iquad2][iv2] = cqv1
				
	print 'done merging all connectivity'
	

def quads_smoothing_near_and_at_chamfer_only(quads,pebbles,pdiameter):
	print 'quads_smoothing_near_and_at_chamfer_only'
	nquads = quads.nquads
	for iquad in range(0,nquads):
		for ivert in range(0,4):
			# laplacian smoothing all vertices of this quads...
			# however shold avoid vertice on chamfer
			
			ifsmooth = True
			
			ifchamfer = False
			
			tag1 =  quads.edge_tag[iquad][ivert]
			tag2 =  quads.edge_tag[iquad][(ivert-1)%4]
			if tag1=='C  ' or tag2=='C  ':
				ifchamfer = True
				
				
			p1 = quads.quad_to_p[iquad][0]
			p2 = quads.quad_to_p[iquad][1]
			p1tag = pebbles.tag[p1]
			p2tag =	pebbles.tag[p2]
			
			if p1tag == 'TOP' or p1tag == 'BOT' or p2tag == 'TOP' or p2tag == 'BOT' : 
				ifsmooth = False
			
			p1xyz = pebbles.xyz[p1]
			p2xyz = pebbles.xyz[p2]
			vxyz = quads.xyz[iquad][ivert]
			d1 = distance(vxyz,p1xyz)
			d2 = distance(vxyz,p2xyz)
			if (d1 > 1.1*(pdiameter/2.0)) and (d2 > 1.1*(pdiameter/2.0)):
				ifsmooth = False
				
			
			if ifsmooth and (not ifchamfer):

				iv1 = (ivert+1)%4
				iv2 = (ivert-1)%4
					
				v1xyz = quads.xyz[iquad][iv1]
				v2xyz = quads.xyz[iquad][iv2]
					
				avgvert = line_split(v1xyz,v2xyz,0.5)
				
				quads.newxyz[iquad][ivert] = avgvert
				
			if ifchamfer:	# move chamfer poits
				#quads.newxyz[iquad][ivert] = quads.xyz[iquad][ivert]
				
				if tag1=='C  ':
					quads.newxyz[iquad][ivert] = quads.xyz[iquad][(ivert+1)%4]
				else:
					quads.newxyz[iquad][ivert] = quads.xyz[iquad][(ivert-1)%4]

				
	for iquad in range(0,nquads):
		for ivert in range(0,4):

			quads.newxyz[iquad][ivert] = quads.xyz[iquad][ivert]
		
			ifsmooth = True
			
			ifchamfer = False
			
			tag1 =  quads.edge_tag[iquad][ivert]
			tag2 =  quads.edge_tag[iquad][(ivert-1)%4]
			if tag1=='C  ' or tag2=='C  ':
				ifchamfer = True
				
				
			p1 = quads.quad_to_p[iquad][0]
			p2 = quads.quad_to_p[iquad][1]
			p1tag = pebbles.tag[p1]
			p2tag =	pebbles.tag[p2]
			
			if p1tag == 'TOP' or p1tag == 'BOT' or p2tag == 'TOP' or p2tag == 'BOT' : 
				ifsmooth = False
			
			p1xyz = pebbles.xyz[p1]
			p2xyz = pebbles.xyz[p2]
			vxyz = quads.xyz[iquad][ivert]
			d1 = distance(vxyz,p1xyz)
			d2 = distance(vxyz,p2xyz)
			if (d1 > 1.1*(pdiameter/2.0)) and (d2 > 1.1*(pdiameter/2.0)):
				ifsmooth = False
			
			
			
			if ifsmooth and (not ifchamfer):
				avgvert = [0.0,0.0,0.0]
				cqv = quads.connect_quads_and_vertices[iquad][ivert] 
				nnbr = len(cqv)

				for nbr_quads_vertices in cqv:
					# simple average now, may use harmonic average later
					iquad2 = nbr_quads_vertices[0]
					ivert2 = nbr_quads_vertices[1]
					
					avgvert[0] = avgvert[0] + quads.newxyz[iquad2][ivert2][0]*(1.0/nnbr)
					avgvert[1] = avgvert[1] + quads.newxyz[iquad2][ivert2][1]*(1.0/nnbr)
					avgvert[2] = avgvert[2] + quads.newxyz[iquad2][ivert2][2]*(1.0/nnbr)
					
								
				for nbr_quads_vertices in cqv:
					# simple average now, may use harmonic average later
					iquad2 = nbr_quads_vertices[0]
					ivert2 = nbr_quads_vertices[1]
					
					quads.newxyz[iquad2][ivert2] = avgvert
			
			if ifchamfer:
				
				avgvert = [0.0,0.0,0.0]
				cqv = quads.connect_quads_and_vertices[iquad][ivert] 
				nnbr = len(cqv)

				for nbr_quads_vertices in cqv:
					# simple average now, may use harmonic average later
					iquad2 = nbr_quads_vertices[0]
					ivert2 = nbr_quads_vertices[1]
					
					avgvert[0] = avgvert[0] + quads.newxyz[iquad2][ivert2][0]*(1.0/nnbr)
					avgvert[1] = avgvert[1] + quads.newxyz[iquad2][ivert2][1]*(1.0/nnbr)
					avgvert[2] = avgvert[2] + quads.newxyz[iquad2][ivert2][2]*(1.0/nnbr)
				
				
				p1 = quads.quad_to_p[iquad][0]
				p1xyz = pebbles.xyz[p1]
				p2 = quads.quad_to_p[iquad][1]
				p2xyz = pebbles.xyz[p2]
				cxyz = line_split(p1xyz,p2xyz,0.5)
				# reproject to chamfer wall now
				d1 = distance(quads.xyz[iquad][ivert],cxyz)
				d2 = distance(avgvert,cxyz)
				
				avgvert = line_split(cxyz,avgvert,d1/d2)
							
				for nbr_quads_vertices in cqv:
					# simple average now, may use harmonic average later
					iquad2 = nbr_quads_vertices[0]
					ivert2 = nbr_quads_vertices[1]

					quads.newxyz[iquad2][ivert2] = avgvert
	
	
	for iquad in range(0,nquads):
		for ivert in range(0,4):
			quads.xyz[iquad][ivert] = quads.newxyz[iquad][ivert]
						
def project_quads_to_sidewall(quads,pebbles,cyl_top,cyl_bot,pdiameter,cyl_radius_no_bl):

	print 'reproject to sidewall'
	nquads = quads.nquads
	for iquad in range(0,nquads):
		p1 = quads.quad_to_p[iquad][0]
		p1xyz = pebbles.xyz[p1]
		p2 = quads.quad_to_p[iquad][1]
		p2xyz = pebbles.xyz[p2]
		
		p1tag = pebbles.tag[p1]
		p2tag =	pebbles.tag[p2]
		
		pradius  = pdiameter/2.0
		
		ifchamfer = False
		dist = distance(p1xyz,p2xyz)
		if (dist < 1.02*pdiameter):
			ifchamfer = True
		
		if (p1tag=='SW ') or (p2tag=='SW '):
			# this quad is on sidewall. need to reproject to sidewall....

			
			for ivert in range(0,4):
				vxyz = quads.xyz[iquad][ivert]	
				xx = vxyz[0]
				yy = vxyz[1]
				zz = vxyz[2]
				
				orgzz = (0,0,zz)
				d = math.sqrt(xx**2.0+yy**2.0)
				
				swxyz = line_split(orgzz,vxyz,cyl_radius_no_bl/d)
				
				cqv = quads.connect_quads_and_vertices[iquad][ivert]
				for nbr_quads_vertices in cqv:
					iquad2 = nbr_quads_vertices[0]
					iv2 = nbr_quads_vertices[1]
					quads.xyz[iquad2][iv2]	= swxyz
				
				#if (ifchamfer):
				#	cxyz = line_split(p1xyz,p2xyz,0.5)
				#	dist = distance(vxyz,cxyz)
				#	
				#	alpha = 1.0
				#	if (dist<pradius*0.5):
				#		alpha = 1.0 # complete project to cylinder
				#	elif (dist>=pradius*0.5) and (dist<pradius*0.8):
				#		alpha = 1.0 - (dist - pradius*0.5)/(pradius*0.5)
				#	elif (dist>=pradius*0.8):
				#		alpha = 0.0
				#	
				#	nvxyz = line_split(vxyz,swxyz,alpha)
		       #
				#	cqv = quads.connect_quads_and_vertices[iquad][ivert]
				#	for nbr_quads_vertices in cqv:
				#		iquad2 = nbr_quads_vertices[0]
				#		iv2 = nbr_quads_vertices[1]
				#		quads.xyz[iquad2][iv2]	= nvxyz
		       #
				#if (zz > (cyl_top-1.75*pdiameter)) or (zz < (cyl_bot+1.75*pdiameter)):
				#	
				#	cqv = quads.connect_quads_and_vertices[iquad][ivert]
				#	for nbr_quads_vertices in cqv:
				#		iquad2 = nbr_quads_vertices[0]
				#		iv2 = nbr_quads_vertices[1]
				#		quads.xyz[iquad2][iv2]	= swxyz
		   
		
def build_connectivity_for_bot_quads(quads,pebbles,start,end,iproc,filetag):
	print 'build_connectivity_for_bot_quads'	
	dummy = []
	
	for iquad in range(start,end):
		ip1= quads.quad_to_p[iquad][0]
		quads_list =  pebbles.bot_quads[ip1] 
		
		for ivert in range(0,4):
			vxyz = quads.xyz[iquad][ivert]
			# store self first
			dummy.append([iquad,ivert,iquad,ivert])
			
			# find all quads,vertice that share the same location
			for iquad2 in quads_list:
				if iquad2 != iquad:
					verts2 = quads.xyz[iquad2]
					for i in range(0,4):
						dist = distance(vxyz,verts2[i])
						if dist < 1e-4:
							dummy.append([iquad,ivert,iquad2,i])
	
	# writing to a file 
	filename = 'dummy_'+filetag+'/dummy_'+str(iproc)
	print 'writing files to '+filename
	dummyfile = open(filename,'w')
	nlines = len(dummy)
	line = str(nlines)+'\n'
	dummyfile.write(line)
	for iline in range(0,nlines):
		line = str(dummy[iline][0]) + ' ' + str(dummy[iline][1]) + ' ' + str(dummy[iline][2]) + ' ' + str(dummy[iline][3]) +'\n'
		dummyfile.write(line)
	dummyfile.close()
	print 'done writing files to '+filename
	
def build_connectivity_for_bot_quads_parallel(quads,pebbles,np,new_connectivity,filetag):
	print 'build_connectivity_for_bot_quads_parallel'
	nquads = quads.nquads
	jobs = []
	start = []
	end = []
	
	nrange = int(math.ceil(nquads/np))

	for ip in range(0,np):
		start.append(ip*nrange)
		
	for ip in range(0,np-1):
		end.append(start[ip+1])	
	end.append(nquads)
	

	if new_connectivity:
	
		os.system('rm -rf dummy_'+filetag)
		os.system('mkdir dummy_'+filetag)
	
		for ip in range(0,np):
			process = multiprocessing.Process(target=build_connectivity_for_bot_quads,args=(quads,pebbles,start[ip],end[ip],ip,filetag))
			jobs.append(process)
	
		for j in jobs:
			j.start()

		for j in jobs:
			j.join()

	# return connectivity
	# read all dummy files from 
	print 'combining all connectivity info from dummy files'
	
	for ip in range(0,np):
		filename = 'dummy_'+filetag+'/dummy_'+str(ip)
		dummyfile = open(filename)
		line = dummyfile.readline()
		nlines =  int(line.split()[0])
		for iline in range(0,nlines):
			line = dummyfile.readline()
			iquad1 = int(line.split()[0])
			iv1    = int(line.split()[1])
			iquad2 = int(line.split()[2])
			iv2    = int(line.split()[3])
			quads.connect_quads_and_vertices[iquad1][iv1].append([iquad2,iv2])
			
		dummyfile.close()
	print 'done combining all connectivity info'
	
	# merge all all connectivity info'
	# make sure there is no problem 
	for iquad1 in range(0,nquads):
		for iv1 in range(0,4):
			cqv1 = quads.connect_quads_and_vertices[iquad1][iv1]
			cqv1 = remove_duplicates(cqv1)
			quads.connect_quads_and_vertices[iquad1][iv1] = cqv1
			if(len(cqv1)<1): 
				print 'ERROR, cqv length 0 for bot quads'
			
	
	for iquad1 in range(0,nquads):
		for iv1 in range(0,4):
			cqv1 = quads.connect_quads_and_vertices[iquad1][iv1]
		
			for nbr_quads_vertices in cqv1:
				iquad2 = nbr_quads_vertices[0]
				iv2 = nbr_quads_vertices[1]
				cqv2 = quads.connect_quads_and_vertices[iquad2][iv2]
				
				if (len(cqv1) != len(cqv2)):
					print 'ERROR, cqv length did not match'
					
				if  [iquad1,iv1] not in cqv2:
					print 'ERROR, cqv1 and  cqv2 not match '
	
	
	print 'done merging all connectivity'
	
def bot_quads_smoothing(quads,pebbles,pdiameter):
	print 'bot_quads_smoothing'
	nquads = quads.nquads
	for iquad in range(0,nquads):
		for ivert in range(0,4):
			# laplacian smoothing all vertices of this quads...
			# however shold avoid vertice on chamfer
			
			notchamfer = True
			
			tag1 =  quads.edge_tag[iquad][ivert]
			tag2 =  quads.edge_tag[iquad][(ivert-1)%4]
			if tag1=='C  ' or tag2=='C  ':
				notchamfer = False
				
			#p1 = quads.quad_to_p[iquad][0]
			#p1xyz = pebbles.xyz[p1]
			#vxyz = quads.xyz[iquad][ivert]
			#d = distance(vxyz,p1xyz)
			#if d < 1.1*(pdiameter/2.0):
			#	notchamfer = False
			
			if notchamfer:
				
				iv1 = (ivert+1)%4
				iv2 = (ivert-1)%4
					
				v1xyz = quads.xyz[iquad][iv1]
				v2xyz = quads.xyz[iquad][iv2]
					
				avgvert = line_split(v1xyz,v2xyz,0.5)
				
				quads.newxyz[iquad][ivert] = avgvert

				
			else:	
				quads.newxyz[iquad][ivert] = quads.xyz[iquad][ivert]

				
	for iquad in range(0,nquads):
		for ivert in range(0,4):

			notchamfer = True
			
			tag1 =  quads.edge_tag[iquad][ivert]
			tag2 =  quads.edge_tag[iquad][(ivert-1)%4]
			if tag1=='C  ' or tag2=='C  ':
				notchamfer = False
			
			#p1 = quads.quad_to_p[iquad][0]
			#p1xyz = pebbles.xyz[p1]
			#vxyz = quads.xyz[iquad][ivert]
			#d = distance(vxyz,p1xyz)
			#if d < 1.1*(pdiameter/2.0):
			#	notchamfer = False
			
			if notchamfer:
				avgvert = [0.0,0.0,0.0]
				cqv = quads.connect_quads_and_vertices[iquad][ivert] 
				nnbr = len(cqv)

				for nbr_quads_vertices in cqv:
					# simple average now, may use harmonic average later
					iquad2 = nbr_quads_vertices[0]
					ivert2 = nbr_quads_vertices[1]
					
					avgvert[0] = avgvert[0] + quads.newxyz[iquad2][ivert2][0]*(1.0/nnbr)
					avgvert[1] = avgvert[1] + quads.newxyz[iquad2][ivert2][1]*(1.0/nnbr)
					avgvert[2] = avgvert[2] + quads.newxyz[iquad2][ivert2][2]*(1.0/nnbr)
					
								
				for nbr_quads_vertices in cqv:
					# simple average now, may use harmonic average later
					iquad2 = nbr_quads_vertices[0]
					ivert2 = nbr_quads_vertices[1]
					
					quads.newxyz[iquad2][ivert2] = avgvert
			
			else:
				quads.newxyz[iquad][ivert] = quads.xyz[iquad][ivert]
	
	for iquad in range(0,nquads):
		for ivert in range(0,4):
			quads.xyz[iquad][ivert] = quads.newxyz[iquad][ivert]
		
		