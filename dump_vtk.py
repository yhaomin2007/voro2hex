# write hex mesh into vtk file.

def	dump_hex_vtk(hex20s,vtkFileName):
		
	vtkFileFile = open(vtkFileName,'w')
	header = '# vtk DataFile Version 2.0 \nUnstructured Grid Example \nASCII \nDATASET UNSTRUCTURED_GRID\n'
	
	vtkFileFile. write(header)
	nhexs = hex20s.nhex

	npoints = nhexs*8	
	line = 'POINTS '+str(npoints) + ' float\n'
	vtkFileFile. write(line)
	for ihex in range(0,nhexs):
		v8 = hex20s.v8[ihex]
		for ip in range (0,8):
			line = str(v8[ip][0]) + ' ' + str(v8[ip][1]) + ' ' +str(v8[ip][2]) + '\n' 
			vtkFileFile. write(line)

	print 'writing ' +str(nhexs)+ ' hexes to vtk fie'
			
	line = 'CELLS '+str(nhexs) + ' '+ str(nhexs*9)+'\n'
	vtkFileFile. write(line)
	ipt = 0
	for ihex in range(0,nhexs):
		line = '8 '
		for ip in range (0,8):
			line = line + str(ipt) + ' '
			ipt = ipt + 1
		line = line + '\n'
		vtkFileFile. write(line)


	line = 'CELL_TYPES  '+str(nhexs) +'\n'
	vtkFileFile. write(line)
	for ihex in range(0,nhexs):
		line = '12\n'
		vtkFileFile. write(line)

	vtkFileFile.close()

def	dump_hex_2nd_vtk(hex20s,vtkFileName):
		
	vtkFileFile = open(vtkFileName,'w')
	header = '# vtk DataFile Version 2.0 \nUnstructured Grid Example \nASCII \nDATASET UNSTRUCTURED_GRID\n'
	
	vtkFileFile. write(header)
	nhexs = hex20s.nhex

	npoints = nhexs*20	
	line = 'POINTS '+str(npoints) + ' float\n'
	vtkFileFile. write(line)
	for ihex in range(0,nhexs):
		v8 = hex20s.v8[ihex]
		for ip in range (0,8):
			line = str(v8[ip][0]) + ' ' + str(v8[ip][1]) + ' ' +str(v8[ip][2]) + '\n' 
			vtkFileFile. write(line)
		
		e12 = hex20s.e12[ihex]
		for ip in range (0,12):
			line = str(e12[ip][0]) + ' ' + str(e12[ip][1]) + ' ' +str(e12[ip][2]) + '\n' 
			vtkFileFile. write(line)

	print 'writing ' +str(nhexs)+ ' hexes to vtk fie'
			
	line = 'CELLS '+str(nhexs) + ' '+ str(nhexs*21)+'\n'
	vtkFileFile. write(line)
	ipt = 0
	for ihex in range(0,nhexs):
		line = '20 '
		for ip in range (0,20):
			line = line + str(ipt) + ' '
			ipt = ipt + 1
		line = line + '\n'
		vtkFileFile. write(line)


	line = 'CELL_TYPES  '+str(nhexs) +'\n'
	vtkFileFile. write(line)
	for ihex in range(0,nhexs):
		line = '25\n'
		vtkFileFile. write(line)

	vtkFileFile.close()
	
def	dump_hex_2nd_vtk_nnrh(hex20s,vtkFileName):
	# dump nek non right hand elements
	print 'calling: dump_hex_2nd_vtk_nnrh'
	
	vtkFileFile = open(vtkFileName,'w')
	header = '# vtk DataFile Version 2.0 \nUnstructured Grid Example \nASCII \nDATASET UNSTRUCTURED_GRID\n'
	
	vtkFileFile. write(header)
	nhexs = hex20s.nhex
	
	nnnrh = 0
	for ihex in range(0,nhexs):
		if hex20s.nnrh[ihex]:
			nnnrh = nnnrh + 1

	npoints = nnnrh*20	
	line = 'POINTS '+str(npoints) + ' float\n'
	vtkFileFile. write(line)
	for ihex in range(0,nhexs):
		if hex20s.nnrh[ihex]:
			v8 = hex20s.v8[ihex]
			for ip in range (0,8):
				line = str(v8[ip][0]) + ' ' + str(v8[ip][1]) + ' ' +str(v8[ip][2]) + '\n' 
				vtkFileFile. write(line)
		
			e12 = hex20s.e12[ihex]
			for ip in range (0,12):
				line = str(e12[ip][0]) + ' ' + str(e12[ip][1]) + ' ' +str(e12[ip][2]) + '\n' 
				vtkFileFile. write(line)

	print 'writing ' +str(nnnrh)+ ' hexes to vtk fie'
			
	line = 'CELLS '+str(nnnrh) + ' '+ str(nnnrh*21)+'\n'
	vtkFileFile. write(line)
	ipt = 0
	for ihex in range(0,nhexs):
		if hex20s.nnrh[ihex]:
			line = '20 '
			for ip in range (0,20):
				line = line + str(ipt) + ' '
				ipt = ipt + 1
			line = line + '\n'
			vtkFileFile. write(line)


	line = 'CELL_TYPES  '+str(nnnrh) +'\n'
	vtkFileFile. write(line)
	for ihex in range(0,nhexs):
		if hex20s.nnrh[ihex]:
			line = '25\n'
			vtkFileFile. write(line)

	vtkFileFile.close()
	
def	dump_hex_2nd_vtk_lar(hex20s,vtkFileName):
	# dump nek large aspec ratio elements
	
	print 'calling: dump_hex_2nd_vtk_lar'
	vtkFileFile = open(vtkFileName,'w')
	header = '# vtk DataFile Version 2.0 \nUnstructured Grid Example \nASCII \nDATASET UNSTRUCTURED_GRID\n'
	
	vtkFileFile. write(header)
	nhexs = hex20s.nhex
	
	nnnrh = 0
	for ihex in range(0,nhexs):
		if hex20s.lar[ihex]:
			nnnrh = nnnrh + 1

	npoints = nnnrh*20	
	line = 'POINTS '+str(npoints) + ' float\n'
	vtkFileFile. write(line)
	for ihex in range(0,nhexs):
		if hex20s.lar[ihex]:
			v8 = hex20s.v8[ihex]
			for ip in range (0,8):
				line = str(v8[ip][0]) + ' ' + str(v8[ip][1]) + ' ' +str(v8[ip][2]) + '\n' 
				vtkFileFile. write(line)
		
			e12 = hex20s.e12[ihex]
			for ip in range (0,12):
				line = str(e12[ip][0]) + ' ' + str(e12[ip][1]) + ' ' +str(e12[ip][2]) + '\n' 
				vtkFileFile. write(line)

	print 'writing ' +str(nnnrh)+ ' hexes to vtk fie'
			
	line = 'CELLS '+str(nnnrh) + ' '+ str(nnnrh*21)+'\n'
	vtkFileFile. write(line)
	ipt = 0
	for ihex in range(0,nhexs):
		if hex20s.lar[ihex]:
			line = '20 '
			for ip in range (0,20):
				line = line + str(ipt) + ' '
				ipt = ipt + 1
			line = line + '\n'
			vtkFileFile. write(line)


	line = 'CELL_TYPES  '+str(nnnrh) +'\n'
	vtkFileFile. write(line)
	for ihex in range(0,nhexs):
		if hex20s.lar[ihex]:
			line = '25\n'
			vtkFileFile. write(line)

	vtkFileFile.close()
	
def	dump_quad_vtk(quads,vtkFileName):
	
	
	vtkFileFile = open(vtkFileName,'w')
	header = '# vtk DataFile Version 2.0 \nUnstructured Grid Example \nASCII \nDATASET UNSTRUCTURED_GRID\n'
	
	vtkFileFile. write(header)
	nquads = quads.nquads

	npoints = nquads*4
	line = 'POINTS '+str(npoints) + ' float\n'
	vtkFileFile. write(line)

	for iquad in range(0,nquads):
		qxyz = quads.xyz[iquad]
		for ip in range (0,4):
			line = str(qxyz[ip][0]) + ' ' + str(qxyz[ip][1]) + ' ' +str(qxyz[ip][2]) + '\n' 
			vtkFileFile. write(line)
			
	print 'writing ' +str(nquads)+ ' quads to vtk fie'
			
	line = 'CELLS '+str(nquads) + ' '+ str(nquads*5)+'\n'
	vtkFileFile. write(line)
	ipt = 0
	for iquad in range(0,nquads):
		line = '4 '
		for ip in range (0,4):
			line = line + str(ipt) + ' '
			ipt = ipt + 1
		line = line + '\n'
		vtkFileFile. write(line)


	line = 'CELL_TYPES  '+str(nquads) +'\n'
	vtkFileFile. write(line)
	for iquad in range(0,nquads):
		line = '9\n'
		vtkFileFile. write(line)

	vtkFileFile.close()
	
def	dump_quad_2nd_vtk(quads,vtkFileName):
	
	
	vtkFileFile = open(vtkFileName,'w')
	header = '# vtk DataFile Version 2.0 \nUnstructured Grid Example \nASCII \nDATASET UNSTRUCTURED_GRID\n'
	
	vtkFileFile. write(header)
	nquads = quads.nquads

	npoints = nquads*8
	line = 'POINTS '+str(npoints) + ' float\n'
	vtkFileFile. write(line)

	for iquad in range(0,nquads):
		qxyz = quads.xyz[iquad]
		for ip in range (0,4):
			line = str(qxyz[ip][0]) + ' ' + str(qxyz[ip][1]) + ' ' +str(qxyz[ip][2]) + '\n' 
			vtkFileFile. write(line)
		qxyz = quads.mid_xyz[iquad]
		for ip in range (0,4):
			line = str(qxyz[ip][0]) + ' ' + str(qxyz[ip][1]) + ' ' +str(qxyz[ip][2]) + '\n' 
			vtkFileFile. write(line)
			
	print 'writing ' +str(nquads)+ ' quads to vtk fie'
			
	line = 'CELLS '+str(nquads) + ' '+ str(nquads*9)+'\n'
	vtkFileFile. write(line)
	ipt = 0
	for iquad in range(0,nquads):
		line = '8 '
		for ip in range (0,8):
			line = line + str(ipt) + ' '
			ipt = ipt + 1
		line = line + '\n'
		vtkFileFile. write(line)


	line = 'CELL_TYPES  '+str(nquads) +'\n'
	vtkFileFile. write(line)
	for iquad in range(0,nquads):
		line = '23\n'
		vtkFileFile. write(line)

	vtkFileFile.close()
	
def	dump_poly_vtk(vor_vert,vor_face,vtkFileName):
	
	
	vtkFileFile = open(vtkFileName,'w')
	header = '# vtk DataFile Version 2.0 \nUnstructured Grid Example \nASCII \nDATASET UNSTRUCTURED_GRID\n'
	
	vtkFileFile. write(header)
	nv = vor_vert.nv
	nf = vor_face.nf
	
	nv_from_face = 0
	for iface in range(0,nf):
		ifghost = vor_face.ifghost[iface]
		ifcollapse = vor_face.ifcollapse[iface]
		if (not ifghost) and (not ifcollapse):
			nvf = len(vor_face.f_to_v[iface])
			nv_from_face = nv_from_face + nvf
	
	line = 'POINTS '+str(nv_from_face) + ' float\n'
	vtkFileFile. write(line)

	for iface in range(0,nf):
		ifghost = vor_face.ifghost[iface]
		ifcollapse = vor_face.ifcollapse[iface]
		if (not ifghost) and (not ifcollapse):
			for ivert in vor_face.f_to_v[iface]:
				vxyz = vor_vert.xyz[ivert]
				line = str(vxyz[0]) + ' ' + str(vxyz[1]) + ' ' +str(vxyz[2]) + '\n' 
				vtkFileFile. write(line)
			
	celllistsize = 0
	nf_r = 0
	for iface in range(0,nf):
		ifghost = vor_face.ifghost[iface]
		ifcollapse = vor_face.ifcollapse[iface]
		if (not ifghost) and (not ifcollapse):
			nvf = len(vor_face.f_to_v[iface])
			nf_r = nf_r + 1
			celllistsize = celllistsize + nvf + 1
			
	print 'writing ' +str(nf_r)+ ' polygons to vtk fie'
			
	line = 'CELLS '+str(nf_r) + ' '+ str(celllistsize)+'\n'
	ivf = 0
	vtkFileFile.write(line)
	for iface in range(0,nf):
		ifghost = vor_face.ifghost[iface]
		ifcollapse = vor_face.ifcollapse[iface]
		if (not ifghost) and (not ifcollapse):
			nvf = len(vor_face.f_to_v[iface])
			line = str(nvf)+' '
			for ip in range (0,nvf):
				#ivf = vor_face.f_to_v[iface][ip]
				line = line + str(ivf) + ' '
				ivf = ivf + 1
			line = line + '\n'
			vtkFileFile.write(line)


	line = 'CELL_TYPES  '+str(nf_r) +'\n'
	vtkFileFile.write(line)
	for iquad in range(0,nf_r):
		line = '7\n'
		vtkFileFile.write(line)

	vtkFileFile.close()
	print 'done:: writing ' +str(nf_r)+ ' polygons to vtk fie'
	
def	dump_poly_vtk_for_chamfer(vor_vert,vor_face,vtkFileName):
	
	
	vtkFileFile = open(vtkFileName,'w')
	header = '# vtk DataFile Version 2.0 \nUnstructured Grid Example \nASCII \nDATASET UNSTRUCTURED_GRID\n'
	
	vtkFileFile. write(header)
	nv = vor_vert.nv
	nf = vor_face.nf
	
	nv_from_face = 0
	for iface in range(0,nf):
		ifghost = vor_face.ifghost[iface]
		ifcollapse = vor_face.ifcollapse[iface]
		ifchamfer = vor_face.ifchamfer[iface] 
		if (not ifghost) and (not ifcollapse) and ifchamfer:
			nvf = len(vor_face.f_to_v[iface])
			nv_from_face = nv_from_face + nvf
	
	line = 'POINTS '+str(nv_from_face) + ' float\n'
	vtkFileFile. write(line)

	for iface in range(0,nf):
		ifghost = vor_face.ifghost[iface]
		ifcollapse = vor_face.ifcollapse[iface]
		ifchamfer = vor_face.ifchamfer[iface] 
		if (not ifghost) and (not ifcollapse) and ifchamfer:
			for ivert in vor_face.f_to_v[iface]:
				vxyz = vor_vert.xyz[ivert]
				line = str(vxyz[0]) + ' ' + str(vxyz[1]) + ' ' +str(vxyz[2]) + '\n' 
				vtkFileFile. write(line)
			
	celllistsize = 0
	nf_r = 0
	for iface in range(0,nf):
		ifghost = vor_face.ifghost[iface]
		ifcollapse = vor_face.ifcollapse[iface]
		ifchamfer = vor_face.ifchamfer[iface] 
		if (not ifghost) and (not ifcollapse) and ifchamfer:
			nvf = len(vor_face.f_to_v[iface])
			nf_r = nf_r + 1
			celllistsize = celllistsize + nvf + 1
			
	print 'writing ' +str(nf_r)+ ' polygons to vtk fie'
			
	line = 'CELLS '+str(nf_r) + ' '+ str(celllistsize)+'\n'
	ivf = 0
	vtkFileFile. write(line)
	for iface in range(0,nf):
		ifghost = vor_face.ifghost[iface]
		ifcollapse = vor_face.ifcollapse[iface]
		ifchamfer = vor_face.ifchamfer[iface] 
		if (not ifghost) and (not ifcollapse) and ifchamfer:
			nvf = len(vor_face.f_to_v[iface])
			line = str(nvf)+' '
			for ip in range (0,nvf):
				#ivf = vor_face.f_to_v[iface][ip]
				line = line + str(ivf) + ' '
				ivf = ivf + 1
			line = line + '\n'
			vtkFileFile. write(line)


	line = 'CELL_TYPES  '+str(nf_r) +'\n'
	vtkFileFile.write(line)
	for iquad in range(0,nf_r):
		line = '7\n'
		vtkFileFile.write(line)

	vtkFileFile.close()