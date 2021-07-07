from data_class import *

def write_rea(reaFile,newReaFile,hex20s):
	# reaFile: 		base rea file name
	# newReaFile: 	new rea file name
	# hex20s: 		nek_hex20s class data
	
	print '==========================='
	print 'start writing new .rea file'
	
	reaFileHolder = open(reaFile,'rw')
	wholeReaFile = reaFileHolder.read()
	
	newReaFileHolder = open(newReaFile,'w')
	newReaFileHolder.seek(0)
	
	newRea = wholeReaFile[:wholeReaFile.find('**MESH')-1]
	
	newRea = newRea + ' **MESH DATA** 6 lines are X,Y,Z;X,Y,Z. Columns corners 1-4;5-8\n'

	nhexs = hex20s.nhex
	print 'for '+str(nhexs) + ' elements'
	
	newRea = newRea + '         '+str(nhexs)+'  3         '+str(nhexs)+'           NEL,NDIM,NELV\n'

	newReaFileHolder.write(newRea)
	newRea = ''

	rjust1 = 10 
	rjust2 = 14
	
	scale = 1.0
	for i in range(0,nhexs):
		v8 = hex20s.v8[i]
		#print 'write element ',i
		# write cells information into rea file. each cell is regarded as an element in nek
		
		newRea = newRea + '            ELEMENT     '+str(i+1).rjust(7)+' [    1a]  GROUP  0\n'
		# adding xyz for points 1-4
		# adding x
		newRea = (newRea + format(v8[0][0]*scale,'f').rjust(rjust1) + format(v8[1][0]*scale,'f').rjust(rjust2)
			+ format(v8[2][0]*scale,'f').rjust(rjust2) + format(v8[3][0]*scale,'f').rjust(rjust2) +'\n')
		# adding y
		newRea = (newRea + format(v8[0][1]*scale,'f').rjust(rjust1) + format(v8[1][1]*scale,'f').rjust(rjust2)
			+ format(v8[2][1]*scale,'f').rjust(rjust2) + format(v8[3][1]*scale,'f').rjust(rjust2) +'\n')
		# adding z
		newRea =  (newRea + format(v8[0][2]*scale,'f').rjust(rjust1) + format(v8[1][2]*scale,'f').rjust(rjust2)
			+ format(v8[2][2]*scale,'f').rjust(rjust2) + format(v8[3][2]*scale,'f').rjust(rjust2) +'\n')
		
		# adding xyz for points 5-8
		# adding x
		newRea = (newRea + format(v8[4][0]*scale,'f').rjust(rjust1) + format(v8[5][0]*scale,'f').rjust(rjust2)
			+ format(v8[6][0]*scale,'f').rjust(rjust2) + format(v8[7][0]*scale,'f').rjust(rjust2) +'\n')
		# adding y
		newRea = (newRea + format(v8[4][1]*scale,'f').rjust(rjust1) + format(v8[5][1]*scale,'f').rjust(rjust2)
			+ format(v8[6][1]*scale,'f').rjust(rjust2) + format(v8[7][1]*scale,'f').rjust(rjust2) +'\n')
		# adding z
		newRea =  (newRea + format(v8[4][2]*scale,'f').rjust(rjust1) + format(v8[5][2]*scale,'f').rjust(rjust2)
			+ format(v8[6][2]*scale,'f').rjust(rjust2) + format(v8[7][2]*scale,'f').rjust(rjust2) +'\n')
		
		# write to file
		newReaFileHolder.write(newRea)
		newRea = ''
	# all elements are written.
	
	newRea =''
	newRea = newRea + '  ***** CURVED SIDE DATA *****\n'
	newReaFileHolder.write(newRea)
	
	ifcurve = True
	
	if ifcurve:  # not quite work, nek cannot read it correctly
		# write curvature infos
		ncurves = nhexs*12
		newRea = ''
		newRea = newRea + '      '+str(ncurves)+' Curved sides follow IEDGE,IEL,CURVE(I),I=1,5, CCURVE\n'
		newReaFileHolder.write(newRea)
	
		
		# this is to match exactly rea file format......very strange...
		rjust3 = 14
		
		if nhexs < 1000:
			rjust1 = 3
			rjust2 = 3
		elif nhexs < 1e6:
			rjust1 = 2
			rjust2 = 6
		else:
			rjust1 = 2
			rjust2 = 12
		#
		
		for ih in range(0,nhexs):
			e12 = hex20s.e12[ih]
			
			for ie in range(0,12):
				newRea = str(ie+1).rjust(rjust1)+str(ih+1).rjust(rjust2) + format(e12[ie][0],'f').rjust(rjust3)+ format(e12[ie][1],'f').rjust(rjust3)+ format(e12[ie][2],'f').rjust(rjust3)
				newRea = newRea + format(0.0,'f').rjust(rjust3)+format(0.0,'f').rjust(rjust3)+' m\n'
				newReaFileHolder.write(newRea)
	
	else:
		newRea = ''
		newRea = newRea + '      0 Curved sides follow IEDGE,IEL,CURVE(I),I=1,5, CCURVE\n'
		newReaFileHolder.write(newRea)

	newRea = ''
	newRea = newRea + ' ***** BOUNDARY CONDITIONS *****\n  ***** FLUID   BOUNDARY CONDITIONS *****\n'
	newReaFileHolder.write(newRea)
	newRea = ''
	
	# write BC tag to element faces
	# no element connectivity information
	for i in range(0,nhexs):
		#print 'write faces for element ',i
		s6 = hex20s.s6[i]
		for j in range(0,6):
		
			#bctag = 'E  '
			#bctag = 'TOP'
			#bctag = 'BOT'
			#bctag = 'SW '
			#bctag = 'PW '
			bctag = s6[j]
			newRea =' ' + bctag
			
			k = 5
			if ((i+1) > 999): k = 6
			if ((i+1) > 9999): k = 7
			if ((i+1) > 99999): k = 8
			if ((i+1) > 999999): k = 9
			if ((i+1) > 9999999): k = 10
			if ((i+1) > 99999999): k = 11
			#if (len(cells) > 100000):
			#	print 'WARNING: too many elements (more than 1E6)'
			#	sys.exit()
			
			newRea = newRea + str(i+1).rjust(k) + str(j+1).rjust(3)
			newRea = newRea + '\n'
			newReaFileHolder.write(newRea)
			newRea = ''
	
	newRea = newRea + '  ***** NO THERMAL BOUNDARY CONDITIONS *****\n'
	
	newRea = newRea + wholeReaFile[wholeReaFile.find('0 PRESOLVE/RESTART OPTIONS  *****')-10:]
	
	newReaFileHolder.write(newRea)
	newRea = ''
	
	reaFileHolder.close()
	newReaFileHolder.close()

	print 'new .rea file done'

def write_rea_cht(reaFile,newReaFile,hex20s,solid_hex20s):
	# reaFile: 		base rea file name
	# newReaFile: 	new rea file name
	# hex20s: 		fluid nek_hex20s class data
	# solid_hex20s: solid nek_hex20s class data
	
	print '==========================='
	print 'start writing new .rea file'
	
	reaFileHolder = open(reaFile,'rw')
	wholeReaFile = reaFileHolder.read()
	
	newReaFileHolder = open(newReaFile,'w')
	newReaFileHolder.seek(0)
	
	newRea = wholeReaFile[:wholeReaFile.find('**MESH')-1]
	
	newRea = newRea + ' **MESH DATA** 6 lines are X,Y,Z;X,Y,Z. Columns corners 1-4;5-8\n'

	fluid_nhexs = hex20s.nhex
	print 'for '+str(fluid_nhexs) + ' fluid elements'
	
	solid_nhexs = solid_hex20s.nhex
	
	print 'and for '+str(solid_nhexs) + ' solid elements'
	
	nhexs = fluid_nhexs + solid_nhexs
	
	newRea = newRea + '         '+str(nhexs)+'  3         '+str(fluid_nhexs)+'           NEL,NDIM,NELV\n'

	newReaFileHolder.write(newRea)
	newRea = ''

	rjust1 = 10 
	rjust2 = 14
	
	scale = 1.0
	# write fluid hex8 
	for i in range(0,fluid_nhexs):
		v8 = hex20s.v8[i]
		#print 'write element ',i
		# write cells information into rea file. each cell is regarded as an element in nek
		
		newRea = newRea + '            ELEMENT     '+str(i+1).rjust(7)+' [    1a]  GROUP  0\n'
		# adding xyz for points 1-4
		# adding x
		newRea = (newRea + format(v8[0][0]*scale,'f').rjust(rjust1) + format(v8[1][0]*scale,'f').rjust(rjust2)
			+ format(v8[2][0]*scale,'f').rjust(rjust2) + format(v8[3][0]*scale,'f').rjust(rjust2) +'\n')
		# adding y
		newRea = (newRea + format(v8[0][1]*scale,'f').rjust(rjust1) + format(v8[1][1]*scale,'f').rjust(rjust2)
			+ format(v8[2][1]*scale,'f').rjust(rjust2) + format(v8[3][1]*scale,'f').rjust(rjust2) +'\n')
		# adding z
		newRea =  (newRea + format(v8[0][2]*scale,'f').rjust(rjust1) + format(v8[1][2]*scale,'f').rjust(rjust2)
			+ format(v8[2][2]*scale,'f').rjust(rjust2) + format(v8[3][2]*scale,'f').rjust(rjust2) +'\n')
		
		# adding xyz for points 5-8
		# adding x
		newRea = (newRea + format(v8[4][0]*scale,'f').rjust(rjust1) + format(v8[5][0]*scale,'f').rjust(rjust2)
			+ format(v8[6][0]*scale,'f').rjust(rjust2) + format(v8[7][0]*scale,'f').rjust(rjust2) +'\n')
		# adding y
		newRea = (newRea + format(v8[4][1]*scale,'f').rjust(rjust1) + format(v8[5][1]*scale,'f').rjust(rjust2)
			+ format(v8[6][1]*scale,'f').rjust(rjust2) + format(v8[7][1]*scale,'f').rjust(rjust2) +'\n')
		# adding z
		newRea =  (newRea + format(v8[4][2]*scale,'f').rjust(rjust1) + format(v8[5][2]*scale,'f').rjust(rjust2)
			+ format(v8[6][2]*scale,'f').rjust(rjust2) + format(v8[7][2]*scale,'f').rjust(rjust2) +'\n')
		
		# write to file
		newReaFileHolder.write(newRea)
		newRea = ''
		
	for i in range(0,solid_nhexs):
		v8 = solid_hex20s.v8[i]
		#print 'write solid element ',i
		# write cells information into rea file. each cell is regarded as an element in nek
		
		newRea = newRea + '            ELEMENT     '+str(i+1+fluid_nhexs).rjust(7)+' [    1a]  GROUP  0\n'
		# adding xyz for points 1-4
		# adding x
		newRea = (newRea + format(v8[0][0]*scale,'f').rjust(rjust1) + format(v8[1][0]*scale,'f').rjust(rjust2)
			+ format(v8[2][0]*scale,'f').rjust(rjust2) + format(v8[3][0]*scale,'f').rjust(rjust2) +'\n')
		# adding y
		newRea = (newRea + format(v8[0][1]*scale,'f').rjust(rjust1) + format(v8[1][1]*scale,'f').rjust(rjust2)
			+ format(v8[2][1]*scale,'f').rjust(rjust2) + format(v8[3][1]*scale,'f').rjust(rjust2) +'\n')
		# adding z
		newRea =  (newRea + format(v8[0][2]*scale,'f').rjust(rjust1) + format(v8[1][2]*scale,'f').rjust(rjust2)
			+ format(v8[2][2]*scale,'f').rjust(rjust2) + format(v8[3][2]*scale,'f').rjust(rjust2) +'\n')
		
		# adding xyz for points 5-8
		# adding x
		newRea = (newRea + format(v8[4][0]*scale,'f').rjust(rjust1) + format(v8[5][0]*scale,'f').rjust(rjust2)
			+ format(v8[6][0]*scale,'f').rjust(rjust2) + format(v8[7][0]*scale,'f').rjust(rjust2) +'\n')
		# adding y
		newRea = (newRea + format(v8[4][1]*scale,'f').rjust(rjust1) + format(v8[5][1]*scale,'f').rjust(rjust2)
			+ format(v8[6][1]*scale,'f').rjust(rjust2) + format(v8[7][1]*scale,'f').rjust(rjust2) +'\n')
		# adding z
		newRea =  (newRea + format(v8[4][2]*scale,'f').rjust(rjust1) + format(v8[5][2]*scale,'f').rjust(rjust2)
			+ format(v8[6][2]*scale,'f').rjust(rjust2) + format(v8[7][2]*scale,'f').rjust(rjust2) +'\n')
		
		# write to file
		newReaFileHolder.write(newRea)
		newRea = ''
		
	# all elements are written.
	
	newRea =''
	newRea = newRea + '  ***** CURVED SIDE DATA *****\n'
	newReaFileHolder.write(newRea)
	
	ifcurve = True
	
	if ifcurve:  # not quite work, nek cannot read it correctly
		# write curvature infos
		ncurves = nhexs*12
		newRea = ''
		newRea = newRea + '      '+str(ncurves)+' Curved sides follow IEDGE,IEL,CURVE(I),I=1,5, CCURVE\n'
		newReaFileHolder.write(newRea)
	
		
		# this is to match exactly rea file format......very strange...
		rjust3 = 14
		
		if nhexs < 1000:
			rjust1 = 3
			rjust2 = 3
		elif nhexs < 1e6:
			rjust1 = 2
			rjust2 = 6
		else:
			rjust1 = 2
			rjust2 = 12
		#
		# fluid element curvature
		for ih in range(0,fluid_nhexs):
			e12 = hex20s.e12[ih]
			
			for ie in range(0,12):
				newRea = str(ie+1).rjust(rjust1)+str(ih+1).rjust(rjust2) + format(e12[ie][0],'f').rjust(rjust3)+ format(e12[ie][1],'f').rjust(rjust3)+ format(e12[ie][2],'f').rjust(rjust3)
				newRea = newRea + format(0.0,'f').rjust(rjust3)+format(0.0,'f').rjust(rjust3)+' m\n'
				newReaFileHolder.write(newRea)
		
		# solid element curvature
		for ih in range(0,solid_nhexs):
			e12 = solid_hex20s.e12[ih]
			
			for ie in range(0,12):
				newRea = str(ie+1).rjust(rjust1)+str(ih+1+fluid_nhexs).rjust(rjust2) + format(e12[ie][0],'f').rjust(rjust3)+ format(e12[ie][1],'f').rjust(rjust3)+ format(e12[ie][2],'f').rjust(rjust3)
				newRea = newRea + format(0.0,'f').rjust(rjust3)+format(0.0,'f').rjust(rjust3)+' m\n'
				newReaFileHolder.write(newRea)
	
	else:
		newRea = ''
		newRea = newRea + '      0 Curved sides follow IEDGE,IEL,CURVE(I),I=1,5, CCURVE\n'
		newReaFileHolder.write(newRea)

	newRea = ''
	newRea = newRea + ' ***** BOUNDARY CONDITIONS *****\n  ***** FLUID   BOUNDARY CONDITIONS *****\n'
	newReaFileHolder.write(newRea)
	newRea = ''
	
	# write BC tag to element faces
	# no element connectivity information

	# fluid bc
	for i in range(0,fluid_nhexs):
		#print 'write faces for element ',i
		s6 = hex20s.s6[i]
		for j in range(0,6):
		
			#bctag = 'E  '
			#bctag = 'TOP'
			#bctag = 'BOT'
			#bctag = 'SW '
			#bctag = 'PW '
			bctag = s6[j]
			newRea =' ' + bctag
			
			k = 5
			if ((i+1) > 999): k = 6
			if ((i+1) > 9999): k = 7
			if ((i+1) > 99999): k = 8
			if ((i+1) > 999999): k = 9
			if ((i+1) > 9999999): k = 10
			if ((i+1) > 99999999): k = 11
			#if (len(cells) > 100000):
			#	print 'WARNING: too many elements (more than 1E6)'
			#	sys.exit()
			
			newRea = newRea + str(i+1).rjust(k) + str(j+1).rjust(3)
			newRea = newRea + '\n'
			newReaFileHolder.write(newRea)
			newRea = ''
	
	newRea = ''
	newRea = newRea + '  ***** THERMAL BOUNDARY CONDITIONS *****\n'
	newReaFileHolder.write(newRea)
	newRea = ''
	
	# fluid element thermal bc
	for i in range(0,fluid_nhexs):
		#print 'write faces for element ',i
		s6 = hex20s.s6[i]
		for j in range(0,6):
		
			#bctag = 'E  '
			#bctag = 'TOP'
			#bctag = 'BOT'
			#bctag = 'SW '
			#bctag = 'PW '
			bctag = s6[j]
			newRea =' ' + bctag
			
			k = 5
			if ((i+1) > 999): k = 6
			if ((i+1) > 9999): k = 7
			if ((i+1) > 99999): k = 8
			if ((i+1) > 999999): k = 9
			if ((i+1) > 9999999): k = 10
			if ((i+1) > 99999999): k = 11
			#if (len(cells) > 100000):
			#	print 'WARNING: too many elements (more than 1E6)'
			#	sys.exit()
			
			newRea = newRea + str(i+1).rjust(k) + str(j+1).rjust(3)
			newRea = newRea + '\n'
			newReaFileHolder.write(newRea)
			newRea = ''
	
	
	# solid bc
	for i in range(0,solid_nhexs):
		#print 'write faces for element ',i
		s6 = solid_hex20s.s6[i]
		for j in range(0,6):
		
			#bctag = 'E  '
			#bctag = 'TOP'
			#bctag = 'BOT'
			#bctag = 'SW '
			#bctag = 'PW '
			bctag = s6[j]
			newRea =' ' + bctag
			
			k = 5
			if ((i+1+fluid_nhexs) > 999): k = 6
			if ((i+1+fluid_nhexs) > 9999): k = 7
			if ((i+1+fluid_nhexs) > 99999): k = 8
			if ((i+1+fluid_nhexs) > 999999): k = 9
			if ((i+1+fluid_nhexs) > 9999999): k = 10
			if ((i+1+fluid_nhexs) > 99999999): k = 11
			#if (len(cells) > 100000):
			#	print 'WARNING: too many elements (more than 1E6)'
			#	sys.exit()
			
			newRea = newRea + str(i+1+fluid_nhexs).rjust(k) + str(j+1).rjust(3)
			newRea = newRea + '\n'
			newReaFileHolder.write(newRea)
			newRea = ''
	
	newRea = newRea + wholeReaFile[wholeReaFile.find('0 PRESOLVE/RESTART OPTIONS  *****')-10:]
	
	newReaFileHolder.write(newRea)
	newRea = ''
	
	reaFileHolder.close()
	newReaFileHolder.close()

	print 'new .rea file done'