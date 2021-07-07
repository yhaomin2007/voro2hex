# modifed version to reduce memory print
# so you can run on other machine other than compute001
#

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
#=================================================================================
starttime = time.time()

nprocs = 16 # number of parallel processors

# negative jacobian elements, need to linearrize here.
nje = [
] 
# low scale jacobian elements
lje = [
]

nje.extend(lje)

if len(nje)>0:
	nje_new = []
	for e in nje:
		nje_new.append(e)
		nje_new.append(e-1)
		nje_new.append(e-2)
	nje = remove_duplicates(nje_new)

enje = []
v8nje = []

reaFile = 'base.rea'
newReaFile = 'pba.rea'

print '==========================='
print 'start writing new .rea file'
	
reaFileHolder = open(reaFile,'rw')
wholeReaFile = reaFileHolder.read()

newReaFileHolder = open(newReaFile,'w')
newReaFileHolder.seek(0)	
newRea = wholeReaFile[:wholeReaFile.find('**MESH')-1]
newRea = newRea + ' **MESH DATA** 6 lines are X,Y,Z;X,Y,Z. Columns corners 1-4;5-8\n'
# get nhexs
tothexs = 0
for ip in range(0,nprocs):
	file =  'proc'+str(ip)+'/hfile.dat'
	print 'reading ' + file
	fileHolder = open(file,'r')
	line = fileHolder.readline()
	nhexs = int(line.split()[0])
	tothexs = tothexs + nhexs
	fileHolder.close()
	
	print 'reading ' + str(nhexs)+ ' hex elements'
	print 'total ' + str(tothexs)+ ' hex elements'

newRea = newRea + '         '+str(tothexs)+'  3         '+str(tothexs)+'           NEL,NDIM,NELV\n'
newReaFileHolder.write(newRea)
newRea = ''		
	
rjust1 = 10 
rjust2 = 14
	
scale = 1.0

th = 0
print 'writing v8 info'
# write v8
for ip in range(0,nprocs):
	file =  'proc'+str(ip)+'/hfile.dat'
	print 'reading ' + file
	fileHolder = open(file,'r')
	line = fileHolder.readline()
	nhexs = int(line.split()[0])
	
	print 'reading ' + str(nhexs)+ ' hex elements'
	
	for ih in range(0,nhexs):
		
		v8 = []
		# read v8
		line = fileHolder.readline()
		for iv in range(0,8):
			vxyz = []
			for i in range(0,3):
				index = iv*3+i
				vxyz.append(float(line.split()[index]))
			v8.append(vxyz)
		
		th = th + 1
		if th in nje:
			print 'store nje ' + str(th) + ' v8 info'
			enje.append(th)
			v8nje.append(v8)
		
		newRea = newRea + '            ELEMENT     '+str(th).rjust(7)+' [    1a]  GROUP  0\n'
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
		
		line = fileHolder.readline()
		line = fileHolder.readline()
		
	fileHolder.close()

	
newRea =''
newRea = newRea + '  ***** CURVED SIDE DATA *****\n'
newReaFileHolder.write(newRea)	
		
ncurves = tothexs*12
newRea = ''
newRea = newRea + '      '+str(ncurves)+' Curved sides follow IEDGE,IEL,CURVE(I),I=1,5, CCURVE\n'
newReaFileHolder.write(newRea)

rjust3 = 14
		
if tothexs < 1000:
	rjust1 = 3
	rjust2 = 3
elif tothexs < 1e6:
	rjust1 = 2
	rjust2 = 6
else:
	rjust1 = 2
	rjust2 = 12


print 'writing e12 info'

th = 0
	
for ip in range(0,nprocs):
	file =  'proc'+str(ip)+'/hfile.dat'
	print 'reading ' + file
	fileHolder = open(file,'r')
	line = fileHolder.readline()
	nhexs = int(line.split()[0])
	
	for ih in range(0,nhexs):
		
		line = fileHolder.readline()
		
		th = th + 1
		e12 = []
		# read e12
		line = fileHolder.readline()
		for ie in range(0,12):
			vxyz = []
			for i in range(0,3):
				index = ie*3+i
				vxyz.append(float(line.split()[index]))
			e12.append(vxyz)
	
		if th in nje:
			print 'linearize nje ' + str(th) + ' element'
			ind1 = enje.index(th)
			v8 = v8nje[ind1]
			e12 = v8_to_e12(v8) # linearize this element
	
		for ie in range(0,12):
			newRea = str(ie+1).rjust(rjust1)+str(th).rjust(rjust2) + format(e12[ie][0],'f').rjust(rjust3)+ format(e12[ie][1],'f').rjust(rjust3)+ format(e12[ie][2],'f').rjust(rjust3)
			newRea = newRea + format(0.0,'f').rjust(rjust3)+format(0.0,'f').rjust(rjust3)+' m\n'
			newReaFileHolder.write(newRea)
	
		line = fileHolder.readline()
	fileHolder.close()

newRea = ''
newRea = newRea + ' ***** BOUNDARY CONDITIONS *****\n  ***** FLUID   BOUNDARY CONDITIONS *****\n'
newReaFileHolder.write(newRea)
newRea = ''


print 'writing s6 info'
i = -1
for ip in range(0,nprocs):
	file =  'proc'+str(ip)+'/hfile.dat'
	print 'reading ' + file
	fileHolder = open(file,'r')
	line = fileHolder.readline()
	nhexs = int(line.split()[0])
	
	for ih in range(0,nhexs):
		
		line = fileHolder.readline()
		line = fileHolder.readline()
		s6 = []
		line = fileHolder.readline()
		for iface in range(0,6):
			ss = int(line.split()[iface])
			if (ss == 0): s6.append('E  ')
			if (ss == 1): s6.append('PW ')
			if (ss == 2): s6.append('C  ')
			if (ss == 3): s6.append('SW ')
			if (ss == 4): s6.append('TOP')
			if (ss == 5): s6.append('BOT')

		i = i + 1
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
	fileHolder.close()
			
newRea = newRea + '  ***** NO THERMAL BOUNDARY CONDITIONS *****\n'
	
newRea = newRea + wholeReaFile[wholeReaFile.find('0 PRESOLVE/RESTART OPTIONS  *****')-10:]

newReaFileHolder.write(newRea)
newRea = ''
	
reaFileHolder.close()
newReaFileHolder.close()

print 'new .rea file done'	

endtime = time.time()
print 'time used:' + str(endtime-starttime)	
			