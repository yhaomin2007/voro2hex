import math
import multiprocessing
import os


def fix_non_right_hand_parallel(hex20,np):
	print 'fix non_right_hand_elements of hex20 elements in parallel'

	nhex = hex20.nhex	
	jobs = []
	start = []
	end = []
	
	nrange = int(math.ceil(nhex/np))

	for ip in range(0,np):
		start.append(ip*nrange)
		
	for ip in range(0,np-1):
		end.append(start[ip+1])	
	end.append(nhex)
	
	for ip in range(0,np):
		process = multiprocessing.Process(target=hex20.fix_non_right_hand_elements_range,args=(start[ip],end[ip]))
		jobs.append(process)
	
	for j in jobs:
		j.start()

	for j in jobs:
		j.join()


def check_max_aspect_ratio_parallel(hex20,np,ftag):
	print 'check aspect_ratio of hex20 elements in parallel'

	nhex = hex20.nhex	
	jobs = []
	start = []
	end = []
	
	for ih in range(0,hex20.nhex):
		hex20.nnrh[ih] = False
	
	nrange = int(math.ceil(nhex/np))

	for ip in range(0,np):
		start.append(ip*nrange)
		
	for ip in range(0,np-1):
		end.append(start[ip+1])	
	end.append(nhex)
	
	
	os.system('rm -rf dummy_'+ftag)
	os.system('mkdir dummy_'+ftag)
	
	for ip in range(0,np):
		process = multiprocessing.Process(target=hex20.check_max_aspect_ratio_range,args=(start[ip],end[ip],ftag,ip))
		jobs.append(process)
	
	for j in jobs:
		j.start()

	for j in jobs:
		j.join()
		
	maxar = 0.0
	
	print 'combining aspect-ratio information'
	for ip in range(0,np):
		filename = 'dummy_'+ftag+'/dummy_'+str(ip)
		dummyfile = open(filename)
		line = dummyfile.readline()
		nlines =  int(line.split()[0])
		for iline in range(0,nlines):
			line = dummyfile.readline()
			ih = int(line.split()[0])
			ar = float(line.split()[1])
			if ar > maxar: maxar = ar
			
			hex20.nnrh[ih] = True # borrow this flag to dump hex 
			dummyfile.close()
	print 'maxar: ' + str(maxar)
	print 'done combining all aspect-ratio info'	

	
def check_non_right_hand_elements_use_nek_method_parallel(hex20,np,ftag):
	print 'check non_right_hand_elements of hex20 elements in parallel'

	nhex = hex20.nhex	
	jobs = []
	start = []
	end = []
	
	for ih in range(0,hex20.nhex):
		hex20.nnrh[ih] = False
	
	nrange = int(math.ceil(nhex/np))

	for ip in range(0,np):
		start.append(ip*nrange)
		
	for ip in range(0,np-1):
		end.append(start[ip+1])	
	end.append(nhex)
	
	os.system('rm -rf dummy_'+ftag)
	os.system('mkdir dummy_'+ftag)
	
	
	for ip in range(0,np):
		process = multiprocessing.Process(target=hex20.check_non_right_hand_elements_use_nek_method_range,args=(start[ip],end[ip],ftag,ip))
		jobs.append(process)
	
	for j in jobs:
		j.start()

	for j in jobs:
		j.join()
		
	print 'combining non-hand-elements information'
	for ip in range(0,np):
		filename = 'dummy_'+ftag+'/dummy_'+str(ip)
		dummyfile = open(filename)
		line = dummyfile.readline()
		nlines =  int(line.split()[0])
		for iline in range(0,nlines):
			line = dummyfile.readline()
			ih = int(line.split()[0])
			hex20.nnrh[ih] = True # borrow this flag to dump hex 
		dummyfile.close()
	print 'done combining all non-hand-elements info'	
