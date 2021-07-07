
from splitting_subroutines import *
from merge_vertices import *
from data_class import *

def fix_negative_jacobian_elements_from_nek(hex20):
	# element list with -jac from nek
	#neg_jacobian_elist = [456516,303598,500166,56102,192156,358406,442204,396942,499416,156976,45706,96710,447622]
	#neg_jacobian_elist = [353362,332940,95370,154784,299612,49238,45188,391664,442308,40238,449694,244830,441220,492152]
	#neg_jacobian_elist = [53684,53685,375806,42419,42420,318540,49112,483717,483718,354335,418757,418758,476144,100631,163626,474880,474881,258954,528317,258955]
	
	#neg_jacobian_elist = [491601,491602,485701,485702,55957,232723,55958,329498,365503,142017,433023,433024,104485,501245,451230,501246]

	#neg_jacobian_elist = [2053695,2053693,2185513]
	
	#neg_jacobian_elist = [1884425,1884426,1888967,2112314,2152169,1888968,2112315,2112326,2112327,1901321,2152170,1887830,1901333,1901334,1901322,1887831,1887914,1887915,2152418,2152419,2152436,2152437,2152439,2152440,1980266,1980267,1888349,1967228,2044124,
	#1888350,1983260,1967229,1967237,1967238,1921667,1986845,2152541,
	#1986695,1986695,2101154,1901897,2171538,1968869,1921668,1934705,
	#1934706,1940508,1886030,1889039,2101155,1892864,1901898,1901909,
	#1901910,2105663,2171537,2173352,2173353,2173355,2044125,1983261,
	#1934708,1934709,2157677,1886031,1892865,1892876,2105664,1968870,1968881,1968882,2156075,1986846,2157678,1889040,2099966,1892877,
	#2173356,2156076,1986696,1889051,2099967,2099978,2099979,1996592,
	#2173334,2152542,1889052,2173335,1996593,2163689,2163690]
	
	#neg_jacobian_elist = [2096958,2096959,2100556,2100557]
	
	#neg_jacobian_elist = [2059258,2059259,2059260,2112361,2112362,2112363]
	
	#neg_jacobian_elist = [2059258,2059259,2059260,2112361,2112362,2112363,2057884,2057885,2057886,2110723,2110724,2110725]
	
	#neg_jacobian_elist = [2212907,2212908,2107723,2107724,2107725,2055430,2055431,2055432]
	
	#neg_jacobian_elist = [2201630,2201631,2187977,2187978,2017297,2017298,2017299,2022326,2022327,2022338,2022339]
	
	#neg_jacobian_elist = [2263294,2263295,2263296,2023565,2023566] # kp1b1
	#neg_jacobian_elist = [2262880,2262881,2262882,2022737,2022738] # kp1b4
	
	neg_jacobian_elist = [45987287,45987288,46100856,45979739,45979740,42116630,42116631,45989723,45989724,46100855,46370897,46370898] # kp3b1
	
	for ie_rea in neg_jacobian_elist:
		ie = ie_rea -1
		s6 = hex20.s6[ie]
		v8 = hex20.v8[ie]
		e12 = v8_to_e12(v8)
		hex20.e12[ie] = e12
		
		if 'PW ' in s6: continue
	
		# now, from ie do backward search, until s6 contains a 'PW ' flag
		while True:
			ie = ie - 1
			s6 = hex20.s6[ie]
			v8 = hex20.v8[ie]
			e12 = v8_to_e12(v8)
			hex20.e12[ie] = e12
			
			if 'PW ' in s6: break
		
def fix_negative_jacobian_elements_from_nek_for_cht_mesh(hex20,hex20_solid):
	# cht mesh needs some special treatment,
	# also need to linearize its contacting solid element

	#neg_jacobian_elist = [2053693,2185513]
	#neg_jacobian_elist = [2053695,2185513]
	neg_jacobian_elist = [2263294,2263295,2263296,2023565,2023566] # kp2b1
	
	for ie_rea in neg_jacobian_elist:
		ie = ie_rea -1
		s6 = hex20.s6[ie]
		v8 = hex20.v8[ie]
		e12 = v8_to_e12(v8)
		hex20.e12[ie] = e12
		
		if 'PW ' in s6:
			# linearize its touching solid elements
			# PW is always at face6 (5 in python)
			vc_f,fn,vquad,equad,ftag = hex20.return_face_info(ie,5)
			
			for ies in range(0,hex20_solid.nhex):
				vc_s,fn,vquad,equad,ftag = hex20_solid.return_face_info(ies,4)
				dist = distance(vc_f,vc_s)
				if(dist<0.001):
					v8 = hex20_solid.v8[ies]
					e12 = v8_to_e12(v8)
					hex20_solid.e12[ies] = e12
					print 'attaching solid element fixed'
					
			continue
	
		# now, from ie do backward search, until s6 contains a 'PW ' flag
		while True:
			ie = ie - 1
			s6 = hex20.s6[ie]
			v8 = hex20.v8[ie]
			e12 = v8_to_e12(v8)
			hex20.e12[ie] = e12
			
			if 'PW ' in s6: 
				# linearize its touching solid elements
				# PW is always at face6 (5 in python)
				vc_f,fn,vquad,equad,ftag = hex20.return_face_info(ie,5)
			
				for ies in range(0,hex20_solid.nhex):
					vc_s,fn,vquad,equad,ftag = hex20_solid.return_face_info(ies,4)
					dist = distance(vc_f,vc_s)
					if(dist<0.001):
						v8 = hex20_solid.v8[ies]
						e12 = v8_to_e12(v8)
						hex20_solid.e12[ies] = e12
						print 'attaching solid element fixed'
				break