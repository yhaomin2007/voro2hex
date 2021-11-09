def detect_jacobian_for_hex20(v8,e12):
	# convert hex20 to nek hex27
	s6,vc = hex20tohex27(v8,e12)

	# map xyz to rst and calcuate dirivative
	xr,xs,xt,yr,ys,yt,zr,zs,zt = get_dirivative_rst(v8,e12,s6,vc)

	# caclualte jacobian matrix determinent of each point in hex27
	jac27 = local_jacobian(xr,xs,xt,yr,ys,yt,zr,zs,zt)

	# detect vanishing jacobian and bad jacobian
	quality_threshold = 7e-3 # needed for nekrs mesh, from my experience
	if_vanishing,if_low_quality = detect_jacobian(jac27,quality_threshold)

	if (if_vanishing):
		print 'vanishing jacobian for this element, this elelemt will be linearized'

	if (if_low_quality):
		print 'low scaled-jacobian for this element, this elelemt will be linearized'

	if_linearize = False
	if if_vanishing or if_vanishing:
		if_linearize = True
		
	return if_linearize

def hex20tohex27(v8,e12):
	#1 convert to hex20
	# pack hex20
	hex27 = []
	for i in range(0,8):
		hex27.append(v8[i])

	for i in range(0,12):
		hex27.append(e12[i])
		
	s6 = []
	# patch hex27, but in nek5000 defintiion
	#2 convert to hex27

	# face 1
	quad = []
	quad.append(v8[0]) #1
	quad.append(v8[1]) #2
	quad.append(v8[5]) #6
	quad.append(v8[4]) #5

	quad.append(e12[0]) #1
	quad.append(e12[9]) #10
	quad.append(e12[4]) #5
	quad.append(e12[8]) #9
	
	fc = transfinite(quad)
	s6.append(fc)
	
	# face 2	
	quad = []
	quad.append(v8[1]) #2
	quad.append(v8[2]) #3
	quad.append(v8[6]) #7
	quad.append(v8[5]) #6

	quad.append(e12[1]) #2
	quad.append(e12[10]) #11
	quad.append(e12[5]) #6
	quad.append(e12[9]) #10
	
	fc = transfinite(quad)
	s6.append(fc)
	
	# face 3
	quad = []
	quad.append(v8[2]) #3
	quad.append(v8[3]) #4
	quad.append(v8[7]) #8
	quad.append(v8[6]) #7

	quad.append(e12[2]) #3
	quad.append(e12[11]) #12
	quad.append(e12[6]) #7
	quad.append(e12[10]) #11
	
	fc = transfinite(quad)
	s6.append(fc)
	
	# face 4
	quad = []
	quad.append(v8[3]) #4
	quad.append(v8[0]) #1
	quad.append(v8[4]) #5
	quad.append(v8[7]) #8

	quad.append(e12[3]) #4
	quad.append(e12[8]) #9
	quad.append(e12[7]) #8
	quad.append(e12[11]) #12
	
	fc = transfinite(quad)
	s6.append(fc)
	
	# face 5
	quad = []
	quad.append(v8[0]) #1
	quad.append(v8[1]) #2
	quad.append(v8[2]) #3
	quad.append(v8[3]) #4

	quad.append(e12[0]) #1
	quad.append(e12[1]) #2
	quad.append(e12[2]) #3
	quad.append(e12[3]) #4
	
	
	fc = transfinite(quad)
	s6.append(fc)

	# face 6
	quad = []
	quad.append(v8[4]) #5
	quad.append(v8[5]) #6
	quad.append(v8[6]) #7
	quad.append(v8[7]) #8

	quad.append(e12[4]) #5
	quad.append(e12[5]) #6
	quad.append(e12[6]) #7
	quad.append(e12[7]) #8
	
	fc = transfinite(quad)
	s6.append(fc)
	
	# volume center point
	quad = []
	quad.append(e12[8]) #9
	quad.append(e12[9]) #10
	quad.append(e12[10]) #11
	quad.append(e12[11]) #12

	quad.append(s6[0])  #1
	quad.append(s6[1])  #2
	quad.append(s6[2])  #3
	quad.append(s6[3])  #4
	
	fc1 = transfinite(quad)

	quad = []
	quad.append(e12[0]) #1
	quad.append(e12[4]) #5
	quad.append(e12[6]) #7
	quad.append(e12[2]) #3

	quad.append(s6[0])  #1
	quad.append(s6[5])  #6
	quad.append(s6[2])  #3
	quad.append(s6[4])  #5
	
	fc2 = transfinite(quad)
	
	
	quad = []
	quad.append(e12[1]) #2
	quad.append(e12[5]) #6
	quad.append(e12[7]) #8
	quad.append(e12[3]) #4

	quad.append(s6[1])  #2
	quad.append(s6[5])  #6
	quad.append(s6[3])  #4
	quad.append(s6[4])  #5
	
	fc3 = transfinite(quad)

	vc = []
	for i in range(0,3):
		vc.append((fc1[i]+fc2[i]+fc3[i])/3.0)

	return s6,vc

def get_dirivative_rst(v8,e12,s6,vc):
	# hex27 node order follow nek5000 definition.
	# caculate x,y,z dirivative to r,s,t.
	
	xr = []
	xs = []
	xt = []
	yr = []
	ys = []
	yt = []
	zr = []
	zs = []
	zt = []

	# shift it hex27 xyz to 3d array

	x = [[[0 for col in range(3)]for row in range(3)] for k in range(3)]
	y = [[[0 for col in range(3)]for row in range(3)] for k in range(3)]
	z = [[[0 for col in range(3)]for row in range(3)] for k in range(3)]

	x[0][0][0] =  v8[0][0]
	x[1][0][0] = e12[0][0]
	x[2][0][0] =  v8[1][0]
	x[0][1][0] = e12[3][0]
	x[1][1][0] =  s6[4][0]
	x[2][1][0] = e12[1][0]
	x[0][2][0] =  v8[3][0]
	x[1][2][0] = e12[2][0]
	x[2][2][0] =  v8[2][0]
	x[0][0][1] = e12[8][0]
	x[1][0][1] =  s6[0][0]
	x[2][0][1] = e12[9][0]
	x[0][1][1] =  s6[3][0]
	x[1][1][1] =     vc[0]
	x[2][1][1] =  s6[1][0]
	x[0][2][1] =e12[11][0]
	x[1][2][1] =  s6[2][0]
	x[2][2][1] =e12[10][0]
	x[0][0][2] =  v8[4][0]
	x[1][0][2] = e12[4][0]
	x[2][0][2] =  v8[5][0]
	x[0][1][2] = e12[7][0]
	x[1][1][2] =  s6[5][0]
	x[2][1][2] = e12[5][0]
	x[0][2][2] =  v8[7][0]
	x[1][2][2] = e12[6][0]
	x[2][2][2] =  v8[6][0]
	
	y[0][0][0] =  v8[0][1]
	y[1][0][0] = e12[0][1]
	y[2][0][0] =  v8[1][1]
	y[0][1][0] = e12[3][1]
	y[1][1][0] =  s6[4][1]
	y[2][1][0] = e12[1][1]
	y[0][2][0] =  v8[3][1]
	y[1][2][0] = e12[2][1]
	y[2][2][0] =  v8[2][1]
	y[0][0][1] = e12[8][1]
	y[1][0][1] =  s6[0][1]
	y[2][0][1] = e12[9][1]
	y[0][1][1] =  s6[3][1]
	y[1][1][1] =     vc[1]
	y[2][1][1] =  s6[1][1]
	y[0][2][1] =e12[11][1]
	y[1][2][1] =  s6[2][1]
	y[2][2][1] =e12[10][1]
	y[0][0][2] =  v8[4][1]
	y[1][0][2] = e12[4][1]
	y[2][0][2] =  v8[5][1]
	y[0][1][2] = e12[7][1]
	y[1][1][2] =  s6[5][1]
	y[2][1][2] = e12[5][1]
	y[0][2][2] =  v8[7][1]
	y[1][2][2] = e12[6][1]
	y[2][2][2] =  v8[6][1]
	
	z[0][0][0] =  v8[0][2]
	z[1][0][0] = e12[0][2]
	z[2][0][0] =  v8[1][2]
	z[0][1][0] = e12[3][2]
	z[1][1][0] =  s6[4][2]
	z[2][1][0] = e12[1][2]
	z[0][2][0] =  v8[3][2]
	z[1][2][0] = e12[2][2]
	z[2][2][0] =  v8[2][2]
	z[0][0][1] = e12[8][2]
	z[1][0][1] =  s6[0][2]
	z[2][0][1] = e12[9][2]
	z[0][1][1] =  s6[3][2]
	z[1][1][1] =     vc[2]
	z[2][1][1] =  s6[1][2]
	z[0][2][1] =e12[11][2]
	z[1][2][1] =  s6[2][2]
	z[2][2][1] =e12[10][2]
	z[0][0][2] =  v8[4][2]
	z[1][0][2] = e12[4][2]
	z[2][0][2] =  v8[5][2]
	z[0][1][2] = e12[7][2]
	z[1][1][2] =  s6[5][2]
	z[2][1][2] = e12[5][2]
	z[0][2][2] =  v8[7][2]
	z[1][2][2] = e12[6][2]
	z[2][2][2] =  v8[6][2]
	
	for irr in range (0,3):
		for iss in range (0,3):
			for itt in range(0,3):
				dxdr = dirivative(x[0][iss][itt],x[1][iss][itt],x[2][iss][itt],irr)
				dxds = dirivative(x[irr][0][itt],x[irr][1][itt],x[irr][2][itt],iss)
				dxdt = dirivative(x[irr][iss][0],x[irr][iss][1],x[irr][iss][2],itt)
				
				dydr = dirivative(y[0][iss][itt],y[1][iss][itt],y[2][iss][itt],irr)
				dyds = dirivative(y[irr][0][itt],y[irr][1][itt],y[irr][2][itt],iss)
				dydt = dirivative(y[irr][iss][0],y[irr][iss][1],y[irr][iss][2],itt)
				
				dzdr = dirivative(z[0][iss][itt],z[1][iss][itt],z[2][iss][itt],irr)
				dzds = dirivative(z[irr][0][itt],z[irr][1][itt],z[irr][2][itt],iss)
				dzdt = dirivative(z[irr][iss][0],z[irr][iss][1],z[irr][iss][2],itt)
				
				xr.append(dxdr)
				xs.append(dxds)
				xt.append(dxdt)
				yr.append(dydr)
				ys.append(dyds)
				yt.append(dydt)
				zr.append(dzdr)
				zs.append(dzds)
				zt.append(dzdt)

	return xr,xs,xt,yr,ys,yt,zr,zs,zt


def local_jacobian(xr,xs,xt,yr,ys,yt,zr,zs,zt):
	# calculate local jacobian matrix determinant for eeach grid point
	jac27 = []	
	for i in range (0,27):
		jac = 0.0
		jac = xr[i]*ys[i]*zt[i]
		jac = jac + xs[i]*yt[i]*zr[i]
		jac = jac + xt[i]*yr[i]*zs[i]
		jac = jac - xr[i]*yt[i]*zs[i]
		jac = jac - xs[i]*yr[i]*zt[i]
		jac = jac - xt[i]*ys[i]*zr[i]
		jac27.append(jac)
		
	return jac27

def detect_jacobian(jac27,quality_threshold):
	if_vanishing = False
	if_low_quality = False
	# 
	sign0 = jac27[0]
	for i in range (1,27):
		sign = jac27[i]
		if sign0*sign <= 0:
			if_vanishing = True
			return if_vanishing,if_low_quality
	
	minJ = min(jac27)
	maxJ = max(jac27) 
	if (minJ/maxJ) < quality_threshold:
		if_low_quality = True
		return if_vanishing,if_low_quality
		
	return if_vanishing,if_low_quality
	
def dirivative(x1,x2,x3,irr):
	#  irr is 0,1,2, 

	if irr==0:
		dxdr = (x2-x1)/1.0
	elif irr == 1:
		dxdr = (x3-x1)/2.0
	elif irr == 2:
		dxdr = (x3-x2)/1.0

	return dxdr


def transfinite(quad):
	# trasfinite interpolation face mid for from quad8
	fc = []
	for i in range(0,3):
		fcx = 0.5*(quad[4][i] + quad[5][i] + quad[6][i] + quad[7][i]) -0.25*(quad[0][i] + quad[1][i] + quad[2][i] + quad[3][i])
		fc.append(fcx)
	
	return fc