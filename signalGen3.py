import csv
import numpy as np
import matplotlib.pyplot as plt
from math import * 
import eqtools
import time #for testing runtime



"""
Sep 22, 2017 - Version 2 adds flux interpolation along with the helper functions magAngle and 
aveEmiss

Sep 30 2017 - Version 3 magAngle now takes arrays, passing arrays is faster because it minimizes
              calls to eqtools
              FluxInterpolate now assumes, theta and psi are in emiss array to minimize calls to 
              eqtools
             - moved myRead to pinholeGUI since file reading of program happens there
               now pinholeGUI is imported for testing in this file
"""

maxSearch = 100


""" We live in a 2D world """
def distance(x,y,x2,y2):
	return sqrt((x-x2)**2 +(y-y2)**2)


"""slopePoint = [ slope, x1, y1] finds the value at x of line defined
   by slopePoint"""
def line(slopePoint,x):
	return slopePoint[0]*(x-slopePoint[1])+slopePoint[2]

""" takes in r and z position as array and eqtools eqDSKReader object. Calculates the angle 
    from the horizontal with the magnetic center as the axis. returned angle is
    from -pi to pi """

def magAngle(r,z,eqObj):

	if (len(r) != len(z)):
		print('magAngle Error: R and Z arrays have different length')

	# find the magnetic center
	magR = eqObj.getMagR()
	magZ = eqObj.getMagZ()

	angles = []

	# go over inputs finding the angles
	for i in range(len(r)):

		x = r[i] - magR
		y = z[i] - magZ

		# We use atan2 because it knows about quadrants and returns a signed 
		# result which makes minimizing theta in fluxInterpolate easier
		angles.append(atan2(y,x))

	
	return angles

"""" returns weighted average by distance"""
def aveEmiss(emiss1,emiss2,dis1,dis2):

	dTotal = dis1+dis2 # total distance

	return (emiss1 * dis2/dTotal) + (emiss2 * dis1/dTotal)

""" Takes in a point and its corresponding psi, goes over emiss finding closest point
    in positive and negative theta direction on the same flux surface. Returns
    the distance weighted average emissivity of the two points. """

def fluxInterpolate(r,z, psi, emiss, maxDist,eqObj):

	closestR1 = r              # set theses in the loop to correct value
	closestZ1 = z 
	closestEmiss1 = 0
	disNegTheta = 1000         # we are minimizing; so good place to start

	# update this somehow to use distance in emiss values Sep 22 2017

	# It takes two to interpolate
	closestR2 = r
	closestZ2 = z
	closestEmiss2 = 0
	disPosTheta= 1000

	theta = magAngle([r],[z],eqObj)[0]

	tolerance = 0.01           # this should be tweeked depending on spacing of data sep 22 2017

	# looking for closest data point with smaller theta and larger theta on same psi in emiss
	for i in range(len(emiss)):
		
		curPsi = emiss[i][3]

		#We can just skip the point if it is not on same psi

		if (curPsi  > (psi + tolerance)) or (curPsi  < (psi - tolerance)):
			continue

		curR = emiss[i][0]
		curZ = emiss[i][1]
		curTheta = emiss[i][4]
		curDis = distance(r,z,curR,curZ)

		# first check if curTheta is in positive or negative theta direction
		if  (curTheta - theta) >=  0 and (curTheta - theta) < pi :

			# Now check if it is closest point 
			if(curDis < disPosTheta):
			
				closestR2 = curR
				closestZ2 = curZ
				closestEmiss2 = emiss[i][2]
				disPosTheta = curDis

		#otherwise we are looking in negative theta direction
		else:


			if(curDis < disNegTheta):

				closestR1 = emiss[i][0]
				closestZ1 = emiss [i][1]
				closestEmiss1 = emiss[i][2]
				disNegTheta = curDis


	# Do a weighted average of closest points
	aEmiss = aveEmiss(closestEmiss1,closestEmiss2, disNegTheta,disPosTheta)

	# If one of the closest data points is too far away, ignore it
	if (disNegTheta > maxDist):
		
		print('One No point within maxDist')

		closestEmiss1 = 0.0
		aEmiss = closestEmiss2 #ignore the far point in average


	if (disPosTheta > maxDist):

		print(' Two No point within maxDist')

		closestEmiss2 = 0.0
		aEmiss = closestEmiss1

	return [closestR1,closestZ1,aEmiss]



"""closestPoint is no longer used as of Version3, keeping becasue
   it is generally useful and maybe used in future """

"""searches emiss until it finds the closest point in emiss
   to x, y, returns that entry of emiss i.e. [r/x ,z/y,emiss Value]. If
   the closestpoint is farther than MaxDist it returns the input
   x,y and emissitivity of 0"""

def closestPoint(x,y, emiss, maxDist):
	t0 = time.time()

	closestX = x # set theses in the loop to correct value
	closestY = y # we set equal to x, y in case of not finding a point
	closestEmiss = 0 

	dist = distance(x,y,emiss[0][0],emiss[0][1])

	for i in range(len(emiss)):
		curDist = distance(x,y,emiss[i][0],emiss[i][1])
		#trying to find the minimum distance from x,y and return it
		if (curDist < dist and curDist < maxDist):
			
			closestX = emiss[i][0]
			closestY = emiss [i][1]
			closestEmiss = emiss[i][2]
			dist = curDist # our curDist is shortest

	# outside this range, we will return defaul values of x,y, 0
	if (dist > maxDist):
		print('No point within maxDist')

	t1 = time.time()

	print('time to run: ' + str(t1-t0))

	return [closestX,closestY,closestEmiss]


"""returns an array of data points on the line defined by x,y and slope
   in the spacing of delta"""
def genGrid(x, y, slope, delta, xMin, xMax):

	numSteps = int(((xMax - xMin)/delta)+2) # needs to be int to pass to range
	# the +2 is there so that we create a grid that fills all of xMin to xMax
	xData = [((i * delta) + xMin) for i in range(numSteps)]

	grid = []

	for i in range(len(xData)):
		yval = line([slope,x,y],xData[i])

		grid.append([xData[i],yval])


	return grid

# creates a line of x,y slope between xmin and xmax and finds
# the emissitivity for every point on that line. Returns an array
# ordered by x of the form [[x,y, emiss],[x,y,emiss]...] where x and y
# are on the line
def lineData(x,y,slope,xMin,xMax, delta, emiss, eqObj):

	grid = genGrid(x,y,slope,delta,xMin,xMax)

	eqCalls = 0

	lineData = []

	for i in range(len(grid)):
		curX = grid[i][0]
		curY = grid[i][1]
		curPsi = eqObj.rz2psinorm(np.array([curX]),np.array([curY]))
		#curEmiss = fluxInterpolate(curX,curY,emiss,eqObj,maxSearch*delta)[2]
		fI = fluxInterpolate(curX,curY,curPsi,emiss,maxSearch*delta,eqObj)
		curEmiss = fI[2]
		lineData.append([curX,curY,curEmiss])


	return lineData


"""toroidel abel transform calculator. Takes in lineEmiss [[x,y,emiss],...]
   and returns the abel transorm [[y,emissTransfrom],...]. This assumes lineEmiss
    is ordered by lineEmiss[i][0]"""
def torAbel(lineEmiss):

	abel = []

	# we assume there is an equal spacing
	dr = lineEmiss[1][0] - lineEmiss [0][0]


	for i in range(len(lineEmiss)-1 ,-1 ,-1):
		integrand = 0.0
		y = lineEmiss[i][0]

		# the i+1 is to avoid the case where r = y and the denominator blows up 
		# in the integral
		
		for j in range(i + 1 , len(lineEmiss)):

			r = lineEmiss[j][0] 

			dr = lineEmiss[j][0] - lineEmiss [j-1][0]
			r = (lineEmiss[j][0] +lineEmiss[j-1][0])/2.0 # use the midpoint

			emissR = (lineEmiss[j][2]+lineEmiss[j-1][2])/2.0 #average

			integrand += (2 * emissR * r * dr /sqrt(r**2 - y**2))
		
		# due to counting backwards insert first element, this is bad performance wise
		abel.insert(0,[y,integrand])

	return abel


""" this is for poloidal detectors were we are just calculating the line integral
    the line integral is already parameterized by linEmiss. This assumes lineEmiss
    is ordered by path of integration"""
def lineInt(lineEmiss):

	integrand  = 0.0 

	#We don't go over the full range because of indexing for calculating dr
	for i in range(len(lineEmiss)-1):

		dr = sqrt((lineEmiss[i+1][0] - lineEmiss [i][0])**2 + 
		(lineEmiss[i+1][1] - lineEmiss [i][1])**2)

		integrand +=  (lineEmiss[i][2]*dr)

	#take care of the last one 
	i = len(lineEmiss) - 1

	dr = sqrt((lineEmiss[i][0] - lineEmiss [i-1][0])**2 + 
		(lineEmiss[i][1] - lineEmiss [i-1][1])**2)

	integrand +=  (lineEmiss[i][2]*dr)

	return integrand

"""assumes that f is of the form [[x,fx],[x,fx]...] and xMin and xMax are within
   the range of x given in f"""
def boxAve(f,xMin,xMax):
	# we assume there is an equal spacing, and lineEmiss is ordered
	dx = f[1][0] - f[0][0]
	integrand  = 0.0
	totalX = 0.0

	if (xMin < f[0][0]) or (f[-1][0] < xMax):
			print('ERROR: xMin or xMax out of range of f')
			return 0.0

	# go over all of f, not super fast if f is large array
	for i in range(len(f)):
		curX = f[i][0]

		if (xMin <= curX) and (curX <= xMax):
			integrand +=   (f[i][1] * dx)
			totalX += dx
	return integrand / totalX




##########         FOR TESTING               ##############
##########         FOR TESTING               ##############
##########         FOR TESTING               ##############
##########         FOR TESTING               ##############

def main():

	"""

	a = genGrid(1,0,0.0,0.01,0.75,1.5)

	test = []

	for i in range(len(a)):
		if (a[i][0] < 1):
			test.append([a[i][0],a[i][1],1])
		else:
			test.append([a[i][0],a[i][1],0])

	print(test)
	b = torAbel(test)
	print(b)

	plt.plot([row[0] for row in b],[row[1] for row in b])
	plt.show()




	a = genGrid(0,0,2,0.01,-0.75,0.75)
	print(a)
	print('\n')

	test = []

	for i in range(len(a)):
		dr  = sqrt((a[i][0])**2 + (a[i][1])**2)

		if (dr <= 0.5):
			test.append([a[i][0],a[i][1],1])
		else:
			test.append([a[i][0],a[i][1],0])

	print(test)
	print('\n')

	print(lineInt(test))
	

	"""
	lEmiss = pinholeGUI4.myRead("SOLPS/lyman_alpha_data.txt")		
	rPos = pinholeGUI4.myRead("SOLPS/lyman_alpha_R.txt")
	zPos = pinholeGUI4.myRead("SOLPS/lyman_alpha_Z.txt")

	edr = eqtools.EqdskReader(gfile = 'g156867.02006',afile = 'a156867.02006')
	edr.plotFlux()

	"""

	emiss = []
	for i in range(len(lEmiss)):
		emiss.append([ rPos[i],zPos[i],lEmiss[i] ])
	
	r = 1.48
	z = 0.91
	fluxInterpolate(r,z,emiss,edr,0.5)
	closestPoint(r,z,emiss,0.5)
	
	
	dta = []
	maxr = 0
	maxz = 0
	maxe = 0 

	for i in range(len(emiss)):
		if emiss[i][0] < 2.4 and emiss[i][0] > 2.00 :
			if emiss[i][1] < -0.24 and emiss[i][1] > -0.35 :
				dta.append([emiss[i][0],emiss[i][1],emiss[i][2]])

				if emiss[i][2] > maxe:
					maxe = emiss[i][2]
					maxr = emiss[i][0]
					maxz = emiss[i][1]

	print('max z:' + str(maxz))
	print('max r:' + str(maxr))

	print(edr.rz2psinorm(np.array([maxr]),np.array([maxz])))
	print(edr.rz2psinorm(np.array([2.24]),np.array([-0.189])))
	p1 = plt.subplot(211)
	
	plt.scatter([row[0] for row in dta],[row[2] for row in dta])
	plt.subplot(212)
	plt.scatter([row[0] for row in dta],[row[1] for row in dta], c= [row[2] for row in dta])
	plt.show()

	"""

	"""

	data = torAbel(lineData(1.1,-0.015,0.0,0.001,1.1,1.15,emiss))

	plt.plot([row[0] for row in data],[row[1] for row in data])
	val = boxAve(data,1.11,1.15)
	yval = [val for i in range(len(data))]
	plt.plot([row[0] for row in data],yval)
	plt.show()
	"""
	


if __name__ == "__main__":
	main()
