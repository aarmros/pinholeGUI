from detectorClass import *
import csv
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.gridspec as gridspec
import signalGen3
from math import *
import eqtools 
import time

import pidly #for testing with IDL

"""Sep 22, 2017  - Version 3 adds eqtools for magnetic reconstruction; mostly just added to function
calls to work as pass through to signalGen2 where actual changes were made"""

# step size for generating line integrals
numStep = 50

balmer2Lyman =1 #correcting from balmer to lyman

""" used for setting the alpha level when plotting signal to noise. 
    Very Simple."""
def aboveOne(num):
	if num > 1.0:
		return 1.0
	else:
		return num

"""angle measured from 3 o'clock towards midnight. Inteded for giving
    coordinates of detector image plane in tokomak coordinates.
    Used for the poloidal view calculation. Note -deltaR due to geometry."""
def rZCoord(r,z,distance,angle):
	deltaR = sin(angle) * distance
	deltaZ = cos(angle) * distance
	return [r-deltaR,z+deltaZ]


"""finds the edges of a box which incloses the emissitivity data. 
   returns in the form [rMin,rMax, zMin,zMax]. emiss is of the form
   [[r0,z0,emiss0],[r1,z1,emiss1],...]"""
def boundaries(emiss):
	rMin = emiss[0][0]
	rMax = emiss[0][0]
	zMin = emiss[0][1]
	zMax = emiss[0][1]

	#Go over all the data updating as you go
	for i in range(len(emiss)):
		curR = emiss[i][0]
		curZ = emiss[i][1]

		if (curR < rMin):
			rMin = curR

		elif (curR > rMax):
			rMax = curR

		elif (curZ < zMin):
			zMin = curZ

		elif ( curZ > zMax):
			zMax = curZ

	return [rMin,rMax, zMin, zMax]

"""returns the r point which hits boundary of the permiter
   defined by boundaries above on the line defined by r,z,slope"""
def findBound(r,z,slope,zMin,zMax, curR):

	# if the line hits the top boundary, we return the r point
	# where it hits the boundary
	if (zMax < signalGen3.line([slope,r,z],curR)):
		return (zMax-z)/slope + r

	# same thing for the bottom boundary
	elif( zMin > signalGen3.line([slope,r,z],curR)):
		return (zMin - z)/slope + r

	#otherwise it hits the boundaries outside our range of interest
	else:
		return curR


""" Why does plotting have to take so many lines :( Anyways, this plots
    for toroidal view"""
def torPlot(r, z, angle, cov, obDist, detArray,lineEmiss,brightness, emiss):

	rMin = r - (sin(angle) * cov/2)
	zMin = z - (cos(angle) * cov/2)

	# creat the figure
	fig1 = plt.figure()
	gridspec.GridSpec(3,4)

	# display SOLPS data 
	p1 = plt.subplot2grid((3,4),(0,0), colspan = 2, rowspan = 3)
	plt.text(0.8,1.4,'Detector Coverage in Red')
	plt.xlabel('R [m]')
	plt.ylabel('Z [m]')
	plt.title('SOLPS data')
	plt.scatter( [row[0] for row in emiss],[row[1] for row in emiss], \
		c = [row[2] for row in emiss])

	# adding in array coverage graphic, overlay in SOLPS data 
	# Note this is all to scale...pretty dope
	detArrayHeight = detArray[0].correctSize(obDist,detArray[0].sHeight, detArray[0].height)

	# A Lengthy note about geometry: the rectangle command takes coordinates which define the 
	# lower left coordinate. So we calculate the lower left coordinate of the detector defined 
	# by r,z,angle and detector Height. We then rotate it to the correct angle. Rectangles 
	# angle is anti clockwise and must be given in degrees
	p1.add_patch(patches.Rectangle((rMin+(detArrayHeight/2*cos(angle)),zMin-(detArrayHeight/2 * sin(angle))),\
	cov, detArrayHeight, angle = 90 - degrees(angle), color ='red',alpha = 0.75))

	# first creating the signal to noise plot
	p2 = plt.subplot2grid((3,4),(2,2), colspan = 2, rowspan = 1)
	plt.plot([row[0] for row in brightness],[1.0 for row in brightness],\
		label = 'sig/noise',color = 'red')
	plt.xlabel('R [m]')
	plt.ylabel('Sig/Noise')

	# Each bar is the width of a detector ob plane, positioned at the r position 
	# in the ob plane. The height is sig/noise.
	for i in range(len(detArray)):
		curR = r + (detArray[i].obPlane(obDist)*sin(angle))
		w = sin(angle)*detArray[i].correctSize(obDist, detArray[i].sWidth, detArray[i].width)

		sig = detArray[i].signal
		vNoise  = detArray[i].noise
		#add in a bar for each detector
		p2.bar(curR-w, sig/vNoise,w,alpha = aboveOne(sig/vNoise))	
		p2.text(curR,0,str(i))


	#Now create the subplot for the brightness
	p3 = plt.subplot2grid((3,4),(1,2), colspan = 2, rowspan = 1)
	plt.plot([row[0] for row in brightness],[row[1] for row in brightness],label = 'brightness')
	plt.ylabel('Brightness \n [photons/(m^2 sr s)]')
	plt.legend()

	#Finally the emissitivity
	p4 = plt.subplot2grid((3,4),(0,2), colspan = 2, rowspan = 1)
	plt.ylabel('Emiss \n [photons/(m^3 sr s)]')
	plt.plot([row[0] for row in lineEmiss],[row[2] for row in lineEmiss],label = 'emiss')
	plt.legend()

	#Display it!
	fig1.tight_layout()
	plt.show()

	return 0 

# Once again plotting is the worst. This time it is for the poloidal plot.
def polPlot(r,z,angle, obDist, detArray, emiss, lineEmissArray):

	sampleDet = detArray[0] #assuming array is all the same

	horizontal = sampleDet.obPlaneLength(obDist)

	# For drawing the object plane, NOTE this is slightly incorrect since
	# we don't factor in the detector height as we do in tor plot
	rMin = rZCoord(r,z,-horizontal/2,angle)[0]
	zMin = rZCoord(r,z,-horizontal/2,angle)[1]


	# creat the figure
	fig1 = plt.figure()
	gridspec.GridSpec(3,4)

	# display SOLPS data 
	p1 = plt.subplot2grid((3,4),(0,0), colspan = 2, rowspan = 3)
	plt.text(0.8,1.3,'Detector Coverage in Red \n Lines of Sight in Green')
	plt.xlabel('R [m]')
	plt.ylabel('Z [m]')
	plt.title('SOLPS data')
	plt.scatter( [row[0] for row in emiss],[row[1] for row in emiss], \
		c = [row[2] for row in emiss])

	# Overlay detector coverage graphic to SOLPS data. To scale!!
	p1.add_patch(patches.Rectangle((rMin,zMin),\
	horizontal, 0.01, angle = 90 + degrees(angle), \
	color ='red',alpha = 0.75,zorder = 3))
	
	# draw the lines of sight of each detector
	for i in range(len(lineEmissArray)):
		p1.plot([row[0] for row in lineEmissArray[i]], \
			[row[1] for row in lineEmissArray[i]],color = 'green')

	# Now create the signal to noise plot
	p2 = plt.subplot2grid((3,4),(2,2), colspan = 2, rowspan = 1)
	plt.xlabel('Distance from Object Plane Center [m]')
	plt.ylabel('Sig/Noise')
	oneLine = []

	for i in range(len(detArray)):
		#calculate the distance and width of each detectors
		dist = detArray[i].obPlane(obDist)  
		w = detArray[i].correctSize(obDist,detArray[i].sWidth, detArray[i].width)
		sig = detArray[i].signal
		vNoise  = detArray[i].noise
		#add in a bar for each detector
		p2.bar(dist - w, sig/vNoise,w,alpha = aboveOne(sig/vNoise))	
		p2.text(dist,0,str(i))

		#this is for drawing a line of 1
		oneLine.append(dist)

	plt.plot([oneLine[i] for i in range(len(oneLine))],[1.0 for i in range(len(oneLine))],\
		label = 'sig/noise',color = 'red')
	
	# Now for the emissivity, we will display what each detector sees
	p3 = plt.subplot2grid((3,4),(0,2), colspan = 2, rowspan = 2)

	for i in range(len(lineEmissArray)):
		p3.plot([row[0] for row in lineEmissArray[i]], [row[2] for row in lineEmissArray[i]],\
			label = str(i))


	plt.ylabel('Emiss \n [photons/(m^3 sr s)]')
	plt.xlabel('R [m]')
	plt.title('Emissitivity for Detector Lines of Sight')
	p3.legend(loc = 'upper center', bbox_to_anchor = (0.95,0.95), ncol = 2)

	#Display it!
	fig1.tight_layout()
	plt.show()


	return 0 

"""Calculates the brightness for each detector in detArray. The detArray is located
   a distance obDist away from the center of the object plane defined by r,z, and angle.
   This is for a toroidal view and uses an inverse abel transform. """

def torBrightness(r,z, angle, emiss, detArray, obDist,eqObj):
	# Just avoiding any divide by zeroes limited to pi rotation due to symmetry
	if(angle <= 0.0) or (angle >= pi):
		angle = pi/2
		print('ERROR: Allowed toroidal angle range from 0 to Pi')

	sampleDet = detArray[0] #assuming array is all the same

	# the distance covered by the array is the distance to the farthest
	# detector object plane (center of detector object plane) + detector object width/2
	length = sampleDet.obPlaneLength(obDist)

	# for finding where to stop our grid
	boundary = boundaries(emiss)

	rMin = r - (sin(angle) * length/2) #farthest coordinates of detector object plane
	zMin = z - (cos(angle) * length/2)
	rMax = boundary[1]

	print('rMax: '+ str(rMax))
	print('rMin: ' + str(rMin))
	print('step: ' + str(length/numStep))

	# Here we are calculating the object plane slope
	slope = (z - zMin)/(r-rMin)

	# if the detector is horizontal, it is much faster to treat it seperately because
	# the torAbel calculation is one dimensional
	if(slope == 0.0):
		#Find all the data in the object plane slope,
		lineEmiss = signalGen3.lineData(r,z,slope,rMin,rMax, length/numStep, emiss,eqObj)

		#to use abel transform,we assume toroidal symmetry
		brightness = signalGen3.torAbel(lineEmiss)

	# for each detector, we need to find the coordinates of the plane it is imaging
	# then use the abel transformed brightness to calculate the signal.

	# We are technically assuming that each detector has a parallel line of sight
	for i in range(len(detArray)):

		#finding the coordinates
		curR = r + (detArray[i].obPlane(obDist)*sin(angle))
		curZ = z + (detArray[i].obPlane(obDist) * cos(angle)) # if horizontal all the same z
		w = sin(angle)*detArray[i].correctSize(obDist, detArray[i].sWidth, detArray[i].width)/2

		print('range: '+str(curR - w)+' to '+str(curR + w))

		# if the detector is horizontal, the Able Transform saves as computation time
		if(slope == 0.0):
			curBright = signalGen3.boxAve(brightness, curR - w, curR + w)

		# if it is at an angle, we have to use the axially symmetric 3D abel transfrom
		else:
			lineEmiss = signalGen3.lineData(curR,curZ,0.0,curR,rMax, length/numStep,emiss,eqObj)
			brightness = signalGen3.torAbel(lineEmiss)
			curBright = brightness[0][1]
			# Precisely, we would box average for each detector over a range of Z at each grid 
			# point; however this is really computationally expensive and uncessary 
			# with the sparsness of the SOLPS data


		print('brigtness' + str(curBright) ) 

		detArray[i].sigCalc(balmer2Lyman*curBright)# this 20 is from assuming solps data is for balmer
		
		
		sig = detArray[i].signal
		vNoise  = detArray[i].noise

		print('signal:'+str(sig/vNoise)+'\n')

	# Now plot it
	torPlot(r, z, angle, length, obDist, detArray,lineEmiss,brightness,emiss)


	return detArray

"""Calculates the brightness for each detector in detArray. The detArray is located
   a distance obDist away from the center of the object plane defined by r,z, and angle.
   This is for a poloidal view and integrates every line of sight"""

def polBrightness(r,z, angle, emiss, detArray, obDist,eqObj):
	sampleDet = detArray[0] #assuming array is all the same

	horizontal = sampleDet.obPlaneLength(obDist)

	#Find the coordinates of r at the edges of the object plane
	rMin = rZCoord(r,z,-horizontal/2,angle)[0]
	rMax = rZCoord(r,z,horizontal/2,angle)[0]


	print('rMax: '+ str(rMax))
	print('rMin: ' + str(rMin))

	# This is the poisition of the pinhole in tokomak space
	# We use it to define the lines of sight of each detector
	pinholeR = obDist * cos(angle) + r
	pinholeZ = obDist * sin(angle) + z

	# for plotting to give to polPlot
	lineEmissArray = []

	# for finding where to stop our grid
	boundary = boundaries(emiss)

	# make sure the grid generated doesn't end too early
	padding = horizontal/2


	# We will now do the line integration for each detector in the array
	for i in range(len(detArray)):
		curCoords = rZCoord(r,z,detArray[i].obPlane(obDist),angle)
		curR = curCoords[0]
		curZ = curCoords[1]

		# each detector's line of sight is defined by the object plane
		# and the pinhole location
		slope = (curZ - pinholeZ)/(curR-pinholeR) 

		# We only want to integrate over a line a little larger(by padding)
		# than the bounds of the emissivity data
		minR = findBound(curR,curZ, slope, boundary[2]+padding,\
			boundary[3]+padding,boundary[0])
		maxR = findBound(curR,curZ, slope, boundary[2]+padding,\
			boundary[3]+padding,boundary[1])
	
		step = (maxR-minR)/numStep #gurantees numStep number of steps
		
		# load up the emissivity data onto our line of sight
		lineEmiss = signalGen3.lineData(curR,curZ, slope, minR,maxR, step, emiss,eqObj)

		#load it for plotting
		lineEmissArray.append(lineEmiss)

		# find the signal, by integrating emissivity
		curBright = signalGen3.lineInt(lineEmiss)

		detArray[i].sigCalc(balmer2Lyman*curBright)


	polPlot(r,z,angle, obDist, detArray, emiss, lineEmissArray)

	return detArray




"""detPara = [ detWidth,detHeight, nDetectors, detGap, responseAXUV,filterTrans,lyALine,filterFWHM,gain,
filterN,voltNoise]  all in meters"""
def runSim(sigZ,sigR, angle, obDist, resolution, sWidth, sHeight, tor, detParam,emiss,eq):

	t0 = time.time()
	detW = detParam[0] 
	imDist = obDist * detW/ resolution
	detArray = []
	nDet = detParam[2]
	# We will assume the pinhole marks x-center for the detector coordinates
	# The detector coordinates are not in tokomak space, see detector class
	# for definition of coordinates for Detectors.

	for i in range(nDet):

		#This formula is a bit confusing but gives correct result for even or odd detector
		#arrays, with center of array at (0,imDist)
		xDist = (-(nDet-1)/2.0 + i) * (detW + detParam[3]) 
		
		# initialize them with the values
		det = Detector(xDist,-imDist,detW, detParam[1], detParam[4], 0.0, sWidth, sHeight,\
	              detParam[5],detParam[6],detParam[7], detParam[8],detParam[9],detParam[10])

		detArray.append(det)

	"""	
	#load our data
	lEmiss = signalGen3.myRead(fileNames[0])		
	rPos = signalGen3.myRead(fileNames[1])
	zPos = signalGen3.myRead(fileNames[2])

	#magnetic reconstruction for flux interpolation
	eq = eqtools.EqdskReader(gfile = 'g' + str(fileNames[3]),afile = 'a'+ str(fileNames[3]))

	#Organize it for how the functions are written
	emiss = []
	for i in range(len(lEmiss)):
		curPsi = eq.rz2psinorm(np.array([rPos[i]]),np.array([zPos[i]]))
		curAng = signalGen3.magAngle(rPos[i],zPos[i],eq)
		emiss.append([ rPos[i],zPos[i],lEmiss[i],curPsi,curAng])

	"""

	#We have two ways of calculating the Brightness for each detector
	# With a toroidal view (tor == True) or poloidal view
	# the toroidal view uses an Abel Inversion. 
	if (tor == True):
		torBrightness(sigR,sigZ, angle, emiss, detArray, obDist,eq)
		t1 = time.time()

		print('time to run: '+ str(t1-t0))
		return 0 

	# the poloidal view does a line integral
	else:
		polBrightness(sigR,sigZ,angle, emiss, detArray, obDist,eq)
		t1 = time.time()
		print('time to run: '+ str(t1-t0))
		return 0

	
##########         FOR TESTING               ##############
##########         FOR TESTING               ##############
##########         FOR TESTING               ##############
##########         FOR TESTING               ##############




def main():

	"""

	# All measurements given in m
	lyALine = 0.000121567 *10**(-3) # wavelength in m
	eLyman = 1.63*10**(-18) # energy of lyman photon J


	# We will assume that the detector quantities are fixed 
	detWidth = 0.75 *10**(-3) # in m 
	detHeight = 4 *10**(-3) ## in m
	detGap = 0.20 *10**(-3)#Gap between detectors in array
	nDetectors = 20 # number of detectors
	detLength = nDetectors*detWidth+(nDetectors-1)*detGap
	responseAXUV = 0.15 #AXUV respon

	filterTrans = 0.05  #We could get this in 10 though with great FWHM
	filterFWHM = 0.00001 *10**(-3)# m, must match units of center
	filterN = 2.6 #filter effective refractive index
	balmerBright = 5*10**18 #taken directly from matlab (ph/s sr m2)
	lymanBright = 20 * balmerBright

	voltNoise = 20*10**(-3)
	gain = 1*10**7 #set by pre-amp V/A

	sigZ = -0.25
	sigR = 2.25
	
	obDist = 1000 *10**(-3)
	resolution = 5 *10**(-3)
	sWidth = 2 *10**(-3)
	sHeight = 8 *10**(-3)
	imDist = obDist * detWidth / resolution
	detArray = []
	# We will assume the pinhole marks x-center
	"""
	#for loading in stuff from idl

	idl = pidly.IDL()

	idl('restore, "~/shot1091022012_1905.sav"')
	lineEmiss = []
	x = []

	for i in range(len(idl.emiss0)):
		lineEmiss.append([idl.emiss0[i],-0.038,idl.e[i]])
		x.append(idl.emiss0[i])

	lineEmiss.append([0.931,-0.038,0])

	print(lineEmiss)
	print('\n')
	
	


	abel = signalGen3.torAbel(lineEmiss)
	"""
	abel = []
	for i in range(len(lineEmiss)-1 ,-1 ,-1):
		integrand = 0.0
		y = lineEmiss[i][0]

		print '\n'

		for j in range(i , len(lineEmiss)):

			r = lineEmiss[j][0] 

			if r != y:  #make sure it is finite
				dr = lineEmiss[len(lineEmiss)-1][0] - lineEmiss [len(lineEmiss)-2][0]
				print('r is:' + str(r))
				print('y is:' + str(y))
				if (j != len(lineEmiss)-1):
					print('not default dr')
					dr = lineEmiss[j][0] - lineEmiss [j-1][0]
					
				integrand += (2 * lineEmiss[j-1][2] * r * dr /sqrt(r**2 - y**2))
				
		
		# due to counting backwards insert first element, this is bad performance wise
		abel.insert(0,[y,integrand])
	


	for i in range(len(lineEmiss)-1 ,-1 ,-1):
		integrand = 0.0
		y = lineEmiss[i][0]

		for j in range(i+2 , len(lineEmiss)):

			r = lineEmiss[j][0] 

			if r != y:  #make sure it is finite
				dr = lineEmiss[j][0] - lineEmiss [j-1][0]
				r = (lineEmiss[j][0] +lineEmiss[j-1][0])/2.0 
				print('r is:' + str(r))
				print('y is:' + str(y))

				emissR = (lineEmiss[j][2]+lineEmiss[j-1][2])/2.0
					
				#integrand += (2 * emissR * r * dr /sqrt(r**2 - y**2))

				integrand += ((2 * lineEmiss[j][2] * r/sqrt(r**2 - y**2)) + (2 * lineEmiss[j-1][2]/sqrt(lineEmiss[j-1][0]**2 - y**2)))*dr/2
				
		
		# due to counting backwards insert first element, this is bad performance wise
		abel.insert(0,[y,integrand])
	"""
	print abel
	
	plt.plot([row[0] for row in abel],[row[1] for row in abel])
	plt.show()

	idl.close()


	"""

	for i in range(nDetectors):
		#This formula is a bit confusing but gives correct result for even or odd detector
		#arrays, with center of array at (0,imDist)
		xDist = (-(nDetectors-1)/2.0 + i) * (detWidth + detGap) 
		
		det = Detector(xDist,-imDist,detWidth, detHeight, responseAXUV, 0.0, sWidth, sHeight,\
	              filterTrans,lyALine,filterFWHM, gain,filterN, voltNoise)

		detArray.append(det)

		

	lEmiss = signalGen2.myRead("SOLPS/lyman_alpha_data.txt")		
	rPos = signalGen2.myRead("SOLPS/lyman_alpha_R.txt")
	zPos = signalGen2.myRead("SOLPS/lyman_alpha_Z.txt")

	emiss = []
	for i in range(len(lEmiss)):
		emiss.append([ rPos[i],zPos[i],lEmiss[i] ])

	
	torBrightness(sigR,sigZ,pi/2, emiss, detArray, obDist)
	"""


	return 0

if __name__ == "__main__":
	main()