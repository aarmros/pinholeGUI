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

#import pidly #for testing with IDL

"""Sep 22, 2017  - Version 3 adds eqtools for magnetic reconstruction; mostly just added to function
calls to work as pass through to signalGen2 where actual changes were made"""

# step size for generating line integrals
numStep = 50

balmer2Lyman =1 #correcting from balmer to lyman

noiseFrac = 1.0 #fraction of carbon III of Lyman
noiseWvlngth = 117.5*10**(-9)#carbon III VUV line is of concern in nm

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
	rMax = r + (sin(angle) * cov/2)
	zMax = z + (cos(angle) * cov/2)

	# creat the figure
	fig1 = plt.figure()
	gridspec.GridSpec(3,4)

	plotEmiss = emiss
	
	# display SOLPS data, we only want data close to area of interest
	
	plotEmiss = []
	for i in range(len(emiss)):
		if emiss[i][0] < (rMax + 0.2):
			if emiss[i][0] > (rMin - 0.2):
				if emiss[i][1] > (zMin - 0.2):
					if emiss[i][1] < (zMax + 0.2):
						plotEmiss.append([emiss[i][0],emiss[i][1],emiss[i][2]])
	


	p1 = plt.subplot2grid((3,4),(0,0), colspan = 2, rowspan = 3)
	plt.text(0.8,1.4,'Detector Coverage in Red')
	plt.xlabel('R [m]')
	plt.ylabel('Z [m]')
	plt.title('SOLPS data')
	plt.scatter( [row[0] for row in plotEmiss],[row[1] for row in plotEmiss], \
		c = [row[2] for row in plotEmiss])
	plt.colorbar()

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

	

	#below is for sig/noise bar chart
	"""
	plt.plot([row[0] for row in brightness if row[0] < rMax],[1.0 for row in brightness if row[0] < rMax],\
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
	"""

	rPlot = []
	lSig = []
	cSig = []

	for i in range(len(detArray)):

		curR = r + (detArray[i].obPlane(obDist)*sin(angle))
		rPlot.append(curR)
		

		sig = detArray[i].signal
		lSig.append(sig)

		vNoise  = detArray[i].cSignal
		cSig.append(vNoise)


	plt.errorbar(rPlot,lSig, yerr = detArray[0].noise,label = 'Ly-a signal')
	plt.plot(rPlot,cSig, label = 'CIII signal')
	plt.ylim(ymin =0.0)
	plt.ylabel('Signal [V]')
	plt.xlabel('R [m]')
	plt.legend()


	#Now create the subplot for the brightness
	p3 = plt.subplot2grid((3,4),(1,2), colspan = 2, rowspan = 1)
	plt.plot([row[0] for row in brightness],[row[1] for row in brightness],label = 'brightness')
	plt.ylabel('Brightness \n [photons/(m^2 sr s)]')
	plt.legend()
	p3.add_patch(patches.Rectangle((rMin+(detArrayHeight/2*cos(angle)),0.1),cov, 1.6*10**21,angle = 0.0 , color = 'red' , alpha = 0.5))

	#Finally the emissitivity
	p4 = plt.subplot2grid((3,4),(0,2), colspan = 2, rowspan = 1)
	plt.ylabel('Emiss \n [photons/(m^3 sr s)]')
	plt.plot([row[0] for row in lineEmiss],[row[2] for row in lineEmiss],label = 'emiss')
	plt.legend()
	p4.add_patch(patches.Rectangle((rMin+(detArrayHeight/2*cos(angle)),0.1),cov, 10**21,angle = 0.0 , color = 'red' , alpha = 0.5))
	fig1.suptitle('gain: '+ str(10**7)+' [V/A] mirror angle: 45 filter shift: +2 [nm] obDist: '+ str(obDist) +'[m] res: '+str(obDist*detArray[0].width/detArray[0].y)+ '[m]')
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

	rMin = round(r - (sin(angle) * length/2),4) #farthest coordinates of detector object plane
	zMin = round(z - (cos(angle) * length/2),4)
	rMax = boundary[1]

	print('rMax: '+ str(rMax))
	print('rMin: ' + str(rMin))
	print('step: ' + str(length/numStep))

	# Here we are calculating the object plane slope
	slope = (z - zMin)/(r-rMin)

	print('angle:' + str(slope))

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
			curBright = signalGen3.boxAve(brightness, round(curR - w,4), round(curR + w,4))

		# if it is at an angle, we have to use the axially symmetric 3D abel transfrom
		else:
			lineEmiss = signalGen3.lineData(curR,curZ,0.0,curR,rMax, length/numStep,emiss,eqObj)
			brightness = signalGen3.torAbel(lineEmiss)
			curBright = brightness[0][1]
			# Precisely, we would box average for each detector over a range of Z at each grid 
			# point; however this is really computationally expensive and uncessary 
			# with the sparsness of the SOLPS data

		#curBright = 1*10**20 # for testing geometry effects of detector angles
		print('brigtness' + str(curBright) ) 

		detArray[i].sigCalc(balmer2Lyman*curBright)
		
		
		sig = detArray[i].signal

		print('noise Brightness: ' + str(curBright*noiseFrac))
		detArray[i].carbonCalc(curBright*noiseFrac,noiseWvlngth)
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
		detArray[i].carbonCalc(curBright*noiseFrac,noiseWvlngth)


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
		# carbon and lyman signal is set to 0 to start
		det = Detector(xDist,-imDist,detW, detParam[1], detParam[4], 0.0, sWidth, sHeight,\
	              detParam[5],detParam[6],detParam[7], detParam[8],detParam[9],detParam[10],detParam[11],0.0)

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
def firstNum(line):
	index = 0
	for ch in line:
		if ch != ' ':
			break
		index += 1
	return index

def myRead(filename):
	data = []

	with open(filename) as f:
		for line in f:
			marker = " 0 "
			index = line.rfind(marker) #this finds the last instance of marker

			cutLine = line.strip('\r\n')[index+len(marker):]

			#now we trim off the remaining spaces
			datum = float( cutLine[firstNum(cutLine):] )
			data.append(datum)

	return data

def main():
	lyALine = 0.000121567 # wavelength in mm
	eLyman = 1.63*10**(-18) # energy of lyman photon J
	mmTm = 1*10**(-6) # correction for units of mm^2 to m^2

	detWidth = 0.75 
	detHeight = 4.05 ## in mm
	detGap = 0.20 #Gap between detectors in array
	nDetectors = 20 # number of detectors
	detLength = nDetectors*detWidth+(nDetectors-1)*detGap
	responseAXUV = 0.15 #AXUV respone
	#"""

	filterTrans = 0.05  #We could get this in 10 though with great FWHM
	filterFWHM = 0.00001 # mm, must match units of center
	filterN = 2.6 #filter effective refractive index
	balmerBright = 5*10**18 #taken directly from matlab (ph/s sr m2)
	lymanBright = 20 * balmerBright

	voltNoise = 20*10**(-3)
	gain = 1*10**7 #set by pre-amp V/A
	unit= 10**(-3)
	"""
	lEmiss = myRead("/home/aaronmr/Desktop/cr162944/162944_ly_D_320.txt")		
	rPos = myRead("/home/aaronmr/Desktop/cr162944/cr_162944_2.txt")
	zPos = myRead("/home/aaronmr/Desktop/cr162944/cz_162944_2.txt")

	plt.scatter( [row[0] for row in emiss],[row[1] for row in emiss], \
		c = [row[2] for row in emiss])
	plt.show()
	"""

	cEmiss = [0.2902*10**17,0.5502*10**17,0.5839*10**17,1.4203*10**17,1.1747*10**17,1.5960*10**17,0.8250*10**17,1.0448*10**17]
	rPos = [2.2033,2.2281,2.2440,2.2598,2.2748,2.2895,2.3032,2.3165]

	emiss = []
	for i in range(len(cEmiss)):
		emiss.append([ rPos[i],0.0,cEmiss[i] ])
	"""[detWidth*unit,detHeight*unit,nDetectors, detGap * unit , responseAXUV,filterTrans,lyALine*unit,\
	filterFWHM*unit,gain,filterN,voltNoise]"""

	detParam = [detWidth*unit,detHeight*unit,nDetectors, detGap * unit , 0.2 ,1.0,465*10**(-9),\
	filterFWHM*unit,gain,filterN,voltNoise]

	r = 2.26
	z =0.0
	obDist = 1.0
	res = 0.005
	angle = pi/2


	detW = detParam[0] 
	imDist = obDist * detW/res
	detArray = []
	nDet = detParam[2]
	print('imDist:' +str(imDist))
	# We will assume the pinhole marks x-center for the detector coordinates
	# The detector coordinates are not in tokomak space, see detector class
	# for definition of coordinates for Detectors.

	for i in range(nDet):

		#This formula is a bit confusing but gives correct result for even or odd detector
		#arrays, with center of array at (0,imDist)
		xDist = (-(nDet-1)/2.0 + i) * (detW + detParam[3]) 
		
		# initialize them with the values
		det = Detector(xDist,-imDist,detW, detParam[1], detParam[4], 0.0, 0.002, 0.008,\
	              detParam[5],detParam[6],detParam[7], detParam[8],detParam[9],detParam[10])

		detArray.append(det)

	grid = signalGen3.genGrid(2.3,0.0,0.0,0.001,2.2,2.4)
	
	lineData =[]
	for i in range(len(grid)):
		asd = signalGen3.closestPoint(grid[i][0],grid[i][1], emiss, 1.0)[2]
		print(asd)
		lineData.append([grid[i][0],grid[i][1],asd])
	
	invert = signalGen3.torAbel(lineData)
	
	
	p1 = plt.subplot2grid((1,3),(0,0))
	p1.plot([row[0] for row in emiss],[row[2] for row in emiss],label = 'emiss')
	plt.xlabel('R [m]')
	plt.ylabel('Emissivity \n [photons/(m^3 sr s)]')
	plt.title('Emissivity')
	p2 = plt.subplot2grid((1,3),(0,1))
	p2.plot([row[0] for row in invert],[row[1] for row in invert],label = 'brightness')
	plt.xlabel('R [m]')
	plt.ylabel('Brightness \n [photons/(m^2 sr s)]')
	plt.title('Brightness')
	p3 = plt.subplot2grid((1,3),(0,2))
	plt.xlabel('R [m]')
	plt.ylabel('Signal [v]')
	plt.title('Det signal')


	for i in range(len(detArray)):

		curR = r + (detArray[i].obPlane(obDist)*sin(angle))
		curZ = z + (detArray[i].obPlane(obDist) * cos(angle)) # if horizontal all the same z
		w = sin(angle)*detArray[i].correctSize(obDist, detArray[i].sWidth, detArray[i].width)/2

		print('curR: '+ str(curR))

		curBright = signalGen3.boxAve(invert, curR - w, curR + w)


		detArray[i].sigCalc(balmer2Lyman*curBright)
		
		sig = detArray[i].signal
		print('signal: ' +str(sig))

		p3.bar(curR-(w*2), sig,w*2)	
		p3.text(curR,0,str(i))


	plt.show()

	return 0

if __name__ == "__main__":
	main()
