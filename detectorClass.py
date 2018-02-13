from math import *



h = 6.63*10**(-34) ## plack's constant J s
c = 3.0*10**(8)  ## m/s


class Detector(object):
	"""Class for each detector. I.E. the array is made up of detector classes
	
	all distances should be in m
	Attributes:
	y = Distance from pinhole towards detector (perpendicular to pinhole) (image Distance)
	x = Distance from pinhole center to detector center (in plane of pinhole)
	width =  Width of detector (resolution direction, plasma r major radius)
	height = Height of detector (vertical in height, plasma z)
	response detector response A/W
	signal = signal that detector is recieving, initializes to 0
	Each detector is associated with a filter and slit so store those variables
	sWidth = width of the slit 
	sHeight = height of the slit
	trans = transmission of filter [[lambda0, lambda0 transmission],..]
	wlngth = wavelength of center of filter
	FWHM = full width half maximum of filter
	gain = gain of detector
	n = index of refraction of detector
	mirrorTrans = mirror transmission at CIII and lyman [CIII trans,lya trans] 
	"""


	def __init__(self, x, y, width, height, response, signal , sWidth, sHeight, \
		         trans, wlngth , FWHM, gain, n, noise, mirrorTrans,cSignal):
		self.x = x
		self.y = y
		self.width = width
		self.height = height
		self.response = response
		self.gain = gain
		self.signal = signal
		self.sWidth = sWidth
		self.sHeight = sHeight
		self.trans = trans
		self.wlngth = wlngth
		self.FWHM = FWHM
		self.gain = gain
		self.n = n
		self.noise = noise
		self.mirrorTrans = mirrorTrans
		self.cSignal =cSignal



	"""returns the angle as measured from the pinhole normal"""
	def angle(self):
		return atan(self.x/self.y)

	# Update to something more accurate - AR Aug 2, 2017
	# linear decrease in transmission as filter is tilted 1% decrease over
	# 15 degrees - rough numbers provided by Acton Optics in email
	def calcTrans(self):
		return self.trans-abs(self.angle())/(pi/12)*0.01

	"""Takes in a wavelength wl and ouputs the transmission at that wl
	through the filter. Takes into account shift in transmission curve of filter
	due to angle of light onto filter. DOES NOT take into account decrease
	of transmission overall due to angle of incidence"""
	def filterTrans(self,wl):

		theta = self.angle()
		if (theta > pi/12):
			print('Warning: Incident Light on Filter Exceeds Limits')


		#this equation assumes that the index of refraction around the filter
		#is 1, probably a pretty goood estimate, calculates the shift in the 
		#central wavelength of the filter
		wlShift =  wl*sqrt(1.0-(abs(sin(theta))/self.n)**2)
		

		# We now assume that the filter roughtly behaves the same as normal incidence 
		# with the central wavelength shifting. This is probably NOT completely true

		#Gives transmission at shifted wl without decrease of overall transmission due to shift
		transmission = self.filterTransCurve(wlShift)

		return transmission


	"""returns the x position of the plane being imaged by the detector"""
	def obPlane(self, obDist):
		return tan(self.angle())*obDist

	# Calculate the finite slit size correction
	# Might be useful to see if diffraction is important
	# shouldn't be but could check it here
	def correctSize(self, obDist, sW, dW):

		imDist = abs(self.y)
		correction = ((obDist + imDist)* sW/obDist)
		magnification = ( ( obDist/ imDist ) * dW)
	    
		diffract = 1.22*imDist*self.wlngth/sW
		"""
		if ((10*diffract) > (correction**(1/2))):
	    	#print ("Diffraciton limited")
	    	print('diffraction limited')
		"""
		return (correction+magnification+diffract)

	# this sets and returns the signal voltage for the detector
	def sigCalc(self,brightness):
			#Now we will calculate etendue using solid angle * area
		dA = self.width * self.height
		
		sA = self.sWidth * self.sHeight 

		omegas = sA/(self.y**2) # solid angle of slit from center of detector
	

		sG = omegas * dA 

		eLyman = h * c/self.wlngth
		t = self.filterTrans(self.wlngth)*(self.mirrorTrans[1]*0.01)
		print('Signal Trans: '+ str(t))
		
		power = brightness * sG * t* eLyman #power deposited)

		# use overestimate for now
		#power = brightness * sG * self.trans* eLyman #power deposited)

		i = power * self.response

		self.signal = i * self.gain

		return self.signal
	

	def carbonCalc(self,brightness, carbonWvlngth):

		# this is transmission of noiseWvlngth, factor 0.01 correcting for percent
		nTrans = self.filterTrans(carbonWvlngth)*(self.mirrorTrans[0]*0.01) 
		print('Carbon Trans: '+ str(nTrans))
		
		dA = self.width * self.height
		
		sA = self.sWidth * self.sHeight 

		omegas = sA/(self.y**2) # solid angle of slit from center of detector
	

		sG = omegas * dA
		eNoise = h * c/carbonWvlngth
		
		power = brightness * sG * nTrans* eNoise

		#assumes the response is constant across wavelength which is true
		# if the wvlngth is close to the center wvlngth
		i = power * self.response

		carbonVolt = i * self.gain
		print("carbon Volt: " + str(carbonVolt))

		self.cSignal = carbonVolt

		return carbonVolt

		
	def filterTransCurve(self,curWl):

		i = 0;
		while (curWl >= self.trans[i][0]):
			i += 1
		
		#linear weighted average, convert from percent to transmission fraction with 0.01
		aveTrans = abs(((self.trans[i][1]*0.01*(curWl-self.trans[i-1][0]))+(self.trans[i-1][1]*0.01*(self.trans[i][0]-curWl)))/(self.trans[i][0]-self.trans[i-1][0]))

		return aveTrans
		



	# the distance covered by the array is the distance to the farthest
	# detector object plane (center of detector object plane) + detector object width/2
	def obPlaneLength(self,obDist):
		return 2 * (abs(self.obPlane(obDist)) +\
		self.correctSize(obDist, self.sWidth, self.width)/2)




# THIS creates an overly centralized peak but not a bad fit to the
# data on the website
# function that returns the value of a guassian at x with parameters
# FWHM, center and height
def createGaussian(FWHM, center,height, x):
	sigma = FWHM / (2 * sqrt(2 * log(2)))  # relate sigma to FWHM
	return height * exp(-(x-center)**2/(2 * sigma**2)) #unnormalized

