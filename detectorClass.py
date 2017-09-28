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
	Each etector is associated with a filter and slit so store those variables
	sWidth = width of the slit 
	sHeight = height of the slit
	trans = transmission of filter
	wlngth = wavelength of center of filter
	FWHM = full width half maximum of filter
	gain = gain of detector
	n = index of refraction of detector 
	"""


	def __init__(self, x, y, width, height, response, signal , sWidth, sHeight, \
		         trans, wlngth , FWHM, gain, n, noise):
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



	"""returns the angle as measured from the pinhol normal"""
	def angle(self):
		return atan(self.x/self.y)

	"""Takes in the center wavelength and index of refraction of filter
	and returns the corresponding cwl"""
	def cwl(self):
		theta = self.angle()
		if (theta > pi/12):
			print('Warning: Incident Light on Filter Exceeds Limits')


		#this equation assumes that the index of refraction around teh filter
		# is 1, probably a pretty goood estimate, calculates the shift in the 
		#central wavelength of the filter
		centwave =  self.wlngth*sqrt(1.0-(abs(sin(theta))/self.n)**2)
		

		# We now assume that the filter roughtly behaves the same as central propogation
		# with the central wavelength shiting. This is probably NOT completely true
		# however it does have the benefit of decreasing transmission as the filter is
		# angled in a reasonable manner, so why not? 

		# use wavelength as what we are absorbing
		# decrease the amplitude as angle increases linearly to 4% at 15 degrees
		# this is rough number provided in email.
		transmission = self.calcTrans()
		return createGaussian(self.FWHM, centwave, transmission, self.wlngth)


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
		# real one power = brightness * sG * self.cwl()* eLyman #power deposited)

		# use overestimate for now
		power = brightness * sG * self.trans* eLyman #power deposited)

		i = power * self.response

		self.signal = i * self.gain

		return self.signal

	# Update to something more accurate - AR Aug 2, 2017
	def calcTrans(self):
		return self.trans-abs(self.angle())/(pi/12)*0.01

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
