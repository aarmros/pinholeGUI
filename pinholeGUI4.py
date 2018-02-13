from math import *
# Tkinter changed capitlization betwen python 2 and 3 hence
# the annoying import.
try:
	from Tkinter import *
except ImportError:
	from tkinter import *

import pinholeDet4
import eqtools
import signalGen3
import numpy as np


""" Sep 22, 2017 - No updates to pinholeGUI3 other then calling updated signalGen2 
                   and pinholeDet3
	Sep 30, 2017 - Data loading now has its own part in GUI, all hapens in this file
	               instead of in signalGen. Cuts down rundtime and calls to eqtools
	             - Minor changes to the GUI interface visually"""

# All measurements given in mm, etendue correction assumes this
lyALine = 0.000121567 # wavelength in mm
eLyman = 1.63*10**(-18) # energy of lyman photon J
mmTm = 1*10**(-6) # correction for units of mm^2 to m^2


# We will assume that the detector quantities are fixed 
# Two detector types  

#for 16 array
"""
detWidth = 2.04 
detHeight = 5 ## in mm
detGap = 0.12 #Gap between detectors in array
nDetectors = 16 # number of detectors
detLength = nDetectors*detWidth+(nDetectors-1)*detGap
responseAXUV = 0.15 #AXUV respone
"""

# For 20 array
#"""

detWidth = 0.75 
detHeight = 4.05 ## in mm
detGap = 0.20 #Gap between detectors in array
nDetectors = 20 # number of detectors
detLength = nDetectors*detWidth+(nDetectors-1)*detGap
responseAXUV = 0.15 #AXUV respone
#"""


filterFWHM = 0.00001 # mm, must match units of center
filterN = 2.6 #filter effective refractive index
balmerBright = 5*10**18 #taken directly from matlab (ph/s sr m2)
lymanBright = 20 * balmerBright

voltNoise = 20*10**(-3)
gain = 1*10**7 #set by pre-amp V/A

#Filter Transmission info
# input as [[wavlength0, trans0],...] with trans in percent below is manufacturer specified
"""filterTrans = [[1.1500000000000001e-07, 0.656443754], [1.16e-07, 1.640656909], [1.17e-07, 2.73955072], [1.1800000000000001e-07, 5.411452932],\
 [1.1900000000000001e-07, 7.609680084], [1.2000000000000002e-07, 7.736856405], [1.21e-07, 7.470670058], [1.22e-07, 7.304309542],\
  [1.23e-07, 6.642470524], [1.24e-07, 5.792063537], [1.2500000000000002e-07, 4.893774912], [1.2600000000000002e-07, 4.035165991], [464, 1.0], [466, 1.0]]"""
# this gives +2 nm
filterTrans = [[1.1500000000000001e-07, 0.0],[1.1600000000000001e-07, 0.0],[1.1700000000000001e-07, 0.656443754], [1.18e-07, 1.640656909], [1.19e-07, 2.73955072], [1.2000000000000001e-07, 5.411452932],\
 [1.2100000000000001e-07, 7.609680084], [1.2200000000000002e-07, 7.736856405], [1.23e-07, 7.470670058], [1.24e-07, 7.304309542],\
  [1.25e-07, 6.642470524], [1.26e-07, 5.792063537], [1.2700000000000002e-07, 4.893774912], [1.2800000000000002e-07, 4.035165991], [464, 1.0], [466, 1.0]]
mirrorT = [26.75,48.13] # transmission at [CIII 117.5, lya 121.5] 




"""Helper funciton finds the index of the first non blank charcter in a str"""
def firstNum(line):
	index = 0
	for ch in line:
		if ch != ' ':
			break
		index += 1
	return index

"""using for reading in files created by SOLPS which have a very wierd formatting.
   This was based off of one SOLPS data file so might need to be changed."""

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


# Calculate the finite slit size correction
# Might be useful to see if diffraction is important
# shouldn't be but could check it here
def correctSize(obDist,imDist,pinWidth,detSize):
    correction = ((obDist + imDist)* pinWidth/obDist)
    magnification = ((obDist/imDist)*detSize)
    
    diffract = 1.22*imDist*lyALine/pinWidth
    if 10*diffract > correction**(1/2):
        print ("Diffraciton limited")

        
    #should this be added as squares or just added- big difference    
    # (correction**2+magnification**2+diffract**2)**(1/2)
    return (diffract+correction+magnification)

# updates the calculated variables formats for display as strings
# this is a workaround because tkinter doesn't allow you to do it 
# directly. So some global variables need to be defined which are 
# undesirable. 

#This is sloppy but works
def reFormat():


	magStr.set(str('%.2e' % mag.get()))
	imDistStr.set(str('%.2e' % imDist.get()))
	overlapStr.set(str('%.1f' % overlap.get()))
	eWidthStr.set(str('%.2e' % eWidth.get()))
	eHeightStr.set(str('%.2e' %eHeight.get()))
	covStr.set(str('%.2e' %cov.get()))
	powerStr.set(str('%.2e' % power.get()))
	signalStr.set(str('%.2e' %signal.get()))
	signalIStr.set(str('%.2e' %signalI.get()))
	sig2NoiStr.set(str('%.2f' % sig2Noi.get()))
	angIncStr.set(str('%.1f' % angInc.get()))

	return 0

def updateMirr():
	if mirror.get() == True:
		cTrans.set(mirrorT[0])
		lTrans.set(mirrorT[1])
	else:
		cTrans.set(100.0) # if no mirror, mirror transmission is 100%
		lTrans.set(100.0)

def updateAng():
	if tor.get() == True:
		ang.set(90)
		rSig.set(2.24)
		zSig.set(-0.25)
	else:
		ang.set(0.0)
		rSig.set(1.6)
		zSig.set(1.0)


#this is where all the updating happens whenever one of the sliders
# or entry is changed. trace far below calls this funciton
def MyUpdate(variable,b,c):


	imDist.set(obDist.get() * detWidth / resolution.get())
	print('image dist'+str(imDist.get()))

	mag.set(obDist.get()/ imDist.get())

	# Object size correction
	eWidth.set( correctSize(obDist.get(),imDist.get(),sWidth.get(), detWidth))
	Mtot = eWidth.get()/ detWidth

	DZ =mag.get()*detHeight
	eHeight.set( correctSize(obDist.get(),imDist.get(),sHeight.get(),detHeight))

	#Now we will calculate etendue using solid angle * area
	dA = detWidth*detHeight
	DA = resolution.get()*DZ #this is calculated slightly different from matlab script
	sA = sWidth.get()*sHeight.get()

	omegas = sA/(imDist.get()**2) # solid angle of slit from center of detector
	omegaD = DA/(obDist.get()**2) # solid angle of ojbect from pinhole center

	sG = omegas*dA*mmTm
	DG = omegaD*sA*mmTm	


	power.set(lymanBright * sG * mirrorT[1]*eLyman) #power deposited)

	signalI.set(power.get()*responseAXUV)

	signal.set(gain * signalI.get())
	sig2Noi.set(signal.get()/voltNoise)

	gapRes = (detGap+detWidth/2)/imDist.get() * obDist.get() - resolution.get()/2
	cov.set(nDetectors * resolution.get() + (nDetectors-1)*gapRes)
	
	#I'm a little unsure of these factors of 2 -AR July 31 2017
	overlap.set(100*(eWidth.get() - resolution.get()-gapRes*2)/(2 * resolution.get()))

	# We change the string to degrees, cause no matter what you say
	# it is easier to understand quickly
	#calculate maximum angle of incidence using detector properties and image distance

	angInc.set(degrees(atan(0.5*detLength/imDist.get())))

	reFormat()
	return 0

def pinholeSim():
	print('len is '+ str(len(emiss)))
	sigZ = zSig.get() #these two are in meters
	sigR = rSig.get()

	unit= 10**(-3) #this file is in mm, everything else in m
	#MAKE SURE TO ONLY PASS THINGS IN Meters TO UNDERLYING FILES
	mirrorT = [cTrans.get(),lTrans.get()]


	detParam = [detWidth*unit,detHeight*unit,nDetectors, detGap * unit , responseAXUV,filterTrans,lyALine*unit,\
	filterFWHM*unit,gain,filterN,voltNoise,mirrorT]

	#fileNames = [brightFile.get(),rFile.get(),zFile.get(),agFile.get()]

	
	pinholeDet4.runSim(sigZ,sigR, radians(ang.get()),obDist.get()*unit,resolution.get()*unit,\
		sWidth.get()*unit,sHeight.get()*unit,tor.get(),detParam,emiss,eq)

	root.lift()
	root.mainloop()
	
	return 0

def loadData():

	#magnetic reconstruction for flux interpolation
	global eq
	eq = eqtools.EqdskReader(gfile = 'g' + str(agFile.get()),afile = 'a'+ str(agFile.get()))


	lEmiss = myRead(brightFile.get())		
	rPos = myRead(rFile.get())
	zPos = myRead(zFile.get())
	psiPos = eq.rz2psinorm(np.array(rPos),np.array(zPos))[0] #eqtools uses numpy arrays
	angPos = signalGen3.magAngle(rPos,zPos,eq)


	#Organize it for how the functions are written
	global emiss
	for i in range(len(lEmiss)):
		emiss.append([ rPos[i],zPos[i],lEmiss[i],psiPos[i],angPos[i]])
		#emiss.append([ rPos[i],zPos[i],1.0*10**20,psiPos[i],angPos[i]]) # for testing geometry effects

	print('Data Loaded')

	root.lift()
	root.mainloop()

	return 0



#This is all for setting up the GUI, it takes a lot of code :(
root = Tk()
root.geometry('1100x450+350+70')


# Resolution GUI elements

resLabel = Label(root, text = 'Ideal Resolution')
resLabel.grid(row = 0, column = 0)
resolution = DoubleVar()
resolution.set(5)
resEntry = Entry (root,text = resolution, bg = 'white')
resEntry.grid(row = 0, column = 1)
resSlider = Scale(root, variable = resolution, orient = 'horizontal',from_ = 0, to = 10)
resSlider.grid(row=2, column = 1)

#object distance GUI elements

obLabel = Label(root, text = 'Object Distance',height = 2)
obLabel.grid(row = 3, column = 0)
obDist= DoubleVar()
obDist.set(1000)
obEntry = Entry (root,textvariable = obDist, bg = 'white')
obEntry.grid(row = 3, column = 1)
obSlider = Scale(root, variable = obDist, orient = 'horizontal',from_ = 5, to = 2000)
obSlider.grid(row=4, column = 1)

#define slit Height GUI elements

sHLabel = Label(root, text = 'Slit Height',height = 2)
sHLabel.grid(row = 5, column = 0)
sHeight= DoubleVar()
sHeight.set(8)
sHEntry = Entry (root,textvariable = sHeight, bg = 'white')
sHEntry.grid(row = 5, column = 1)
sHSlider = Scale(root, variable = sHeight, orient = 'horizontal',from_ = 0, to = 20)
sHSlider.grid(row=6, column = 1)

# define slit Width GUI elements

sWLabel = Label(root, text = 'Slit Width ',height = 2)
sWLabel.grid(row = 7, column = 0)
sWidth= DoubleVar()
sWidth.set(2)
sWEntry = Entry (root,textvariable = sWidth, bg = 'white')
sWEntry.grid(row = 7, column = 1)
sWSlider = Scale(root, variable = sWidth, orient = 'horizontal',from_ = 0, to = 10)
sWSlider.grid(row=8, column = 1)


rLabel = Label(root, text = 'Center of Object Plane r(m):')
rLabel.grid(row = 10, column = 0)
rSig = DoubleVar()
rSig.set(1.91)
rEntry = Entry (root,text = rSig, bg = 'white')
rEntry.grid(row = 10, column = 1)

zLabel = Label(root, text = 'Center of Object Plane z(m):')
zLabel.grid(row = 11, column = 0)
zSig = DoubleVar()
zSig.set(-0.72)
zEntry = Entry (root,text = zSig, bg = 'white')
zEntry.grid(row = 11, column = 1)

angLabel = Label(root, text = 'Detector Orientation from vertical (degrees):')
angLabel.grid(row = 12, column = 0)
ang = DoubleVar()
ang.set(90)
angEntry = Entry (root,text = ang, bg = 'white')
angEntry.grid(row = 12, column = 1)

tor = BooleanVar()
tor.set(True)
torButton = Checkbutton(root, text = 'Toroidal View [unchecked poloidal]', variable = tor,
	onvalue = True, offvalue = False, command = lambda: updateAng())
torButton.grid(row = 14, column = 0)

mirror = BooleanVar()
mirror.set(True)
mirrorButton = Checkbutton(root, text = 'Include Mirror', variable = mirror,
	onvalue = True, offvalue = False,command = lambda: updateMirr())
mirrorButton.grid(row = 15, column = 0)

cTrans = DoubleVar()
cTrans.set(mirrorT[0])
lTrans = DoubleVar()
lTrans.set(mirrorT[1])


button = Button(root,text = "Run", command = lambda: pinholeSim())
button.grid(row = 14, column = 1)


# this updates are calculation anything is changed
obDist.trace('w',MyUpdate)
resolution.trace('w',MyUpdate)
sWidth.trace('w',MyUpdate)
sHeight.trace('w',MyUpdate)


#Now for the outputs of this GUI
imText = Label(root, text = ' Ideal Image Distance ')
imDist = DoubleVar()
imDistStr = StringVar()
imLab = Label(root, textvariable = imDistStr)
imText.grid(row = 3, column = 2)
imLab.grid(row = 3, column = 3)

magText = Label(root, text = ' Ideal Magnification ')
mag = DoubleVar()
magStr = StringVar()
magLab = Label(root, textvariable = magStr)
magText.grid(row = 0, column = 2)
magLab.grid(row = 0, column = 3)

covText = Label(root, text = ' Horizontal Coverage ')
cov = DoubleVar()
covStr = StringVar()
covLab = Label(root, textvariable = covStr)
covText.grid(row = 7, column = 2)
covLab.grid(row = 7, column = 3)

eWText = Label(root, text = ' Effective Object Width ')
eWidth = DoubleVar()
eWidthStr = StringVar()
eWLab = Label(root, textvariable = eWidthStr)
eWText.grid(row = 5, column = 2)
eWLab.grid(row = 5, column = 3)

eHText = Label(root, text = ' Effective Object Height')
eHeight = DoubleVar()
eHeightStr = StringVar()
eHLab = Label(root, textvariable = eHeightStr)
eHText.grid(row = 6, column = 2)
eHLab.grid(row = 6, column = 3)

overlapText = Label(root, text = ' Measurement Overlap (%)')
overlap = DoubleVar()
overlapStr = StringVar()
overlapLab = Label(root, textvariable = overlapStr)
overlapText.grid(row = 4, column = 2)
overlapLab.grid(row = 4, column = 3)

filtLab = Label(root, text = 'Filter Transmission')
filtLab.grid(row = 0, column = 4)
filtLab2 = Label(root, text = str(filterTrans))
filtLab2.grid(row = 0, column = 5)

lymLab = Label(root, text = 'Lyman Brightness (ph/s sr m2)')
lymLab2 = Label(root, text = str('%.2e' % lymanBright))
lymLab.grid(row = 1, column = 4)
lymLab2.grid(row = 1, column = 5)

gainLab = Label(root, text = 'Gain (V/A)')
gainLab.grid(row = 2, column = 4)
gainLab2 = Label(root, text = str('%.2e' % gain))
gainLab2.grid(row = 2, column = 5)

vNoiseLab = Label(root, text = 'Voltage Noise (V)' )
vNoiseLab.grid(row = 3, column = 4)
vNoiseLab2 = Label(root, text = str(voltNoise))
vNoiseLab2.grid(row = 3, column = 5)

powerText = Label(root, text = ' Power Deposited (W) ')
power = DoubleVar()
powerStr = StringVar()
powerLab = Label(root, textvariable = powerStr)
powerText.grid(row = 4, column = 4)
powerLab.grid(row = 4, column = 5)

sigText= Label(root, text = ' Signal Voltage (V) ')
signal = DoubleVar()
signalStr = StringVar()
sigLab = Label(root, textvariable = signalStr)
sigText.grid(row = 5, column = 4)
sigLab.grid(row = 5, column = 5)

sigIText= Label(root, text = ' Signal Current (A)')
signalI = DoubleVar()
signalIStr = StringVar()
sigILab = Label(root, textvariable= signalIStr)
sigIText.grid(row = 6, column = 4)
sigILab.grid(row = 6, column = 5)

s2NText = Label(root, text = ' Signal to Noise (V /V)')
sig2Noi = DoubleVar()
sig2NoiStr = StringVar()
s2NLab = Label(root, textvariable= sig2NoiStr)
s2NText.grid(row = 7, column = 4)
s2NLab.grid(row = 7, column = 5)

angIncText = Label(root, text='Maximum Angle of Incidence (deg)')
angInc = DoubleVar()
angIncStr = StringVar()
angIncLab = Label(root, textvariable = angIncStr)
angIncText.grid(row=8, column = 2)
angIncLab.grid(row  = 8, column = 3)

info = Label(root, text = 'ALL DISTANCES ARE IN mm')
info.grid(row=8, column = 4 )

brightLabel = Label(root, text = 'Brightness File Name:')
brightLabel.grid(row = 10, column = 2)
brightFile = StringVar()
brightFile.set("SOLPS/lyman_alpha_data.txt")
brightEntry = Entry (root,text = brightFile, bg = 'white')
brightEntry.grid(row = 10, column = 3)

rLabel = Label(root, text = 'R File Name:')
rLabel.grid(row = 11, column = 2)
rFile = StringVar()
rFile.set("SOLPS/lyman_alpha_R.txt")
brightEntry = Entry (root,text = rFile, bg = 'white')
brightEntry.grid(row = 11, column = 3)

zLabel = Label(root, text = 'Z file Name:')
zLabel.grid(row = 12, column = 2)
zFile = StringVar()
zFile.set("SOLPS/lyman_alpha_Z.txt")
zEntry = Entry (root,text = zFile, bg = 'white')
zEntry.grid(row = 12, column = 3)

agLabel = Label(root, text = 'a/g file Ending:')
agLabel.grid(row = 13, column = 2)
agFile = StringVar()
agFile.set("156867.02006")
agEntry = Entry (root,text = agFile, bg = 'white')
agEntry.grid(row = 13, column = 3)

emiss = []
#eq = eqtools.EqdskReader(gfile = 'g' + str(agFile.get()),afile = 'a'+ str(agFile.get()))

loadButton = Button(root,text = "Load", command = lambda: loadData())
loadButton.grid(row = 14, column = 3)


#button = Button(root, text = "Update", command = Insert)


#button.grid(row = 5, column = 1)


root.mainloop()

