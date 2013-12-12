#!/usr/bin/python 

# DEPENDENCIES: PYDICOM, SCIPY etc..

import dicom, numpy, pylab, sys
from scipy import ndimage
from scipy import special
import matplotlib.pyplot as plt

def getCM(a,threshold=100):
	# Get the center of mass of an ndarray object
	# Returns a tuple that is the CM pixel coordinate
	M = 0
	xc = 0
	yc = 0
	asize = a.shape
	xsize = asize[0]
	ysize = asize[1]
	for x in range(0,xsize):
		for y in range(0,ysize):
			if a[x,y] >= threshold:
				M = M+1
				xc = xc + x
				yc = yc + y
	return (xc/M,yc/M)

def getROI(a,boxsize=60,threshold=100,debug=False):
	# Get region of interest around center of mass of ndarray object
	# Returns cropped array that is the region of interest
	halflen = boxsize/2
	CM = getCM(a,threshold)
	xmin = CM[0]-halflen
	xmax = CM[0]+halflen
	ymin = CM[1]-halflen
	ymax = CM[1]+halflen
	if (debug == True):
		#Debug
		print CM
		pylab.imshow(a[xmin:xmax,ymin:ymax], cmap=pylab.cm.bone)
		plt.colorbar()
		pylab.show()
		# End debug
	return a[xmin:xmax,ymin:ymax]

def avgROI(a,subset=(80,110),boxsize=60,threshold=100):
	# Get average image around wire by averaging using getROIs
	# Return array that is the average image
	ymin=subset[0]
	ymax=subset[1]
	ylen=subset[1]-subset[0]+1
	avg = 0
	acropped = a[124:186,:,220:340]
	for y in range(ymin,ymax):
		if (y==120):
			db = True
		else:
			db = False
		avg = avg + getROI(acropped[:,y,:],boxsize,threshold,debug=db)
		print "Finished slice ", y-ymin+1, "/", ymax-ymin, "."
	return avg/ylen

def killBG(a,threshold=65):
	# Kill background by setting all values under threshold to threshold, and then "zero" the background	
	# Return array that is the image without unwanted background
	afg = a
	afgsize = afg.shape
	xsize = afgsize[0]
	ysize = afgsize[1]
	for x in range(0,xsize):
		for y in range(0,ysize):
			if afg[x,y] <= threshold:
				afg[x,y] = threshold
	return afg-threshold

def MTFcorrect(a):
	# Correct MTF according to Bischof (1977)
	# DOESNT WORK, FIX THIS! UNIT SCALING NEEDS TO BE DONE
	asize = a.shape
	xsize = asize[0]
	zsize = asize[1]
	diam = 0.035 # in mm
	pi = 3.1415926536
	b = numpy.zeros(shape=(61,61))
	for x in range(0,xsize):
		for z in range(0,zsize):	
			r = numpy.sqrt((freqspacingx*(x-30))**2.0+(freqspacingz*(z-30))**2.0)
			b[x,z] = (2*special.jn(1,pi*r*diam))/(pi*r*diam)
			if ( not (x==30 and z==30)):	# This if statement gets rid of the discontinuity at the zero frequency caused by division by zero (BesselJ)
				a[x,z] = (a[x,z]*pi*r*diam)/(2*special.jn(1,pi*r*diam))
	return a

# Path of DICOM file; read it into pydicom image object
path = "/home/yohan/Projects/mbp1024/data/othergroup/ctlab2/WirePhantomCrossSlice/group4.dcm"
img = dicom.read_file(path)
print "Dicom file read."

# Tuple with first and last slice in which the wire point is visible
wireslices = (80, 95) #Make sure wire is visible in this range, i.e. SNR is high or getCM doesn't find wire center and you get a division by zero error

# Pixel spacing in mm
pixspacingx = 0.395
pixspacingz = 0.5
freqspacingx = 1./(61.*pixspacingx)
freqspacingz = 1./(61.*pixspacingz)

# Make ndarray object containing CT HU data from pydicom object 
pixbytes = img.PixelData
pixarray = img.pixel_array
print "Pixel values extracted."

# Convert to 32-bit signed integer so that sum in averaging process doesn't exceed bounds
pixarray = pixarray.astype(numpy.int32)
print "Array type converted to 32-bit signed integer."

# Compute average image around wire
print "Starting averaging process..."
avgimg = avgROI(pixarray,subset=wireslices,threshold=100)
print "Averaging process finished."

# Compute PSF by killing background below threshold
psf = killBG(avgimg,threshold=69)
print "Background noise removed. PSF computed."

# Compute FFT of image
fftimg = abs(numpy.fft.fft2(psf))
print "Fourier transform of PSF computed."

# Compute "Raw" MTF of image by shifting and normalizing FFT
# old normalization: rawmtf = 2.*numpy.fft.fftshift(fftimg) / (fftimg.shape[0]**2.)
rawmtf = numpy.fft.fftshift(fftimg)/fftimg[0,0]
print "Raw MTF computed."

# Compute corrected MTF
mtf = MTFcorrect(rawmtf)
print "Corrected MTF computed."

# Plot 2D PSF
# UNCOMMENT TO DISPLAY PSF
"""
pylab.imshow(psf, cmap=pylab.cm.bone).set_extent([-30*pixspacingx,30*pixspacingx,-30*pixspacingz,30*pixspacingz])
plt.colorbar()
plt.xlabel("Distance in x (mm)",fontdict={'fontsize':16})
plt.ylabel("Distance in z (mm)",fontdict={'fontsize':16})
plt.xticks( numpy.arange(-10,10,5) )
plt.yticks( numpy.arange(-15,15,5) )
plt.title("Point Spread Function")
pylab.show()
"""
# Plot 2D MTF

#pylab.imshow(mtf, cmap=pylab.cm.bone).set_extent([-61./(2*pixspacingx),61./(2*pixspacingx),-61./(2*pixspacingz),61./(2*pixspacingz)])
pylab.imshow(mtf, cmap=pylab.cm.bone).set_extent([-1./(2*pixspacingx),1./(2*pixspacingx),-1./(2*pixspacingz),1./(2*pixspacingz)])
plt.colorbar()
plt.xlabel("Frequency in x (inverse mm)",fontdict={'fontsize':16})
plt.ylabel("Frequency in z (inverse mm)",fontdict={'fontsize':16})
#plt.xticks( numpy.arange(-70,70, 10 ) )
#plt.yticks( numpy.arange(-60,60, 10 ) )
plt.xticks( numpy.arange(-1,1, 0.5 ) )
plt.yticks( numpy.arange(-1,1, 0.5 ) )
plt.title("Modulation Transfer Function")
pylab.show()


#OLD STUFF: 1D PSF/MTF CODE

# Splice MTF into horizontal and vertical slices centered at zero value
#hslice = rawmtf[:,50]
#vslice = rawmtf[50,:]

# Plot 1D
#plt.plot(hslice)
#plt.show()


	
