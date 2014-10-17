#!/usr/bin/python

# DEPENDENCIES: PYDICOM, SCIPY, MATPLOTLIB

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

def getROI(a,boxsize=100,threshold=100):
	# Get region of interest around center of mass of ndarray object
	# Returns cropped array that is the region of interest
	halflen = boxsize/2
	CM = getCM(a,threshold)
	xmin = CM[0]-halflen
	xmax = CM[0]+halflen
	ymin = CM[1]-halflen
	ymax = CM[1]+halflen
	return a[xmin:xmax,ymin:ymax]

def avgROI(a,subset=(61,242),boxsize=100,threshold=100):
	# Get average image around wire by averaging using getROIs
	# Return array that is the average image
	zmin=subset[0]
	zmax=subset[1]
	zlen=subset[1]-subset[0]+1
	avg = 0
	acropped = a[:,210:330,130:350]
	for z in range(zmin,zmax):
		avg = avg + getROI(acropped[z,:,:],boxsize,threshold)
		print "Finished slice ", z-zmin+1, "/", zmax-zmin, "."
	return avg/zlen

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
	# This actually does work 
	asize = a.shape
	xsize = asize[0]
	ysize = asize[1]
	diam = 0.035 # in mm
	pi = 3.1415926536
	b = numpy.zeros(shape=(101,101))
	for x in range(0,xsize):
		for y in range(0,ysize):	
			r = numpy.sqrt((freqspacing*(x-50))**2.0+(freqspacing*(y-50))**2.0)
			b[x,y] = (2*special.jn(1,pi*r*diam))/(pi*r*diam)
			if ( not (x==50 and y==50)):	# This if statement gets rid of the discontinuity at the zero frequency caused by division by zero (BesselJ)
				a[x,y] = (a[x,y]*pi*r*diam)/(2*special.jn(1,pi*r*diam))
	return a

# Path of DICOM file; read it into pydicom image object
path = "/home/yohan/Projects/mbp1024/data/group2/ctlab2/WirePhantomInSlice/1.2.392.200036.9116.2.5.1.48.1221390955.1384925454.66344.dcm"
img = dicom.read_file(path)
print "Dicom file read."

# Tuple with first and last slice in which the wire point is visible
wireslices = (61, 242)

# Pixel spacing in mm
pixspacing = 0.249
freqspacing = 1./(101.*pixspacing)

# Make ndarray object containing CT HU data from pydicom object 
pixbytes = img.PixelData
pixarray = img.pixel_array
print "Pixel values extracted."

# Convert to 32-bit signed integer so that sum in averaging process doesn't exceed bounds
pixarray = pixarray.astype(numpy.int32)
print "Array type converted to 32-bit signed integer."

# Compute average image around wire
print "Starting averaging process..."
avgimg = avgROI(pixarray,subset=wireslices)
print "Averaging process finished."

# Compute PSF by killing background below threshold
psf = killBG(avgimg,threshold=67)
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
pylab.imshow(psf, cmap=pylab.cm.bone).set_extent([-50*pixspacing,50*pixspacing,-50*pixspacing,50*pixspacing])
plt.colorbar()
plt.xlabel("Distance in x (mm)",fontdict={'fontsize':16})
plt.ylabel("Distance in y (mm)",fontdict={'fontsize':16})
plt.xticks( numpy.arange(-12,12,2) )
plt.yticks( numpy.arange(-12,12,2) )
plt.title("Point Spread Function")
pylab.show()
"""
# Plot 2D MTF

#pylab.imshow(mtf, cmap=pylab.cm.bone).set_extent([-101./(2*pixspacing),101./(2*pixspacing),-101./(2*pixspacing),101./(2*pixspacing)])
pylab.imshow(mtf, cmap=pylab.cm.bone).set_extent([-1./(2*pixspacing),1./(2*pixspacing),-1./(2*pixspacing),1./(2*pixspacing)])
plt.colorbar()
plt.xlabel("Frequency in x (inverse mm)",fontdict={'fontsize':16})
plt.ylabel("Frequency in y (inverse mm)",fontdict={'fontsize':16})
#plt.xticks( numpy.arange(-200,200, 50 ) )
#plt.yticks( numpy.arange(-200,200, 50 ) )
plt.xticks( numpy.arange(-2,2, 0.50 ) )
plt.yticks( numpy.arange(-2,2, 0.50 ) )
plt.title("Modulation Transfer Function")
pylab.show()


#OLD STUFF: 1D PSF/MTF CODE

# Splice MTF into horizontal and vertical slices centered at zero value
#hslice = rawmtf[:,50]
#vslice = rawmtf[50,:]

# Plot 1D
#plt.plot(hslice)
#plt.show()


	
