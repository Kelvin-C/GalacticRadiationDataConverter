from mayavi import mlab
import numpy as np
import healpy as hp
import matplotlib.pyplot as plt

plt.close()

#load a fits file
input_filename = "COM_CompMap_Synchrotron-commander_0256_R2.00.fits"

#intmap is the first column of hdu1 i.e. is intensity map
intmap = hp.read_map(input_filename) # * 1e3 # muK
nside = hp.npix2nside(len(intmap))

#no. of bins
nbr_bins=300 

#compute histogram of data; normed converts data to pdf so overall integral is 1
#freq is the no. of data points in the bin
#bins_int_dpts are the bins interval data points
freq, bins_int_dpts = np.histogram(intmap, nbr_bins ,normed=False) 

#cumulative distribution function
#cumfreq is cumulative frequency
cumfreq = freq.cumsum() 
#normalise to pdf
#cumprob = cumfreq / cumfreq[-1]

#use linear interpolation of cdf to find new pixel values
im2 = np.interp(intmap,bins_int_dpts[:-1],cumfreq)

#print im2.reshape(intmap.shape), cumprob


#plt histogram
plt.subplot(3, 1, 1)
plt.hist(intmap, nbr_bins, normed=True)#, range = (0,0.1e9))
plt.title("Histogram")
plt.xlabel("Intensity")
plt.ylabel("Frequency")
#plt.axis([0, 0.01e8, 0, 0.000005])


#plt pdf
plt.subplot(3, 1, 2)
plt.hist(intmap, nbr_bins, normed=True, cumulative=True)
plt.title("Cumulative probability")
plt.xlabel("Intensity")
plt.ylabel("Cumulative probability")

#plt after-histogram
plt.subplot(3, 1, 3)
plt.hist(im2, nbr_bins, normed=True)#, range = (0,0.1e9))
plt.title("Histogram")
plt.xlabel("Intensity")
plt.ylabel("Frequency")
#plt.axis([0, 0.01e8, 0, 0.000005])


plt.tight_layout()
plt.show()


