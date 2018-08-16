from astropy.io import fits
import astropy.table as table
import healpy as hp
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.pyplot import cm
import os
import Tkinter as Tk
import tkFileDialog

def histeq(polarisation_fraction, nbr_bins):
    #compute histogram of data; normed converts data to pdf so overall integral is 1
    #freq is the no. of data points in the bin
    #bins_int_dpts are the bins interval data points
    freq, bins_int_dpts = np.histogram(polarisation_fraction, nbr_bins ,normed=False) 
    
    #cumulative distribution function
    #cumfreq is cumulative frequency
    cumfreq = freq.cumsum() 
    #normalise to pdf
    cumprob = cumfreq / cumfreq[-1]
    
    #use linear interpolation of cdf to find new pixel values
    eqdata = np.interp(polarisation_fraction,bins_int_dpts[:-1],cumfreq)
    return eqdata, cumprob
    
def scaling(method_name, polarisation_fraction, nbr_bins, nbr_times):
    if method_name == 'hist':
        eqdata = polarisation_fraction
        for i in range(nbr_times):    
            eqdata, cumprob = histeq(eqdata, nbr_bins)
            if i == 0:
                eqdata1 = eqdata
    elif method_name == 'linear':
        eqdata = polarisation_fraction
        eqdata1 = eqdata
    elif method_name == 'squared':
        eqdata = polarisation_fraction**2
        eqdata1 = eqdata
    elif method_name == 'root_2':
        eqdata = np.sqrt(polarisation_fraction)
        eqdata1 = eqdata
    elif method_name == 'root_3':
        eqdata = (polarisation_fraction)**(1/3.)
        eqdata1 = eqdata
    elif method_name == 'root_4':
        eqdata = (polarisation_fraction)**(1/4.)
        eqdata1 = eqdata
    elif method_name == 'root_5':
        eqdata = (polarisation_fraction)**(1/5.)
        eqdata1 = eqdata
    elif method_name == 'root_6':
        eqdata = (polarisation_fraction)**(1/6.)
        eqdata1 = eqdata
    elif method_name == 'log':
        eqdata = np.log(polarisation_fraction)
        eqdata1 = eqdata
    else:
        raise TypeError("Wrong method_name")
    return eqdata1, eqdata  
    
aaa = 'K'
original_directory = os.getcwd()
    
if aaa == 'K':
    directory = "/Users/Administer/Desktop/Planck data"
    #directory = '/home/kelvin/Desktop/Python'
    os.chdir(directory)

input_filename = 'COM_CompMap_DustPol-commander_1024_R2.00.fits'
intensity_fits = 'COM_CompMap_ThermalDust-commander_2048_R2.00.fits'


#Checks if the fits file exists. If it doesn't, you choose the folder that contains the fits file.
while os.path.isfile(input_filename) == False:
    Tk.Tk().withdraw()
    directory = tkFileDialog.askdirectory()
    os.chdir(directory)

QUtable = fits.open(input_filename)
intensitytable = fits.open(intensity_fits)

I = hp.read_map(intensity_fits, field=0, nest=False)
Q = hp.read_map(input_filename, field=0, nest=False)
U = hp.read_map(input_filename, field=1, nest=False)


nside = 2**5
I = hp.ud_grade(I, nside)
Q = hp.ud_grade(Q, nside)
U = hp.ud_grade(U, nside)
polarisation_fraction = np.sqrt(Q**2 + U**2)/I

nbr_bins = 1000 #no. of bins
nbr_times = 2 #no of times of histogram equalising
for method_name in ['hist']:
    r = 70 #unit in mm
    scale_range = 10 #unit in mm
    threeD_model = False
    fourviews = False
    print 'Reformat:'
    print 'nside reduced to = %s' %(nside)
    print 'new no. of data points =', 12*nside**2
    if method_name == 'hist':
        print 'no. of bins =', nbr_bins
        
    eqdata1, eqdata = scaling(method_name, polarisation_fraction, nbr_bins, nbr_times)
    theta, phi = hp.pix2ang(nside, range(polarisation_fraction.size), nest=False)
    
    #Radius scaling
    r1 = eqdata - np.array(min(eqdata))
    r2 = scale_range * r1/max(r1) + np.array([r]*eqdata.size)
    
    x = -r2*np.sin(theta)*np.cos(phi) #flip probably due to differnt coord systems
    y = r2*np.sin(theta)*np.sin(phi)
    z = r2*np.cos(theta)

    ary = np.array([x,y,z,r2])
    ary = np.transpose(ary)
    os.chdir(original_directory)
    np.savetxt('DUST_bighist_nside%s_r%s_rscale%s.xyz' %(nside,r, scale_range), ary)  
    
plt.figure(2)
hp.mollview(map=Q, fig=2, title="Mollview - Q", nest=False, norm="hist" , sub=231, xsize=2000) #, min=-5e-7, max=5e7, xsize=2000)
hp.mollview(map=U, fig=2, title="Mollview - U", nest=False, norm="hist", sub=232, xsize=2000) #, min=-5e-7, max=5e7, xsize=2000)
hp.mollview(map=polarisation_fraction, fig=2, title="Mollview - Polarisation Fraction", nest=False, norm="hist", sub=233, xsize=2000) #, min=-5e-7, max=5e7, xsize=2000)
hp.mollview(map=Q, fig=2, title="Mollview - Q", nest=False, norm=" " , sub=234, xsize=2000) #, min=-5e-7, max=5e7, xsize=2000)
hp.mollview(map=U, fig=2, title="Mollview - U", nest=False, norm=" ", sub=235, xsize=2000) #, min=-5e-7, max=5e7, xsize=2000)
hp.mollview(map=polarisation_fraction, fig=2, title="Mollview - Polarisation Fraction", nest=False, norm=" ", sub=236, xsize=2000) #, min=-5e-7, max=5e7, xsize=2000)
plt.suptitle("Polarised Dust")
plt.show()














