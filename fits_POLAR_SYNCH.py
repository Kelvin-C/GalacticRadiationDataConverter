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


aaa = 'a'
original_directory = os.getcwd()
    
if aaa == 'K':
    directory = '/home/kelvin/Desktop/Python'
    os.chdir(directory)


#plt.close('all')
input_filename = 'COM_CompMap_SynchrotronPol-commander_0256_R2.00.fits'
intensity_fits = 'COM_CompMap_Synchrotron-commander_0256_R2.00.fits'


#Checks if the fits file exists. If it doesn't, you choose the folder that contains the fits file.
while os.path.isfile(input_filename) == False:
    Tk.Tk().withdraw()
    directory = tkFileDialog.askdirectory()
    os.chdir(directory)

QUtable = fits.open(input_filename)
intensitytable = fits.open(intensity_fits)

I1 = hp.read_map(intensity_fits, field=1, nest=False)
Q = hp.read_map(input_filename, field=0, nest=False)
U = hp.read_map(input_filename, field=1, nest=False)

nside = 2**3
I1 = hp.ud_grade(I1, nside)
Q = hp.ud_grade(Q, nside)
U = hp.ud_grade(U, nside)
I = np.sqrt(Q**2 + U**2) #Assume 100% polarised

polarisation_fraction = np.sqrt(Q**2 + U**2)/I1

E0x = np.sqrt((I + Q)/2.)
E0y = np.sqrt((I - Q)/2.)
ratio = U/(2.*E0x*E0y)
for i in range(len(ratio)):
    if ratio[i] > 1:
        ratio[i] = 1
    elif ratio[i] < -1:
        ratio[i] = -1
phasediff = np.arccos(ratio) #between E0x and E0y
theta_E = np.arctan(E0y/E0x*np.cos(phasediff)) #Polarisation angle

#QUtable.info()
#intensitytable.info()

#QUtable = table.Table.read(input_filename)
#intensitytable = table.Table.read(intensity_fits, hdu=1)
#print QUtable
#print "-"*100
#print intensitytable

#Q = intensitytable[0][:].data
#U = intensitytable[1][:].data

#Q = QUtable[0][:].data
#U = QUtable[1][:].data


#NOW add hist data
def histeq(intmap, nbr_bins):
    #compute histogram of data; normed converts data to pdf so overall integral is 1
    #freq is the no. of data points in the bin
    #bins_int_dpts are the bins interval data points
    freq, bins_int_dpts = np.histogram(intmap, nbr_bins ,normed=False) 
    
    #cumulative distribution function
    #cumfreq is cumulative frequency
    cumfreq = freq.cumsum() 
    #normalise to pdf
    cumprob = cumfreq / cumfreq[-1]
    
    #use linear interpolation of cdf to find new pixel values
    eqdata = np.interp(intmap,bins_int_dpts[:-1],cumfreq)
    return eqdata, cumprob
    
def scaling(method_name, intmap, nbr_bins, nbr_times):
    if method_name == 'hist':
        eqdata = intmap
        for i in range(nbr_times):    
            eqdata, cumprob = histeq(eqdata, nbr_bins)
            if i == 0:
                eqdata1 = eqdata
        
    elif method_name == 'linear':
        eqdata = intmap
        eqdata1 = eqdata
        
    elif method_name == 'squared':
        eqdata = intmap**2
        eqdata1 = eqdata
    elif method_name == 'square root':
        eqdata = np.sqrt(intmap)
        eqdata1 = eqdata
    elif method_name == 'log':
        eqdata = np.log(intmap)
        eqdata1 = eqdata
    return eqdata1, eqdata  
    
nbr_bins = 1000 #no. of bins
nbr_times = 2 #no of times of histogram equalising
method_name = 'hist'
r = 40 #unit in mm
scale_range = 6 #unit in mm #60

eqdata = I
for i in range(nbr_times):    
    eqdata, cumprob = histeq(eqdata, nbr_bins)
    if i == 0:
        eqdata1 = eqdata
        
#Radius scaling
r1 = eqdata - np.array(min(eqdata))
r2 = scale_range * r1/max(r1) + np.array([r]*eqdata.size)

theta, phi = hp.pix2ang(nside, range(I.size), nest=False)
x_eq = -r2*np.sin(theta)*np.cos(phi) #flip probably due to differnt coord systems
y_eq = r2*np.sin(theta)*np.sin(phi)
z_eq = r2*np.cos(theta)
#end of hist data


#Dist between 2 neighbouring points #can be simplified later
x = r*np.sin(theta)*np.cos(phi)
y = r*np.sin(theta)*np.sin(phi)
z = r*np.cos(theta)
ary = np.array([x,y,z])
ary = np.transpose(ary)
dist = np.linalg.norm(ary[0]-ary[1])

"""
#Use eq x,y,z to create ridges
ary = np.array([x_eq,y_eq,z_eq])
ary = np.transpose(ary)
x = x_eq
y = y_eq
z = z_eq
"""

"""
plt.figure(1)
hp.mollview(map=Q, fig=1, title="Mollview - Q", nest=False, norm="hist" , sub=221, xsize=2000) #, min=-5e-7, max=5e7, xsize=2000)
hp.mollview(map=U, fig=1, title="Mollview - U", nest=False, norm="hist", sub=222, xsize=2000) #, min=-5e-7, max=5e7, xsize=2000)
hp.mollview(map=polarisation_fraction, fig=1, title="Mollview", nest=False, norm="hist", sub=223, xsize=2000) #, min=-5e-7, max=5e7, xsize=2000)
"""

"""
plt.figure(2)
#travel with time
t = np.arange(1,10,0.1)
for i in np.arange(5):
    plt.plot(E0x[i]*np.cos(-t),E0y[i]*np.cos(-t+theta[i]))

plt.figure(3)
plt.plot(U[:10], Q[:10])
"""


theta_hat = np.array([np.cos(theta)*np.cos(phi),np.cos(theta)*np.sin(phi), -np.sin(theta)])
phi_hat = np.array([-np.sin(phi), np.cos(phi), np.zeros(phi.size)])
theta_hat = np.transpose(theta_hat) #array([[x1,y1,z1],[x2,y2,z2],...])
phi_hat = np.transpose(phi_hat)
theta_E = np.array([[temp] for temp in theta_E])
E0x_to_xyz = -phi_hat * np.cos(theta_E)
E0y_to_xyz = -theta_hat*np.sin(theta_E)
TotalE_in_xyz = E0x_to_xyz + E0y_to_xyz

n = list(np.zeros(len(ary)))
for i in range(len(ary)):
    n[i] = np.cross(ary[i], TotalE_in_xyz[i])
n = np.array(n)
n_hat = n/np.linalg.norm(n)

r_scale = 0.97 #Reduces radius of points between original points
dist_factor = 0.5
dist_factor2 = 5
no_of_forward = 6
no_of_backward = 6
forward = list(np.zeros(no_of_forward))
backward = list(np.zeros(no_of_backward))
    
sidef = list(np.zeros(2*no_of_forward))
sideb = list(np.zeros(2*no_of_backward))

for i in range(no_of_forward):
    forward[i] = ary + dist*dist_factor/no_of_forward*(i+1)*TotalE_in_xyz
    f1 = (forward[i] + n_hat*dist_factor2*dist)*r_scale
    f2 = (forward[i] - n_hat*dist_factor2*dist)*r_scale
    sidef[2*i] = f1
    sidef[2*i+1] = f2
    
for i in range(no_of_backward):
    backward[i] = ary - dist*dist_factor/no_of_backward*(i+1)*TotalE_in_xyz
    b1 = (backward[i] + n_hat*dist_factor2*dist)*r_scale
    b2 = (backward[i] - n_hat*dist_factor2*dist)*r_scale
    sideb[2*i] = b1
    sideb[2*i+1] = b2   
    
o1 = (ary + n_hat*dist_factor*dist)*r_scale
o2 = (ary - n_hat*dist_factor*dist)*r_scale

TotalE_in_xyz = np.transpose(TotalE_in_xyz) #array([lst x, lst y, lst z])
forward = np.array(forward)
forward = forward.reshape(no_of_forward*len(ary),3) 
forward = np.array([i/np.linalg.norm(i)*r for i in forward]) #plot points on sphere surface
forward = np.transpose(forward)
backward = np.array(backward)
backward = backward.reshape(no_of_backward*len(ary),3)
backward = np.array([i/np.linalg.norm(i)*r for i in backward]) #plot points on sphere surface
backward = np.transpose(backward)
sidef = np.array(sidef)
sidef = sidef.reshape(no_of_forward*2*len(ary),3)
sidef = np.array([i/np.linalg.norm(i)*r*r_scale for i in sidef]) #plot points on sphere surface
sidef = np.transpose(sidef)
sideb = np.array(sideb)
sideb = sideb.reshape(no_of_backward*2*len(ary),3)
sideb = np.array([i/np.linalg.norm(i)*r*r_scale for i in sideb]) #plot points on sphere surface
sideb = np.transpose(sideb)
o1 = np.array([i/np.linalg.norm(i)*r*r_scale for i in o1]) #plot points on sphere surface
o2 = np.array([i/np.linalg.norm(i)*r*r_scale for i in o2]) #plot points on sphere surface
o1 = np.transpose(o1)
o2 = np.transpose(o2)

#u,v,w are components of arrow vectors
u = -1* TotalE_in_xyz[0] #flip probably due to differnt coord systems
v = TotalE_in_xyz[1]
w = TotalE_in_xyz[2]

length = 4

all_xs = -1*np.array(list(x) + list(forward[0]) + list(backward[0]) + list(sidef[0]) + list(sideb[0]) + list(o1[0]) + list(o2[0]))
all_ys = list(y) + list(forward[1]) + list(backward[1]) + list(sidef[1]) + list(sideb[1]) + list(o1[1]) + list(o2[1])
all_zs = list(z) + list(forward[2]) + list(backward[2]) + list(sidef[2]) + list(sideb[2]) + list(o1[2]) + list(o2[2])
all_data = np.array([all_xs, all_ys, all_zs])
all_data = np.transpose(all_data)

os.chdir(original_directory)
np.savetxt('SynchPOLAR_nb%s_nside%s_rscale%s_df1_%s_df2_%s_r%s_sr%s_method%s.xyz' %(((no_of_forward+no_of_backward)*3 +2),nside, r_scale, dist_factor,dist_factor2, r, scale_range, method_name), all_data)  

#plot all data points, neighbour points and vector field
fig = plt.figure(5)
ax = fig.add_subplot(111, projection='3d')
#ax.plot_wireframe(x,y,z)
"""
ax.scatter(-1*x, y, z)
ax.scatter(-1*forward[0], forward[1], forward[2],c='r')
ax.scatter(-1*sidef[0], sidef[1], sidef[2], c='m')
ax.scatter(-1*backward[0], backward[1], backward[2],c='g')
ax.scatter(-1*sideb[0], sideb[1], sideb[2],c= 'c')
"""

#Quiver 
ax = fig.gca(projection='3d')
ax.quiver(-1*x, y, z, u, v, w,                 # data
          length=length,                    # arrow length
          color='Tomato'                    # arrow colour
          )            
ax.set_title('3D Vector Field - Synchrotron')             # title
#ax.view_init(elev=18, azim=30)              # camera elevation and angle
ax.dist= 7    

plt.show()
