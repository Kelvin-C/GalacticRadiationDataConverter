from astropy.io import fits
import astropy.table as table
import healpy as hp
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

plt.close('all')

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
    
def pivoting_ball_radius(r, theta, phi, scale_range):
    x = r*np.sin(theta)*np.cos(phi)
    y = r*np.sin(theta)*np.sin(phi)
    z = r*np.cos(theta)
    ary = np.array([x,y,z])
    ary = np.transpose(ary)
    base = np.linalg.norm(ary[0]-ary[1])
    radius = np.sqrt(base**2 + scale_range**2)
    angle = 180. / np.pi * np.arctan(scale_range/base)
    print 'base = ', base
    if radius > 2*base:
        print '-'*40
        print 'Error: radius needs to be less than 2 * base'
        print '2*base =', 2*base
        print 'maximum scale range =', np.sqrt((2*base)**2 - base**2)
        print '-'*40
    return radius, angle
    
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


input_filename = "COM_CompMap_Synchrotron-commander_0256_R2.00.fits"
print '-'*40
print 'Input file details:'
intmap = hp.read_map(input_filename)
print '-'*40

nside = 2**8 #default 2**8
nbr_bins = 1000 #no. of bins
nbr_times = 2 #no of times of histogram equalising
for method_name in ['hist','linear','squared','square root','log']:
    #method_name = 'square '
    ball_pivoting = False
    r = 30 #unit in mm
    scale_range = 60 #unit in mm
    threeD_model = False
    fourviews = False
    print 'Reformat:'
    print 'nside reduced to = %s' %(nside)
    print 'new no. of data points =', 12*nside**2
    if method_name == 'hsit':
        print 'no. of bins =', nbr_bins
    
    intmap = hp.ud_grade(intmap, nside)
    eqdata1, eqdata = scaling(method_name, intmap, nbr_bins, nbr_times)
    theta, phi = hp.pix2ang(nside, range(intmap.size), nest=False)
    
    #Radius scaling
    r1 = eqdata - np.array(min(eqdata))
    r2 = scale_range * r1/max(r1) + np.array([r]*eqdata.size)
    
    rcolourindex = [[],[],[],[],[]]
    maxr = max(r2)
    minr = min(r2)
    diffr = maxr-minr
    for i in range(len(r2)):
        if r2[i] <= (r + 0.2*diffr):
            rcolourindex[0] += [i]
        elif r2[i] <= (r + 0.4*diffr): 
            rcolourindex[1] += [i]
        elif r2[i] <= (r + 0.6*diffr): 
            rcolourindex[2] += [i]
        elif r2[i] <= (r + 0.8*diffr): 
            rcolourindex[3] += [i]
        elif r2[i] <= (r + 1.*diffr): 
            rcolourindex[4] += [i]
    
    x = -r2*np.sin(theta)*np.cos(phi) #flip probably due to differnt coord systems
    y = r2*np.sin(theta)*np.sin(phi)
    z = r2*np.cos(theta)
    
    coords = [[],[],[],[],[]]
    for i in range(len(rcolourindex)):
        for j in rcolourindex[i]:
            coords[i] += [[x[j], y[j], z[j]]]
    
    ary = np.array([x,y,z])
    ary = np.transpose(ary)
    np.savetxt('nside%s_nbr%s__nbrtimes%s_r%s_sr%s_method%s.xyz' %(nside, nbr_bins, nbr_times, r, scale_range, method_name), ary)  
    
    """
    PHI, THETA = np.meshgrid(phi, theta)
    grid_pix = hp.ang2pix(nside, THETA, PHI, nest=True)
    grid_map = eqdata[grid_pix]
    x2 = r*np.sin(THETA)*np.cos(PHI)
    y2 = -r*np.sin(THETA)*np.sin(PHI)
    z2 = r*np.cos(THETA) 
    mlab.figure(1, bgcolor=(1, 1, 1), fgcolor=(0, 0, 0), size=(50, 50))
    mlab.clf()
    mlab.mesh(x2, y2, z2, scalars=grid_map, colormap="jet")#, vmin=vmin, vmax=vmax) hel
    mlab.savefig("test7.jpg")
    """
    
    
    if method_name == 'hist' and ball_pivoting == True:
        #for triangular sufrace reconstruction only
        radius, angle = pivoting_ball_radius(r, theta, phi, scale_range)
        print 'pivoting_ball_radius =', radius
        
        piv2_lst =[]
        piv2_lst2 =[]
        
        for i in range(len(ary)-1): #next point distance
            temp = np.linalg.norm(ary[i+1]-ary[i])
            piv2_lst += [temp]
        #piv2_lst.sort()
        print 'pivoting_ball_radius2 lies between: ',  2*min(piv2_lst) #max(piv2_lst)
        
        for i in range(len(ary)-2): #next point distance
            temp = np.linalg.norm(ary[i+2]-ary[i])
            piv2_lst2 += [temp]
        #piv2_lst2.sort()
        print 'pivoting_ball_radius3 lies between: ', min(piv2_lst2), max(piv2_lst)
        
        print 'angle =', angle
    
    if method_name == 'hist':
        #plt hist graphs (apply to histeq technique only)
        #plt histogram of original intensity data
        plt.subplot(4, 1, 1)
        plt.hist(intmap, nbr_bins, normed=True)#, range = (0,0.1e9))
        plt.title("Histogram of original data")
        plt.xlabel("Intensity")
        plt.ylabel("Frequency")
        #plt.axis([0, 0.01e8, 0, 0.000005])
        
        #convert histogram into cumulative frequency distribution
        plt.subplot(4, 1, 2)
        plt.hist(intmap, nbr_bins, normed=False, cumulative=True)
        plt.title("Cumulative frequency")
        plt.xlabel("Intensity")
        plt.ylabel("Cumulative frequency")
        
        #Histogram equalised data
        plt.subplot(4, 1, 3)
        plt.hist(eqdata1, nbr_bins, normed=True)#, range = (0,0.1e9))
        plt.title("Histogram of 1st equalised data")
        plt.xlabel("Cumulative frequency")
        plt.ylabel("Frequency")
        #plt.axis([0, 0.01e8, 0, 0.000005])
        
        #plt after-histogram #plotline
        plt.subplot(4, 1, 4)
        plt.hist(eqdata, nbr_bins, normed=True)#, range = (0,0.1e9))
        plt.title("Histogram of nth equalised data")
        plt.xlabel("Cumulative frequency")
        plt.ylabel("Frequency")
        plt.tight_layout()
        #plt.savefig("defhisteq_explain2.jpg", dpi=900)
    
        
    if threeD_model == True:
        #plt 3D model
        colours = ['blue', 'green' , 'yellow', 'orange', 'red']
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        for i in range(len(coords)):
            coords[i] = np.transpose(coords[i])
            if len(coords[i]) == 0:
                pass
            else:
                Axes3D.scatter(ax, coords[i][0], coords[i][1], zs=coords[i][2], c=colours[i])
        plt.show()
        
    if fourviews == True:
    #plt all views
        plt.figure(3)
        plt.figure(3).suptitle("Synchotron View", fontsize="x-large")
        
        hp.mollview(map=eqdata, fig=3, title="Mollview", nest=False, norm=" ", sub=221, xsize=2000) #, min=-5e-7, max=5e7, xsize=2000)
        hp.orthview(map=eqdata, fig=3, rot=None, coord=None, unit='', xsize=800, half_sky=False, title='Orthographic view', nest=False, min=None, max=None, flip='astro', remove_dip=False, remove_mono=False, gal_cut=0, format='%g', format2='%g', cbar=True, cmap=None, notext=False, norm=" ", hold=False, margins=None, sub=222, return_projected_map=False)
        hp.gnomview(map=eqdata, fig=3, rot=None, coord=None, unit='', xsize=200, ysize=None, reso=1.5, title='Gnomonic view', nest=False, remove_dip=False, remove_mono=False, gal_cut=0, min=None, max=None, flip='astro', format='%.3g', cbar=True, cmap=None, norm=" ", hold=False, sub=223, margins=None, notext=False, return_projected_map=False)
        hp.cartview(map=eqdata, fig=3, rot=None, zat=None, coord=None, unit='', xsize=800, ysize=None, lonra=None, latra=None, title='Cartesian view', nest=False, remove_dip=False, remove_mono=False, gal_cut=0, min=None, max=None, flip='astro', format='%.3g', cbar=True, cmap=None, norm=" ", aspect=None, hold=False, sub=224, margins=None, notext=False, return_projected_map=False)
    
    plt.show()


