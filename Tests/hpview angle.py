import healpy as hp
import astropy.io.fits as ap
import astropy.table as table
import numpy as np
import matplotlib.pyplot as plt
import pylab
from matplotlib.pyplot import Slider, Button

filename = 'commander_synch.fits'
data = hp.fitsfunc.read_map(filename)
b = table.Table.read(filename, hdu=1)

fig1 = plt.figure(1, figsize = (16,9))

hp.orthview(data, fig=1, sub=211, norm='hist', half_sky=True)

axcolor = 'lightgoldenrodyellow'
axlatitude = plt.axes([0.2, 0.35, 0.65, 0.03], axisbg=axcolor)
axlongitude = plt.axes([0.2, 0.4, 0.65, 0.03], axisbg=axcolor)
axangle = plt.axes([0.2, 0.45, 0.65, 0.03], axisbg=axcolor)

latitude_slider = Slider(axlatitude, 'Latitude', -90, 90, valinit=0)
longitude_slider = Slider(axlongitude, 'Longitude', -90, 90, valinit=0)
angle_slider = Slider(axangle, 'Angle', -90, 90, valinit=0)

def update(val):
    latitude = latitude_slider.val
    longitude = longitude_slider.val
    angle = angle_slider.val
    hp.orthview(data, fig=1, sub=211, half_sky = True, rot=(longitude, latitude, angle), norm='hist')

latitude_slider.on_changed(update)
longitude_slider.on_changed(update)
angle_slider.on_changed(update)

plt.show()