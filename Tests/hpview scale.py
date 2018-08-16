import healpy as hp
import astropy.io.fits as ap
import astropy.table as table
import numpy as np
import matplotlib.pyplot as plt
import pylab
from mayavi import mlab
from matplotlib.pyplot import Slider, Button

filename = 'commander_synch.fits'
data = hp.fitsfunc.read_map(filename)

fig1 = plt.figure(1, figsize = (16,9))

scale_max = max(data)
scale_min = 0

hp.orthview(data, fig=1, min=scale_min, max=scale_max, sub=221)#, norm='hist')
hp.orthview(data, fig=1, min=scale_min, max=scale_max, sub=222, norm='hist')

axcolor = 'lightgoldenrodyellow'
axmax = plt.axes([0.1, 0.4, 0.3, 0.03], axisbg=axcolor)
axmin = plt.axes([0.1, 0.45, 0.3, 0.03], axisbg=axcolor)

scale_max_slider = Slider(axmax, 'Max', scale_min, scale_max, valinit=scale_max)
scale_min_slider = Slider(axmin, 'Min', scale_min, scale_max, valinit=scale_min)

def update(val):
    scale_max = scale_max_slider.val
    scale_min = scale_min_slider.val
    hp.orthview(data, fig=1, min=scale_min, max=scale_max, sub=221)#, norm='hist')

scale_max_slider.on_changed(update)
scale_min_slider.on_changed(update)


plt.show()