import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

n_total_pixels = []
n_pixels_order = []
for i in range(3):
    n = i
    total_no_of_lines = 2**(n+2) - 1
    total_no_of_middle_lines = 2**(n+1) +1
    no_of_pixels_in_equator = 2**(n+2)
    
    no_of_remaining_lines_in_each_hemisphere = 2**(n+1) - 2**n -1
    lst1 = []
    for i in range(no_of_remaining_lines_in_each_hemisphere):
        temp = 4 +i*4
        lst1 += [temp]
    lst2 = lst1[::-1]
    pixels_order_for_each_n = lst1 + [no_of_pixels_in_equator]*total_no_of_middle_lines + lst2
    n_pixels_order += [pixels_order_for_each_n]
    total_no_of_pixels_in_sphere = sum(lst1)*2 + total_no_of_middle_lines*no_of_pixels_in_equator
    n_total_pixels += [total_no_of_pixels_in_sphere]
    
    
n_pixels = np.array(n_total_pixels)
print "Total no. of pixles for different n:", n_total_pixels
print "-"*60
print "Check if obeys 12*4**n rule:", np.log((n_pixels/12))/np.log(4)
print "-"*60
print "No. of pixels per line for different n:", n_pixels_order

total_pixel = n_total_pixels[-1]
pixel_list = n_pixels_order[-1]
total_lines = len(pixel_list)
ary = np.zeros((3, total_pixel))

r = 1
ary[0] = r

del_theta = np.pi/total_lines
j = 0
del_phi = 0
for i in range(total_lines):
    if j%2 == 0:
        initialphi = -0.5*del_phi
    del_phi = 2*np.pi/pixel_list[i]
    theta = 0.5*del_theta + i*del_theta
    for n in range(pixel_list[i]):
        phi = 2*np.pi - n*del_phi
        ary[1][j], ary[-1][j] = theta, phi
        j += 1

#ary = np.transpose(ary)

x = ary[0]*np.sin(ary[1])*np.cos(ary[2])
y = ary[0]*np.sin(ary[1])*np.sin(ary[2])
z = ary[0]*np.cos(ary[1])


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
Axes3D.scatter(ax, x,y,zs=z)

plt.show()