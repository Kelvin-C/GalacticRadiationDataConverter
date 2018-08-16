import numpy as np
import matplotlib.pyplot as plt

x = np.linspace(-200, 200, 400)
y = np.linspace(-200, 200, 400)
z = np.linspace(-200, 200, 400)

r = np.sqrt(x**2 + y**2 + z**2)

maxheight = 65.
minradius = 40.
maxradius = maxheight + minradius
#r = np.linspace(minradius, minradius+maxheight, 400)
radii = [minradius, minradius + maxheight/2., maxradius]
frequency = np.pi/maxheight
R = np.cos(frequency * (r - radii[-1])) * (127) + 128;
G = np.cos(frequency * (r - radii[1])) * (127) + 128;
B = np.cos(frequency * (r - radii[0])) * (127) + 128;

print "R = cos(3.14/%s * (sqrt(x^2 + y^2 + z^2) - %s)) * (127) + 128" %(maxheight, radii[-1])
print "G = cos(3.14/%s * (sqrt(x^2 + y^2 + z^2) - %s)) * (127) + 128" %(maxheight, radii[1])
print "B = cos(3.14/%s * (sqrt(x^2 + y^2 + z^2) - %s)) * (127) + 128" %(maxheight, radii[0])

def plott():
    plt.plot(r,R, 'r')
    plt.plot(r,G, 'g')
    plt.plot(r,B, 'b')
    plt.axis([minradius, maxradius, 0, 255])
    plt.show()

#plott()