import numpy as np
import matplotlib.pyplot as plt

x = np.arange(20)
y1 = np.log(x)
y2 = x**(1/2.)
y3 = x**(1/3.)
y4 = x**(1/4.)
y5 = x**(1/5.)
y6 = x**(1/6.)

colors = ['r', 'y','g','b','c','m']
ys = [y1,y2,y3,y4,y5,y6]
for i in range(len(ys)):
    plt.plot(x, ys[i], colors[i])
plt.legend()
plt.show()

