import matplotlib.pyplot as plt
import numpy as np

data = np.loadtxt("pulse.csv", delimiter=",", skiprows=1)
t = data[:,0]          # time [s]
y = data[:,1]          # 0 = off, 1 = on

plt.plot(t, y)
plt.show()