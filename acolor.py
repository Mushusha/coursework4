import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy.interpolate import griddata
from matplotlib.colors import ListedColormap

x = np.loadtxt('build/coord_x.txt')
y = np.loadtxt('build/coord_y.txt')

z = np.loadtxt('build/disp_sum.txt')


xi = np.linspace(min(x), max(x), 1000)
yi = np.linspace(min(y), max(y), 1000)
zi = griddata((x, y), z, (xi[None, :], yi[:, None]), method='cubic')  # интерполяция

# Построение цветовой карты
plt.figure(figsize=(12, 6))
plt.pcolormesh(xi, yi, zi, cmap='rainbow', shading='auto')
plt.colorbar(label='')
plt.xlabel('X')
plt.ylabel('Y')
plt.show()
