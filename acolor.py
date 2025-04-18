import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from scipy.interpolate import griddata
from matplotlib.colors import ListedColormap

x = np.loadtxt('build/coord_x.txt')
y = np.loadtxt('build/coord_y.txt')

z = np.loadtxt('build/disp_y.txt')


xi = np.linspace(min(x), max(x), 100)
yi = np.linspace(min(y), max(y), 100)
zi = griddata((x, y), z, (xi[None, :], yi[:, None]), method='cubic')  # интерполяция

# Построение цветовой карты
plt.figure(figsize=(8, 6))
plt.pcolormesh(xi, yi, zi, cmap='coolwarm', shading='auto')
plt.colorbar(label='')
plt.xlabel('X')
plt.ylabel('Y')
plt.show()
