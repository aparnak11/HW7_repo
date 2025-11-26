#!/usr/bin/env python3
import pyvista as pv
import numpy as np

""" image = pv.read("field.vti")
print(image.dimensions)  # should be (cols, rows, 1)
image.plot(scalars="T", cpos="xy") """

# ---- Load from C++ ----
data = np.loadtxt("pulse.csv", delimiter=",", skiprows=1)
t = data[:,0]          # time [s]
y = data[:,1]          # 0 = off, 1 = on

# ---- Build polyline for PyVista ----
points = np.column_stack((t, y, np.zeros_like(t)))
npts = len(points)
lines = np.hstack(([npts], np.arange(npts)))

plt = pv.Plotter()
curve = pv.PolyData(points)
curve.lines = lines

# ---- Render ----
plt.add_mesh(curve, color="white", line_width=3)
plt.show_axes()
plt.camera_position = "xy"
plt.show()

