#!/usr/bin/env python3
import pyvista as pv
image = pv.read("field.vti")
print(image.dimensions)  # should be (cols, rows, 1)
image.plot(scalars="T", cpos="xy")
