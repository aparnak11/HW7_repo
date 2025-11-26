#!/usr/bin/env python3
import pyvista as pv
image = pv.read("field.vti")
image.plot(scalars="T", cpos="xy")