import pyvista as pv

grid = pv.read("part1_heat_xt.vti")
p = pv.Plotter()
p.add_mesh(grid, cmap="inferno")
p.show()
