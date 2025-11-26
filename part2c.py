import numpy as np
from scipy.sparse import lil_matrix, identity
from scipy.sparse.linalg import spsolve

# -------------------------------
# Part 2c: Crank-Nicolson diffusion
# Closed cylinder: inlet only
# -------------------------------

D = 0.1            # diffusion coefficient
dt = 0.01          # time step
nsteps = 5000

# grid: z in [0,1], r in [0,0.5]
nz = 51            # number of points in z
nr = 26            # number of points in r
dz = 1.0 / (nz - 1)
dr = 0.5 / (nr - 1)

z = np.linspace(0.0, 1.0, nz)
r = np.linspace(0.0, 0.5, nr)

def idx(i, j):
    """Map (i,j) -> 1D index. i: z index [0..nz-1], j: r index [0..nr-1]."""
    return j * nz + i

N = nz * nr

# Build Laplacian matrix L in axisymmetric coordinates
L = lil_matrix((N, N))

for j in range(nr):
    for i in range(nz):
        p = idx(i, j)
        rj = max(r[j], 1e-8)  # avoid division by zero at r=0

        # radial part
        if 0 < j < nr - 1:
            coef_rp =  1.0 / (2.0 * rj * dr)
            coef_rm = -1.0 / (2.0 * rj * dr)
            coef_r2 = 1.0 / dr**2

            L[p, p] += -2.0 * coef_r2
            L[p, idx(i, j+1)] += coef_rp + coef_r2
            L[p, idx(i, j-1)] += coef_rm + coef_r2
        else:
            L[p, p] += -2.0 / dr**2
            if j == 0 and nr > 1:
                L[p, idx(i, j+1)] += 2.0 / dr**2
            elif j == nr - 1 and nr > 1:
                L[p, idx(i, j-1)] += 2.0 / dr**2

        # axial part
        if 0 < i < nz - 1:
            L[p, p] += -2.0 / dz**2
            L[p, idx(i+1, j)] += 1.0 / dz**2
            L[p, idx(i-1, j)] += 1.0 / dz**2
        else:
            L[p, p] += -2.0 / dz**2
            if i == 0 and nz > 1:
                L[p, idx(i+1, j)] += 2.0 / dz**2
            elif i == nz - 1 and nz > 1:
                L[p, idx(i-1, j)] += 2.0 / dz**2

L = L.tocsr()

# Identity
I = identity(N, format="csr")
A = (I - 0.5 * D * dt * L).tocsr()
B = (I + 0.5 * D * dt * L).tocsr()

# Initial condition: n = 0 everywhere
n_vec = np.zeros(N, dtype=float)

def apply_dirichlet(A, B, rhs, value_mask, values):
    """
    Apply Dirichlet conditions by overwriting rows in A and B.
    """
    A = A.tolil()
    B = B.tolil()
    for p in np.where(value_mask)[0]:
        A.rows[p] = [p]
        A.data[p] = [1.0]

        B.rows[p] = [p]
        B.data[p] = [0.0]

        rhs[p] = values[p]
    return A.tocsr(), B.tocsr(), rhs

# Dirichlet mask and values for CLOSED cylinder (only inlet)
dirichlet_mask = np.zeros(N, dtype=bool)
dirichlet_values = np.zeros(N, dtype=float)

# Inlet: z = 0, 0 <= r <= 0.04
for j in range(nr):
    if r[j] <= 0.04 + 1e-12:
        i = 0
        p = idx(i, j)
        dirichlet_mask[p] = True
        dirichlet_values[p] = 10.0

A_base = A.copy()
B_base = B.copy()

snapshot_steps = [100, 400, 1000, 2000, 5000]
snapshots = {}

for k in range(1, nsteps+1):
    rhs = B_base.dot(n_vec)
    A_step, B_step, rhs = apply_dirichlet(A_base, B_base, rhs,
                                          dirichlet_mask, dirichlet_values)

    n_new = spsolve(A_step, rhs)

    # Neumann BCs by mirroring
    n_grid = n_new.reshape((nr, nz))

    # r = 0, r = 0.5
    n_grid[0, :] = n_grid[1, :]
    n_grid[-1, :] = n_grid[-2, :]

    # z = 0,1: wherever not Dirichlet, mirror interior
    for j in range(nr):
        p0 = idx(0, j)
        if not dirichlet_mask[p0]:
            n_grid[j, 0] = n_grid[j, 1]

        p1 = idx(nz-1, j)
        if not dirichlet_mask[p1]:
            n_grid[j, -1] = n_grid[j, -2]

    n_vec = n_grid.reshape(N)

    if k in snapshot_steps:
        snapshots[k] = n_grid.copy()

print("Finished Crank-Nicolson diffusion (closed cylinder).")

# -----------------------------------------
# VTI writers
# -----------------------------------------

def write_vti_2d_scalar(filename, data, dx, dy, origin=(0.0, 0.0, 0.0),
                        spacing_z=1.0, scalar_name="n"):
    ny, nx = data.shape
    nz = 1
    ox, oy, oz = origin

    with open(filename, "w") as f:
        f.write('<?xml version="1.0"?>\n')
        f.write('<VTKFile type="ImageData" version="0.1" byte_order="LittleEndian">\n')
        f.write(
            f'  <ImageData WholeExtent="0 {nx-1} 0 {ny-1} 0 {nz-1}" '
            f'Origin="{ox} {oy} {oz}" '
            f'Spacing="{dx} {dy} {spacing_z}">\n'
        )
        f.write(
            f'    <Piece Extent="0 {nx-1} 0 {ny-1} 0 {nz-1}">\n'
            f'      <PointData Scalars="{scalar_name}">\n'
            f'        <DataArray type="Float32" Name="{scalar_name}" format="ascii">\n'
        )

        for j in range(ny):
            line_vals = []
            for i in range(nx):
                line_vals.append(str(float(data[j, i])))
            f.write("          " + " ".join(line_vals) + "\n")

        f.write(
            "        </DataArray>\n"
            "      </PointData>\n"
            "      <CellData/>\n"
            "    </Piece>\n"
            "  </ImageData>\n"
            "</VTKFile>\n"
        )

def write_vti_2d_scalar_zr(filename, data, dz, dr, scalar_name="n"):
    nr, nz = data.shape
    write_vti_2d_scalar(
        filename,
        data=data,    # (ny, nx) = (nr, nz)
        dx=dz,
        dy=dr,
        origin=(0.0, 0.0, 0.0),
        spacing_z=1.0,
        scalar_name=scalar_name
    )

for step, field in snapshots.items():
    fname = f"part2_closed_t{step}.vti"
    write_vti_2d_scalar_zr(fname, field, dz, dr, scalar_name="Density")
    print("Wrote", fname)
