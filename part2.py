import numpy as np
from scipy.sparse import lil_matrix, identity
from scipy.sparse.linalg import spsolve

# -------------------------------
# Part 2: Crank-Nicolson diffusion
# Open cylinder: inlet + outlet window
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
            # (1/r) * d/dr term
            coef_rp =  1.0 / (2.0 * rj * dr)
            coef_rm = -1.0 / (2.0 * rj * dr)
            # second-derivative radial term
            coef_r2 = 1.0 / dr**2

            # self
            L[p, p] += -2.0 * coef_r2
            # j+1
            L[p, idx(i, j+1)] += coef_rp + coef_r2
            # j-1
            L[p, idx(i, j-1)] += coef_rm + coef_r2
        else:
            # temporary one-sided approx; Neumann imposed later by mirroring
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
            # temporary; Neumann/Dirichlet handled separately
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

# Helper to apply Dirichlet BCs
def apply_dirichlet(A, B, rhs, value_mask, values):
    """
    value_mask : boolean array of length N, where True means Dirichlet node.
    values     : array of Dirichlet values (same length N, arbitrary where mask is False)
    """
    A = A.tolil()
    B = B.tolil()
    for p in np.where(value_mask)[0]:
        A.rows[p] = [p]
        A.data[p] = [1.0]

        B.rows[p] = [p]
        B.data[p] = [0.0]

        rhs[p] = values[p]     # enforce n^{k+1}_p = value
    return A.tocsr(), B.tocsr(), rhs

# Dirichlet mask and values for OPEN cylinder
dirichlet_mask = np.zeros(N, dtype=bool)
dirichlet_values = np.zeros(N, dtype=float)

# Inlet: z = 0, 0 <= r <= 0.04 => r index j where r[j] <= 0.04
for j in range(nr):
    if r[j] <= 0.04 + 1e-12:
        i = 0
        p = idx(i, j)
        dirichlet_mask[p] = True
        dirichlet_values[p] = 10.0   # concentration 10

# Outlet window: z = 1, 0.1 <= r <= 0.2
for j in range(nr):
    if 0.1 - 1e-12 <= r[j] <= 0.2 + 1e-12:
        i = nz - 1
        p = idx(i, j)
        dirichlet_mask[p] = True
        dirichlet_values[p] = 0.0    # concentration 0

# Copy A,B; these will be reused each step
A_base = A.copy()
B_base = B.copy()

# Time integration with Crank–Nicolson
snapshot_steps = [100, 400, 1000, 2000, 5000]
snapshots = {}

for k in range(1, nsteps+1):
    # RHS = B * n^k
    rhs = B_base.dot(n_vec)

    # Apply Dirichlet for this step
    A_step, B_step, rhs = apply_dirichlet(A_base, B_base, rhs,
                                          dirichlet_mask, dirichlet_values)

    # Solve A n^{k+1} = rhs
    n_new = spsolve(A_step, rhs)

    # Enforce Neumann BCs by mirroring interior nodes
    n_grid = n_new.reshape((nr, nz))  # shape (nr, nz), [j, i] = [r, z]

    # r = 0 boundary: mirror from j=1
    n_grid[0, :] = n_grid[1, :]
    # r = 0.5 boundary: mirror from j=nr-2
    n_grid[-1, :] = n_grid[-2, :]

    # z = 0 and z = 1: where not Dirichlet, treat as Neumann
    for j in range(nr):
        # z=0
        p0 = idx(0, j)
        if not dirichlet_mask[p0]:
            n_grid[j, 0] = n_grid[j, 1]
        # z=1
        p1 = idx(nz-1, j)
        if not dirichlet_mask[p1]:
            n_grid[j, -1] = n_grid[j, -2]

    n_vec = n_grid.reshape(N)

    if k in snapshot_steps:
        snapshots[k] = n_grid.copy()   # shape (nr, nz)

print("Finished Crank-Nicolson diffusion (open cylinder).")

# -----------------------------------------
# VTI writers
# -----------------------------------------

def write_vti_2d_scalar(filename, data, dx, dy, origin=(0.0, 0.0, 0.0),
                        spacing_z=1.0, scalar_name="n"):
    """
    Generic 2D numpy array 'data' (ny, nx) -> VTK ImageData (.vti).
    """
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
    """
    data : shape (nr, nz) with:
        axis 0 = r index (vertical)
        axis 1 = z index (horizontal)
    """
    nr, nz = data.shape
    # rows = r, cols = z
    write_vti_2d_scalar(
        filename,
        data=data,    # (ny, nx) = (nr, nz)
        dx=dz,        # spacing in z (horizontal axis)
        dy=dr,        # spacing in r (vertical axis)
        origin=(0.0, 0.0, 0.0),
        spacing_z=1.0,
        scalar_name=scalar_name
    )

# Write snapshots to VTI
for step, field in snapshots.items():
    fname = f"part2_open_t{step}.vti"
    write_vti_2d_scalar_zr(fname, field, dz, dr, scalar_name="Density")
    print("Wrote", fname)
