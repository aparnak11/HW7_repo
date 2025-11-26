import numpy as np

# -----------------------------
# Part 1: 1D FTCS heat equation
# -----------------------------

# Physical and numerical parameters
L = 0.01          # domain length [m]
nx = 100          # number of spatial points
nt = 15000        # number of time steps
dt = 1.0e-4       # time step [s]

kappa = 65.0      # thermal conductivity [W/(m·K)]
alpha = 1.7e-5    # thermal diffusivity [m^2/s]
hg = 1.5e4        # convective heat transfer coefficient [W/(m^2·K)]
Tg = 2500.0       # gas stagnation temperature [K]

T_init = 300.0    # initial wall temperature [K]

# Pulse parameters (periodic firing)
pulse_start = 0.0     # first pulse start time [s]
pulse_width = 0.15    # pulse duration [s]
pulse_period = 0.2    # pulse repetition period [s]

dx = L / (nx - 1)
x = np.linspace(0.0, L, nx)
t = np.arange(nt) * dt

r_stab = alpha * dt / dx**2
print("Stability coefficient r =", r_stab)

# Allocate temperature array
# We'll store the full space-time field as Ttime[nt, nx]
T = np.full(nx, T_init, dtype=float)
Ttime = np.zeros((nt, nx), dtype=float)

for n in range(nt):
    # --- store current time slice ---
    Ttime[n, :] = T

    # --- determine if thruster is firing ---
    current_time = n * dt
    tau = (current_time - pulse_start) % pulse_period
    thruster_on = (current_time >= pulse_start) and (tau < pulse_width)

    # --- compute ghost nodes for this time step ---
    # right side: always zero-flux, so Tghost_right = T_{nx-2}
    T_ghost_right = T[-2]

    # left side:
    if thruster_on:
        # convection BC: k dT/dx = h_g (Tg - T0)
        T0 = T[0]
        T_ghost_left = T[1] - 2.0 * dx * hg / kappa * (Tg - T0)
    else:
        # insulated: zero-flux so mirror
        T_ghost_left = T[1]

    # --- FTCS update ---
    T_new = np.empty_like(T)

    # interior nodes: i = 1..nx-1
    for i in range(1, nx-1):
        T_new[i] = T[i] + r_stab * (T[i+1] - 2.0*T[i] + T[i-1])

    # left boundary: i = 0, use ghost_left and neighbor 1
    i = 0
    T_new[i] = T[i] + r_stab * (T[1] - 2.0*T[i] + T_ghost_left)

    # right boundary: i = nx-1, use neighbor nx-2 and ghost_right
    i = nx - 1
    T_new[i] = T[i] + r_stab * (T_ghost_right - 2.0*T[i] + T[nx-2])

    # advance in time
    T = T_new

# -----------------------------------------
# Helper: write 2D scalar field to VTI file
# -----------------------------------------

def write_vti_2d_scalar(filename, data, dx, dy, origin=(0.0, 0.0, 0.0),
                        spacing_z=1.0, scalar_name="Temperature"):
    """
    Write a 2D numpy array 'data' (ny, nx) as a VTK ImageData (.vti).

    - x-direction: columns (fastest axis)
    - y-direction: rows (next axis)
    - z-direction: single layer

    'dx' = spacing in x, 'dy' = spacing in y (here: dt).
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

        # VTK ImageData expects i (x) fastest, then j (y), then k (z)
        # Here nz = 1
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

# Write x–t contour as ImageData:
# x -> horizontal, t -> vertical
write_vti_2d_scalar(
    "part1_heat_xt.vti",
    data=Ttime,     # shape (nt, nx)
    dx=dx,          # spatial spacing
    dy=dt,          # "time spacing" in seconds
    origin=(0.0, 0.0, 0.0),
    spacing_z=1.0,
    scalar_name="Temperature"
)
print("Wrote VTI file: part1_heat_xt.vti")
