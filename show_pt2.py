import glob
import pyvista as pv

# ----------------------------------------------------------
# Choose which set to view and which snapshot index
# ----------------------------------------------------------
# which_set = "open"   # for part2_open_tXXXX.vti
which_set = "open"
# which_set = "closed" # for part2_closed_tXXXX.vti

# index into the sorted list of files (0 = first, 1 = second, etc.)
index = 0

if which_set == "open":
    pattern = "part2_open_t*.vti"
elif which_set == "closed":
    pattern = "part2_closed_t*.vti"
else:
    raise ValueError("which_set must be 'open' or 'closed'")

files = sorted(glob.glob(pattern))

if not files:
    raise FileNotFoundError(
        f"No files found matching pattern '{pattern}'. "
        "Make sure you've run the CN script and are in the right folder."
    )

print("Found VTI files:")
for i, fname in enumerate(files):
    print(f"  [{i}] {fname}")

if index < 0 or index >= len(files):
    raise IndexError(
        f"index={index} out of range for {len(files)} files. "
        "Pick an index from the list above."
    )

filename = files[index]
print(f"\nLoading: {filename}")

# ----------------------------------------------------------
# Load and plot
# ----------------------------------------------------------
grid = pv.read(filename)

plotter = pv.Plotter()
plotter.add_mesh(grid, cmap="viridis")
plotter.add_scalar_bar(title="Density")
plotter.show(title=filename)
