"""
Currently, it only includes an external electric field's potential difference in the Z axis.
TODO: add the axis option
"""

import gridData, argparse, numpy as np

parser = argparse.ArgumentParser(description="Add the effect of an external potential to the .dx file")
parser.add_argument('-i', '--in_file', type=str, required=True, help='[String] Path, either absolute or relative, to the dx file')
parser.add_argument('-f', '--file', type=str, required=False, help='[String] Path to a detailed parameters file, in case a subset of the DX volume is to be calculated. See example')
parser.add_argument('-o', '--out', type=str, required=False, default='outputDX.dx', help='[String] Output path of the resulting modified DX file. ./outputDX.dx by default')

arg = parser.parse_args()
dx = gridData.Grid(arg.in_file)

if arg.file is not None:
    in_file = open(arg.file, 'r')
    for line in in_file:
        l = line.strip().split()
        if 'rad' in l[0].lower():
            rad = float(l[1])
        elif 'center' in l[0].lower():
            center = np.array(l[1:], dtype=float)
        elif 'len' in l[0].lower():
            length = float(l[1])
    in_file.close()

data = dx.grid
x_grid, y_grid, z_grid = np.meshgrid(*dx.edges)
z_shape = data.shape[2]

# TODO: use grid edges for positioning of cylinder

if 'rad' not in locals():
    if 'length' not in locals():
        avg = np.zeros(z_shape, dtype=float)
        std = np.zeros(z_shape, dtype=float)
        for plane in range(z_shape):
            avg[plane] = data[:, :, plane].mean()
            std[plane] = data[:, :, plane].std()
    else:
        lowerZ = center[2] - length/2
        upperZ = center[2] + length/2
        idx = np.argwhere((lowerZ < z_grid) & (upperZ > z_grid))

        avg = np.zeros_like(idx, dtype=float)
        std = np.zeros_like(idx, dtype=float)
        for p, plane in enumerate(idx):
            avg[p] = data[:, :, plane].mean()
            std[p] = data[:, :, plane].std()

else:
    if 'length' not in locals():
        idx = np.argwhere((x_grid - center[0]) ** 2 + (y_grid - center[1]) ** 2 <= rad ** 2)

        avg = np.zeros(z_shape, dtype=float)
        std = np.zeros(z_shape, dtype=float)
        # TODO: Check the out-of-bounds error
        for p, plane in enumerate(idx):
            avg[plane] = data[plane].mean()
            std[plane] = data[plane].std()
    else:
        lowerZ = center[2] - length/2
        upperZ = center[2] + length/2
        idx = np.argwhere(((x_grid - center[0]) ** 2 + (y_grid - center[1]) ** 2 <= rad ** 2) & (lowerZ <= z_grid) & (upperZ >= z_grid))
        idx_z = np.unique(idx[:, 2])

        avg = np.zeros_like(idx_z, dtype=float)
        std = np.zeros_like(idx_z, dtype=float)
        for p, plane in enumerate(idx_z):
            flag = plane == idx[:,2]
            idx_dummy = idx[flag]
            avg[p] = data[idx_dummy[:, 0], idx_dummy[:, 1], idx_dummy[:,2]].mean()
            std[p] = data[idx_dummy[:, 0], idx_dummy[:, 1], idx_dummy[:,2]].std()

with open(arg.out, 'w+') as f:
    f.write('#Z position \t Mean value \t Standard deviation\n')
    z = z_grid[0, 0, idx_z] if 'idx_z' in locals() else np.unique(z_grid)
    for line in range(avg.shape[0]):
        f.write('{0:8.3f} \t {1:8.5f} \t {2:8.5f}\n'.format(z[line] + dx.delta[2]/2, avg[line], std[line]))

