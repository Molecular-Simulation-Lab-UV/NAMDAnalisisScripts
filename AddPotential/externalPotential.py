"""
Currently, it only includes an external electric field's potential difference in the Z axis.
TODO: add the axis option
"""

import gridData, argparse, numpy

parser = argparse.ArgumentParser(description = "Add the effect of an external potential to the .dx file")
parser.add_argument('-i', '--in_file', type = str, required = True, help = '[String] Path, either absolute or relative, to the dx file')
parser.add_argument('-f', '--field', type = float, required = True, help = '[Float] Value, in mV, of the external potential. E.g: -100')
parser.add_argument('-t', '--temp', type = float, required = False, default = 300, help = '[Float] Temperature, in Kelvin, at which the simulation is run. Used to calculate potential in k_BT/e. Default = 300')
parser.add_argument('-o', '--out', type = str, required = False, default = 'outputDX.dx', help = '[String] Output path of the resulting modified DX file. ./outputDX.dx by default')

arg = parser.parse_args()
dx = gridData.Grid(arg.in_file)
finalDX = dx.grid.copy()

cF = 1/(0.0862*arg.temp) # Conversion factor @ arg.temp Kelvin

if arg.field > 0:
    # FALTA TRANSFORMARLO A UNIDADES NAMD
    extPotential = numpy.linspace(0, arg.field*cF, dx.grid.shape[-1])
    finalDX = finalDX + extPotential
    
elif arg.field < 0:
    extPotential = numpy.arange(arg.field*cF, 0, dx.grid.shape[-1])
    finalDX = finalDX + extPotential

else:
    print('Nothing to do here')
    exit()

# Make new grid and write data out
dxOut = gridData.Grid(finalDX, edges=dx.edges, origin=dx.origin, delta=dx.delta)
dxOut.export(arg.out)
