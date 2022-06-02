import numpy
import matplotlib.pyplot as plt
import linecache

def readDX(filePath):

    """
    This script reads the header of the .dx file and stores the important variables:
    - numStep:      number of grid elements in format Nx Ny Nz.
    - origin:       physical origin of the system. One of the corners of the system box.
    - delta:        array containing the 3 "step sizes", the distance between each grid point in each direction.
    - lineNumber:   integer holding the line number where the actual data begins
     """

    file = open(filePath, 'r')
    header = file.readlines()[:12]
    deltas = []
    for i, line in enumerate(header):
        if "gridpositions" in line:
            dummy = line.split()
            numSteps = numpy.array(dummy[-3:]).astype('int64')

        elif "origin" in line:
            dummy = line.split()
            origin = numpy.array(dummy[-3:]).astype('float64')

        elif "delta" in line:
            dummy = line.split()
            for val in dummy:
                try:
                    float(val)
                    if float(val) != 0:
                        deltas.append(float(val))
                except:
                    continue

        elif "data" in line:
            lineNumber = i + 1
            break

    deltas = numpy.array(deltas)
    file.close()
    return origin, deltas, numSteps, lineNumber

def calcMean(filePath):
    """
    Script averaging the data in the dx file for each plane z = something. It starts at z = origin[2], advances with step delta[2], and stops at origin[2] + numSteps[2]*delta[2]
    The script returns the average value of the dx file for all points in a X-Y plane, for a given z.    
    """
    origin, deltas, numSteps, lineNum = readDX(filePath)
    values = numpy.zeros(numSteps[-1])
    dataAmount = numSteps[0]*numSteps[1]*numSteps[2]
    rawData = numpy.loadtxt(filePath, skiprows=lineNum, max_rows = int(dataAmount//3)) # Read the data of every full row (full = 3 values)
    if float(rawData.shape[0]) != dataAmount/3: # Check whether there is an incomplete row (2 or 1 value in the last row)
        data = numpy.reshape(rawData, (rawData.shape[0]*rawData.shape[1])) # Reshape the array to make a 1D vector with the all the full rows' data
        data = numpy.concatenate((data, numpy.array(linecache.getline(filePath, lineNum + dataAmount//3 + 1).split()).astype('float64'))) # Append the last incomplete row's data, if it exists
    else:
        data = numpy.reshape(rawData, dataAmount) # If there is no incomplete row, work as is
    for i in range(numSteps[0]*numSteps[1] - 1):
        values += data[i*numSteps[2]:(i+1)*numSteps[2]] # Add every x,y point for the same value of z
    avgs = values/(numSteps[0]*numSteps[1])
    return avgs


p = '/home/nespinoza/CINV_Research/Colaboraciones_Aportes/Guido/dxProcessing/mapa2us.dx'
avs = calcMean(p)
dum1, dum2, Steps, dum4 = readDX(p)
plt.plot(avs, numpy.linspace(dum1[2], (dum2[2]*Steps[2] - abs(dum1[2])), Steps[2]))
plt.show()
