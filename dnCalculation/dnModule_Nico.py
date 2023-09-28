
#################
#### IMPORTS ####
#################

import prody
import numpy as np
import matplotlib.pyplot as plt


###################
#### CONSTANTS ####
###################

volume_w = 18 / 6.02214e23


###################
#### FUNCTIONS ####
###################

def dnMatrixCalculation(trajectory, structure, selName, zmin, zmax, radius2, binSize, refName):
    """
    All-in-one binning + dn matrix calculation for each frame in a trajectory.

    Parameters
    ----------
    trajectory : prody trajectory object
    structure : prody atomgroup object
    selName: A selection of atoms to be tracked. In ion channels it's the ion of interest.
        In AQPs, "name OH2" should do.
    zmin & zmax : z coordinates of the cylinder limits, asumming its center in (0, 0, 0)
    radius2 : square of the cylinder's radius
    binSize : size of the cylinder's bins, all equal, in [A] and along the length of the cylinder (z-direction)
    refName : Same as used for alignment. It's needed for wrapping the simulation system.

    Returns
    -------
    dnMatrix : A len(trajectory) x N_bins (number of bins) array containing historic
        information about dn at each point in time during the simulation.

    """
    binArray = np.arange(zmin, zmax + binSize, binSize) # This implementation considers fixed bins' array size and cylinder measurements
    dnMatrix = np.zeros((len(trajectory), len(binArray)-1))
    # First frame, to have a frame i-1 with which we can operate.
    f0 = trajectory.next()
    prody.wrapAtoms(structure, unitcell = f0.getUnitcell()[:3], center = prody.calcCenter(structure.select(refName)))
    f0.superpose()
    # Initializers for positions, indices and whether it's in bin or not at frame 0
    oldPos = structure.select(f'{selName} and (x^2 + y^2) < {radius2} and z > {zmin} and z < {zmax}').getCoords()[:,-1]
    oldInd = structure.select(f'{selName} and (x^2 + y^2) < {radius2} and z > {zmin} and z < {zmax}').getIndices()
    oldInBin = np.argwhere((oldPos[:,np.newaxis] >= binArray[np.newaxis,:-1]) & (oldPos[:,np.newaxis] < binArray[np.newaxis,1:]))
    for frame_n, frame in enumerate(trajectory, start = 1):
        prody.wrapAtoms(structure, unitcell = frame.getUnitcell()[:3], center = prody.calcCenter(structure.select(refName)))
        frame.superpose()
        sel = structure.select(f'{selName} and (x^2 + y^2) < {radius2} and z > {zmin} and z < {zmax}')
        pos = sel.getCoords()[:,-1] # Grab selection position, z-coordinate
        ind = sel.getIndices() # Grab selection (atom) indices
        # Defines which INDEX of the pos/ind array goes into which bin.
        inBin = np.argwhere((pos[:,np.newaxis] >= binArray[np.newaxis,:-1]) & (pos[:,np.newaxis] < binArray[np.newaxis,1:])) 
        horizontalStack = np.vstack((ind, inBin[:,1])).T # Combine atom index information with bin's index information
        oldHorizontalStack = np.vstack((oldInd, oldInBin[:,1])).T # Same as above, for previous frame (old) data
        # Mask the combined information array to leave only the ARRAY INDICES of the atomic indices + bin indices; (atom_index, bin_index)
        # For example, if only atom index 1000 is in the same bin as it was in the previous frame, mask will be (1000, index_of_bin)
        mask = np.flatnonzero((oldHorizontalStack == horizontalStack[:,None]).all(-1).any(-1))
        oldMask = np.flatnonzero((horizontalStack == oldHorizontalStack[:,None]).all(-1).any(-1)) # Same, but for the previous (old) frame
        dnMatrix[frame_n, inBin[mask,1]] = pos[mask] - oldPos[oldMask]/binSize # Now we use the masks to go find the positions of the same atom, that is in the same bin, in two consecutive frames, and calculate dn.
        # It is complex, but since two frames may not have the same amount of water molecules in the cylinder and it may be in disarray, direct comparison would be incomplete or straight up wrong.

        # Update everything for the next iteration
        oldPos = pos
        oldInd = ind
        oldInBin = inBin

    return dnMatrix

def CalcMSD(arreglo, window, binNumber):
    MSD_N2 = np.zeros((window, binNumber, binNumber)) # Los ejes 1 y 2 forman la matriz simétrica de correlaciones entre los bin. La diagonal constituye los pfs por bin.
    for row in range(len(arreglo) - window):
        for i in range(binNumber):
            for j in range(binNumber):
                MSD_N2[1:, i, j] += np.cumsum(arreglo[row+1:row+window, i], axis = 0) * np.cumsum(arreglo[row+1:row+window, j], axis = 0)

    return MSD_N2/(len(arreglo) - window)


# Variante propuesta por ChatGPT. No la he probado, la función CalcMSD de más arriba hecha por Agustín ya corre rápidamente, no vale la pena perder tiempo en optimizarla.
# def CalcMSD(arreglo, window, binNumber):
#     MSD_N2 = np.zeros((window, binNumber, binNumber))
#     cumsum = np.cumsum(arreglo, axis=0)
#     cumsum_sq = np.cumsum(arreglo ** 2, axis=0)
#     for row in range(len(arreglo) - window):
#         row_end = row + window
#         MSD_N2[1:] += (cumsum[row_end:row_end+window, :, np.newaxis]
#                        - cumsum[row:row+window-1, :, np.newaxis]) * (cumsum[:, np.newaxis, :] - cumsum[row+1:row+window+1, np.newaxis, :])
#         MSD_N2[1:] += (cumsum_sq[row_end:row_end+window, :, np.newaxis]
#                        - cumsum_sq[row:row+window-1, :, np.newaxis]) * (cumsum_sq[:, np.newaxis, :] - cumsum_sq[row+1:row+window+1, np.newaxis, :])
#     return MSD_N2 / (len(arreglo) - window)




def estimate_coef(x, y):
    # number of observations/points
    n = np.size(x)

    # mean of x and y vector
    m_x = np.mean(x)
    m_y = np.mean(y)

    # calculating cross-deviation and deviation about x
    SS_xy = np.sum(y*x) - n*m_y*m_x
    SS_xx = np.sum(x*x) - n*m_x*m_x

    # calculating regression coefficients
    b_1 = SS_xy / SS_xx
    b_0 = m_y - b_1*m_x

    return (b_0, b_1)

def plot_regression_line(x, y, b):
    # plotting the actual points as scatter plot
    plt.scatter(x, y, color = "m",
                marker = "o", s = 30)

    # predicted response vector
    y_pred = b[0] + b[1]*x

    # plotting the regression line
    plt.plot(x, y_pred, color = "g")

    # putting labels
    plt.xlabel('x')
    plt.ylabel('y')

    # function to show plot
    plt.show()

