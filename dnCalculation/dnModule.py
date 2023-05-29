
#################
#### IMPORTS ####
#################

from tqdm import tqdm
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

def getBinsData(trajectory, structure, zmin, zmax, radius, binSize):
    """"""
    trajData = []
    for frame_n, frame in enumerate(trajectory):
        # print("Analizando frame ", frame_n)
        frame.superpose()
        radius2 = radius ** 2
        # currentIndexes = structure.select(f'name OH2 and (x^2 + y^2)<{radius2} and z>{zmin} and z<{zmin + binSize}')
        # print(f"{currentIndexes.getIndices()}")
        frameData = []
        for bin_n, current_bin in enumerate(np.arange(zmin, zmax, binSize)):
            binData = {} #inicializo un diccionario que contendrá los índices de las moléculas de agua y sus coordenadas como key:values.
            binMin = current_bin #define el minimo del bin
            binMax = current_bin + binSize  #define el máximo del bin
            currentAtoms = structure.select(f'name OH2 and (x^2 + y^2)<{radius2} and z>{binMin} and z<{binMax}')
            # print(currentAtoms.getIndices())
            try: #currentAtoms.getIndices().any(): #equivale a "si hay átomos en la selección:"
                for atom, coords in zip(currentAtoms.getIndices(), currentAtoms.getCoords()):
                    # print(currentAtoms.getIndices())
                    # print(atom, coords[2])
                    binData[atom] = coords[2]
            except: pass
            frameData.append(binData)

        trajData.append(frameData)
    return trajData


def dnMatrixCalculation(BinsData, BinsNumber):
    dnMatrix = np.zeros((len(BinsData), BinsNumber))
    data = BinsData.copy()
    for frame_n in range(1, len(data)):
        for bin_n in range(BinsNumber):
            BinWaters = data[frame_n][bin_n]
            BinWatersBefore = data[frame_n - 1][bin_n]
            WatersIDs = BinWaters.keys()
            for water in WatersIDs:
                if water in BinWatersBefore:
                    dnMatrix[frame_n][bin_n] = BinWaters[water] - BinWatersBefore[water]
    return dnMatrix


def CalcMSD(arreglo, window, binNumber):
    MSD_N2 = np.zeros((window, binNumber, binNumber)) # Los ejes 1 y 2 forman la matriz simétrica de correlaciones entre los bin. La diagonal constituye los pfs por bin.
    for row in range(len(arreglo) - window):
        for i in range(binNumber):
            for j in range(binNumber):
                MSD_N2[1:, i, j] += np.cumsum(arreglo[row+1:row+window, i], axis = 0) * np.cumsum(arreglo[row+1:row+window, j], axis = 0)

    return MSD_N2/(len(arreglo) - window)



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

