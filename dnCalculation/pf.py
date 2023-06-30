import numpy as np
import matplotlib.pyplot as plt
import dnModule


window = 1000
binNumber = 20

dn_array = np.load('./testNico.out.npy')
print(dn_array.shape)
MSD_arr = dnModule.CalcMSD(dn_array, window, binNumber)
pf_matrix = np.zeros((binNumber, binNumber))
for i in range(binNumber):
    for j in range(binNumber):
        b = dnModule.estimate_coef(np.arange(window), MSD_arr[:, i, j])
        pf_matrix[i, j] = (b[1] * 1e12/2) * dnModule.volume_w
np.save(f"pfMatrixRef.npy", pf_matrix)