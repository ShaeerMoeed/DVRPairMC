import numpy as np
import os

def dvr_k_prop(beta, P, l):

    num_grid_pts = 2 * l + 1
    precision = np.float64
    K_DVR = np.zeros((num_grid_pts, num_grid_pts), dtype=precision)
    for i in range(num_grid_pts):
            for j in range(num_grid_pts):
                if i == j:
                    K_DVR[i,i] = (l**2 + l)/3
                else:
                    index_diff = i-j
                    sign = ((-1.0)**(index_diff))
                    angular_arg = np.pi * index_diff/num_grid_pts
                    K_DVR[i,j] =  sign * (np.cos(angular_arg)/(2*(np.sin(angular_arg))**2))

    tau = beta/P                
    prop_k1_diag = np.zeros(np.shape(K_DVR))
    prop_k2_diag = np.zeros(np.shape(K_DVR))
    evals_k, evecs_k = np.linalg.eigh(K_DVR)

    for i in range(len(evals_k)):
        prop_k1_diag[i,i] = np.exp(-tau * evals_k[i])
        prop_k2_diag[i,i] = np.exp(tau * evals_k[i])

    prop_k_1 = np.matmul(np.matmul(evecs_k, prop_k1_diag), np.transpose(evecs_k))
    prop_k_2 = np.matmul(np.matmul(evecs_k, prop_k2_diag), np.transpose(evecs_k))
    print(np.multiply(prop_k_1, prop_k_2))
    prop_k_2 = np.divide(np.ones(np.shape(K_DVR)), prop_k_2)

    print(prop_k_1)
    print(prop_k_2)

    return K_DVR

dvr_k_prop(0.1, 10, 2)

