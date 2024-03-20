import numpy as np
import multiprocessing 
from matplotlib import pyplot as plt
from functools import partial
import os

precision_dict = {0:np.float64, 1:np.float32}

def dvr_k(l, precision_enum):

    num_grid_pts = 2 * l + 1
    precision = precision_dict[precision_enum]
    K_DVR = np.zeros((num_grid_pts, num_grid_pts), dtype=precision)
    l_new = l
    num_g_p_new = 2 * l_new + 1
    for i in range(num_grid_pts):
            for j in range(num_grid_pts):
                if i == j:
                    K_DVR[i,i] = (l_new**2 + l_new)/3
                else:
                    index_diff = i-j
                    sign = ((-1.0)**(index_diff))
                    angular_arg = np.pi * index_diff/num_g_p_new
                    K_DVR[i,j] =  sign * (np.cos(angular_arg)/(2*(np.sin(angular_arg))**2))

    return K_DVR

def dvr_v(l, g, precision_enum):
     
    num_grid_pts = 2 * l + 1
    precision = precision_dict[precision_enum]
    V_DVR = np.zeros((num_grid_pts**2, num_grid_pts**2), dtype=precision)
    for i in range(num_grid_pts):
        for j in range(num_grid_pts):
            grid_pt_1 = 2 * np.pi * i/num_grid_pts
            grid_pt_2 = 2 * np.pi * j/num_grid_pts
            s1 = np.sin(grid_pt_1)
            s2 = np.sin(grid_pt_2)
            c1 = np.cos(grid_pt_1)
            c2 = np.cos(grid_pt_2)
            V_DVR[i * num_grid_pts + j, i * num_grid_pts + j] += g * ((s1 * s2) - (2 * c1 * c2))

    return V_DVR

def propagator_k_marx(l, P, precision_num, output_diag):

    precision = precision_dict[precision_num]
    num_grid_pts = 2 * l + 1
    delta_phi = 2.0 * np.pi/num_grid_pts
    free_rho_marx = np.zeros((num_grid_pts, num_grid_pts), dtype=precision)
    tau = 1/(0.1 * P)
    m_max = l
    for i in range(num_grid_pts):
        for j in range(num_grid_pts):
            dphi = float(i - j) * delta_phi
            integral = 0.
            for m in range(m_max+1):
                integral += np.exp(-1./(4.*tau)*(dphi+2.*np.pi*float(m))**2)
            for m in range(1,m_max+1):
                integral+=np.exp(-1./(4.*tau)*(dphi+2.*np.pi*float(-m))**2)
            integral*=np.sqrt(1./(4.*np.pi*tau))
            free_rho_marx[i,j] = integral

    if output_diag: 
        print("Free rho marx = ")
        print(free_rho_marx)

    rho_2_k = np.kron(free_rho_marx, free_rho_marx)

    return free_rho_marx, rho_2_k


def ed_dvr(g, l):

    k = dvr_k(l, 0)
    v = dvr_v(l, g, 0)
    h = np.kron(k, np.eye(2 * l + 1)) + np.kron(np.eye(2 * l + 1), k) + v
    evals, evecs = np.linalg.eigh(h)
    return evals[0]

def pair_prop_test(g, l, P):

    k = dvr_k(l, 0)
    v = dvr_v(l, g, 0)
    tau = 1.0/(0.1 * P)
    T = 1.0/tau
    h = np.kron(k, np.eye(2 * l + 1)) + np.kron(np.eye(2 * l + 1), k) + v
    evals, evecs = np.linalg.eigh(h)
    prop_h_2_diag = np.zeros(np.shape(h))
    for i in range(len(evals)):
        prop_h_2_diag[i,i] = np.exp(-tau * evals[i])
    
    prop_h_2 = np.matmul(evecs, np.matmul(prop_h_2_diag, np.transpose(evecs)))

    prop_k_diag = np.zeros(np.shape(k))
    evals_k, evecs_k = np.linalg.eigh(k)
    for i in range(len(evals_k)):
        prop_k_diag[i,i] = np.exp(-tau * evals_k[i])

    prop_k = np.matmul(np.matmul(evecs_k, prop_k_diag), np.transpose(evecs_k))

    prop_2_k = np.kron(prop_k, prop_k)
    num_grid_pts = 2 * l + 1

    prop_pair = np.divide(prop_h_2, prop_2_k)
    for i in range(num_grid_pts**2):
        for j in range(num_grid_pts**2):
            if prop_pair[i,j] < 0.0:
                print("Negative pair density element for T = ", T)
                return 0
            if prop_2_k[i, j] < 0.0:
                print("Negative free density element for T = ", T)
                return 0
            
    delta_phi = 2.0 * np.pi/num_grid_pts

    rho_tau = np.multiply(prop_2_k, prop_pair)
    rho_beta=rho_tau.copy()

    for k in range(P-1):
        rho_beta=(delta_phi**2)*np.matmul(rho_beta,rho_tau)

    E0_nmm=0.0
    potential_vec = np.diag(v)
    rho_dot_V=np.matmul(rho_beta,potential_vec)
    Z0=0.0
    for i in range(num_grid_pts**2):
        E0_nmm += rho_dot_V[i]
        for ip in range(num_grid_pts**2):
            Z0 += rho_beta[i,ip]
    E0_nmm/=Z0

    return E0_nmm

def write_densities(parent_dir, l, g_list, test_props, T):

    k = dvr_k(l, 0)

    for i in range(len(g_list)):
        g = g_list[i]
        v = dvr_v(l, g, 0)
        tau = 1.0/T
        h = np.kron(k, np.eye(2 * l + 1)) + np.kron(np.eye(2 * l + 1), k) + v
        evals, evecs = np.linalg.eigh(h)
        prop_h_2_diag = np.zeros(np.shape(h))
        for i in range(len(evals)):
            prop_h_2_diag[i,i] = np.exp(-tau * evals[i])
        
        prop_2_h = np.matmul(evecs, np.matmul(prop_h_2_diag, np.transpose(evecs)))

        prop_k_diag = np.zeros(np.shape(k))
        evals_k, evecs_k = np.linalg.eigh(k)
        for i in range(len(evals_k)):
            prop_k_diag[i,i] = np.exp(-tau * evals_k[i])

        prop_k = np.matmul(np.matmul(evecs_k, prop_k_diag), np.transpose(evecs_k))

        prop_2_k = np.kron(prop_k, prop_k)
        num_grid_pts = 2 * l + 1

        prop_pair = np.divide(prop_2_h, prop_2_k)
        density_negative = False
        if test_props:
            for i in range(num_grid_pts**2):
                for j in range(num_grid_pts**2):
                    if prop_pair[i,j] < 0.0:
                        print("Negative pair density element for T = {}, g = {}".format(T,g))
                        density_negative = True
                    if prop_2_k[i, j] < 0.0:
                        print("Negative free density element for T = {}, g = {}".format(T,g))
                        density_negative = True

        if not density_negative:
            write_file(parent_dir, "Pair", prop_pair, l, T, g)
            write_file(parent_dir, "Free", prop_k, l, T, g)

    return 0

def write_file(parent_dir, density_name, density_matrix, l, T, g):

    file_name = "{}_Propagator_DVR_l_{}_g_{}_T_{}.csv".format(density_name, l, g, T)
    filepath = os.path.join(parent_dir, file_name)

    num_grid_pts = 2*l+1
    if density_name == "Pair":
        dims = num_grid_pts**2
    elif density_name == "Free":
        dims = num_grid_pts

    with open(filepath, 'w') as f:
        header = "{} Density DVR\n".format(density_name)
        header += "l = {}".format(l) + " g = {}".format(g) + " T = {}".format(T) + "\n"
        f.write(header)
        for i in range(dims):
            write_str = str(density_matrix[i,0])
            for j in range(1, dims):
                write_str += "," + str(density_matrix[i,j])
            f.write(write_str + "\n")

    return 0

if __name__ == "__main__":

    T_list = [0.5, 1.0, 1.5, 2.0, 6.0]
    tau_list = [1.0/(T) for T in T_list]
    l = 10
    g_list = [0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0]
    N = 2
    pool = multiprocessing.Pool()
    parent_dir = "/Users/shaeermoeed/Github/DVRPairMC/ProbabilityDensities"
    pool.map(partial(write_densities, parent_dir, l, g_list, True), T_list)
    '''
    e0_list = pool.map(partial(pair_prop_test, g, l), P_list)
    print(e0_list)
    ed_e0 = -0.5292920195323415
    plt.figure()
    plt.plot(tau_list, [ed_e0]*len(P_list))
    plt.scatter(tau_list, e0_list)
    plt.show()
    '''

