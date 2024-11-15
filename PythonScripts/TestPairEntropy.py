import numpy as np
import multiprocessing 
from matplotlib import pyplot as plt
from functools import partial
import os
from scipy.optimize import curve_fit

def dvr_k(l):

    num_grid_pts = 2 * l + 1
    precision = np.float64
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

def dvr_v(l, g):
     
    num_grid_pts = 2 * l + 1
    precision = np.float64
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


def ed_dvr_entropy(g, l):

    k = dvr_k(l)
    v = dvr_v(l, g)
    num_grid_pts = 2 * l + 1
    identity_operator = np.eye(num_grid_pts)

    k_4_body = np.kron(k, np.kron(np.kron(identity_operator, identity_operator), identity_operator)) + \
                np.kron(identity_operator, np.kron(np.kron(k, identity_operator), identity_operator)) + \
                np.kron(identity_operator, np.kron(np.kron(identity_operator, k), identity_operator)) + \
                np.kron(identity_operator, np.kron(np.kron(identity_operator, identity_operator), k))
    
    v_4_body = np.kron(v, np.kron(identity_operator, identity_operator)) + \
                np.kron(identity_operator, np.kron(v, identity_operator)) + \
                np.kron(identity_operator, np.kron(identity_operator, v))
    
    h_4_body = k_4_body + v_4_body
    evals, evecs = np.linalg.eigh(h_4_body)
    e0_vec = evecs[:,0]
    print("E0 = ", evals[0])
    
    s2 = 0.0
    for i1 in range(num_grid_pts):
        for i2 in range(num_grid_pts):
            i = i1 * num_grid_pts + i2
            for j1 in range(num_grid_pts):
                for j2 in range(num_grid_pts):
                    j = j1 * num_grid_pts + j2
                    term_1 = e0_vec[i*num_grid_pts**2 + j]
                    for k1 in range(num_grid_pts):
                        for k2 in range(num_grid_pts):
                            k = k1 * num_grid_pts + k2
                            term_2 = e0_vec[k*num_grid_pts**2 + j]
                            for l1 in range(num_grid_pts):
                                for l2 in range(num_grid_pts):
                                    l = l1 * num_grid_pts + l2
                                    s2 += e0_vec[i * num_grid_pts**2 + l] * e0_vec[k * num_grid_pts**2 + l] * term_1 * term_2

 
    return s2

def entropy_NMM(l, g, T, P):

    k_1_body = dvr_k(l)
    v_2_body = dvr_v(l, g)
    tau = 1.0/(T * P)
    h_2_body = np.kron(k_1_body, np.eye(2 * l + 1)) + np.kron(np.eye(2 * l + 1), k_1_body) + v_2_body
    evals_2_body, evecs_2_body = np.linalg.eigh(h_2_body)

    print("Ground State Energy ED (N=2) = ", evals_2_body[0])
    prop_h_2_diag = np.zeros(np.shape(h_2_body))
    for i in range(len(evals_2_body)):
        prop_h_2_diag[i,i] = np.exp(-tau * evals_2_body[i])
    
    prop_h_2 = np.matmul(evecs_2_body, np.matmul(prop_h_2_diag, np.transpose(evecs_2_body)))

    prop_k_diag = np.zeros(np.shape(k_1_body))
    evals_k, evecs_k = np.linalg.eigh(k_1_body)
    for i in range(len(evals_k)):
        prop_k_diag[i,i] = np.exp(-tau * evals_k[i])

    prop_k = np.matmul(np.matmul(evecs_k, prop_k_diag), np.transpose(evecs_k))
    
    prop_2_k = np.zeros(np.shape(h_2_body))
    num_grid_pts = 2*l+1
    for i in range(len(evals_k)):
        for j in range(len(evals_k)):
            for k in range(len(evals_k)):
                for r in range(len(evals_k)):
                    prop_2_k[i*num_grid_pts+j, k*num_grid_pts+r] = prop_k[i,k] * prop_k[j,r]

    '''
    effective_T = str(1.0/tau)[:3]
    free_rho_filename = "Free_Propagator_DVR_l_{}_g_{}_T_{}.csv".format(str(l), str(g)[:3], effective_T)
    free_rho_path = os.path.join("/Users/shaeermoeed/Github/DVRPairMC/ProbabilityDensities", free_rho_filename)
    prop_free_file = np.loadtxt(free_rho_path, delimiter=",", skiprows=2)
    print(np.array_equiv(prop_free_file, prop_k))
    '''
    print("Constructing Rho Pair")
    prop_pair = np.divide(prop_h_2, prop_2_k)

    '''
    effective_T = str(1.0/tau)[:3]
    pair_filename = "Pair_Propagator_DVR_l_{}_g_{}_T_{}.csv".format(str(l), str(g)[:3], effective_T)
    pair_rho_path = os.path.join("/Users/shaeermoeed/Github/DVRPairMC/ProbabilityDensities", pair_filename)
    prop_pair_file = np.loadtxt(pair_rho_path, delimiter=",", skiprows=2)
    print(np.array_equiv(prop_pair_file, prop_pair))
    '''

    print("Constructing Rho tau")
    rho_tau = np.ones((num_grid_pts**4, num_grid_pts**4))
    for i in range(num_grid_pts):
        for j in range(num_grid_pts):
            for k in range(num_grid_pts):
                for l in range(num_grid_pts):
                    for ip in range(num_grid_pts):
                        for jp in range(num_grid_pts):
                            for kp in range(num_grid_pts):
                                for lp in range(num_grid_pts):
                                    index = (i * (num_grid_pts**3)) + (j * num_grid_pts**2) + (k * num_grid_pts) + l
                                    index_p = (ip * (num_grid_pts**3)) + (jp * num_grid_pts**2) + (kp * num_grid_pts) + lp
                                    rho_tau[index, index_p] *= prop_k[i,ip] * prop_k[j,jp] * prop_k[k,kp] * prop_k[l,lp] * prop_pair[i * num_grid_pts + j, ip * num_grid_pts + jp] * prop_pair[j * num_grid_pts + k, jp * num_grid_pts + kp] * prop_pair[k * num_grid_pts + l, kp * num_grid_pts + lp] 
    
    print("NMM")
    rho_beta_over_2 = np.copy(rho_tau)
    delta_phi = 2.0 * np.pi/num_grid_pts
    for i in range(int(P/2) - 1):
        #rho_beta_over_2 = (delta_phi**4) * np.matmul(rho_beta_over_2, rho_tau)
        rho_beta_over_2 = np.matmul(rho_beta_over_2, rho_tau)
    #rho_beta = (delta_phi**4) * np.matmul(rho_beta_over_2, rho_beta_over_2)
    rho_beta = np.matmul(rho_beta_over_2, rho_beta_over_2)

    print("Calculating rho_AB")
    psi_T_psi_T = np.ones(np.shape(rho_beta_over_2))
    #rho_AB =(delta_phi**8) * np.matmul(rho_beta_over_2, np.matmul(psi_T_psi_T, rho_beta_over_2))
    rho_AB = np.matmul(rho_beta_over_2, np.matmul(psi_T_psi_T, rho_beta_over_2))

    print("Calculating Z_0")
    #Z_0 = (delta_phi**8) * np.sum(rho_beta) 
    Z_0 = np.sum(rho_beta) 
    
    """ Z_0 = 0.0
    for i in range(num_grid_pts**4):
        for j in range(num_grid_pts**4):
            term_1 = 0.0
            term_2 = 0.0
            Z_0 += rho_beta[i,j]
            for k in range(num_grid_pts**4):
                term_1 += rho_beta_over_2[i, k]
                term_2 += rho_beta_over_2[j, k]
            rho_AB[i,j] = term_1 * term_2
     """
    
    print("Calculating S2")
    S2 = 0.0
    for i1 in range(num_grid_pts):
        for i2 in range(num_grid_pts):
            i = i1 * num_grid_pts + i2
            for j1 in range(num_grid_pts):
                for j2 in range(num_grid_pts):
                    j = j1 * num_grid_pts + j2
                    for k1 in range(num_grid_pts):
                        for k2 in range(num_grid_pts):
                            k = k1 * num_grid_pts + k2
                            term = rho_AB[i * (num_grid_pts**2) + j, k * (num_grid_pts**2) + j]
                            for l1 in range(num_grid_pts):
                                for l2 in range(num_grid_pts):
                                    l = l1 * num_grid_pts + l2
                                    S2 += term * rho_AB[k * (num_grid_pts**2) + l, i * (num_grid_pts**2) + l]
    S2 /= (Z_0**2)
    #S2 *= (delta_phi ** 8)

    S2 = -np.log(S2)            
    print("P = ", P)
    print("S2 = ", S2)

    return S2

def func_cubic(x, a, b):

    return a*x*x*x + b

def func_quadratic(x, a, b):

    return a*x*x + b
    
def extrapolate_results(x_pts, y_pts, func, num_pts_out):
        
        y_pts_to_fit = np.exp(-y_pts)

        coefficients, cov = curve_fit(func, x_pts, y_pts_to_fit)
        err = np.sqrt(np.diag(cov))[1]

        fit_x = []
        fit_y = []
        for i in range(num_pts_out):
            x = max(x_pts) * (1.0 - (float(i)/float(num_pts_out-1.0)))
            y = func(x, coefficients[0], coefficients[1])
            fit_x.append(x)
            fit_y.append(y)

        err = err/fit_y[-1]

        return fit_x, -np.log(fit_y), err

if __name__ == "__main__":

    '''
    P_list = [80, 100]
    pool = multiprocessing.Pool()
    pool.map(partial(entropy_NMM, 4, 1.0, 0.1), P_list)
    '''
    
    entropy_matrix_results = "/Users/shaeermoeed/Github/DVRPairMC/entropy_ed_test.txt"
    with open(entropy_matrix_results,'r') as f:
        line_count = 0
        lines = f.readlines()
        ed_line = lines[1]
        data_ed = ed_line.strip().split("=")
        S2_ed = float(data_ed[1])
        dmrg_line = lines[2]
        data_dmrg = dmrg_line.strip().split("=")
        S2_dmrg = float(data_dmrg[1])

    P_vals, S2_NMM = np.loadtxt(entropy_matrix_results, delimiter=",", skiprows=5, unpack=True)
    tau_list = [1.0/(0.1 * P) for P in P_vals]
    tau_list = tau_list
    S2_NMM = S2_NMM
    print(S2_NMM)
    print(tau_list)
    fit_tau_cubic, fit_S2_cubic, S2_err_cubic = extrapolate_results(tau_list, S2_NMM, func_cubic, 100)
    fit_tau_quadratic, fit_S2_quadratic, S2_err_quadratic = extrapolate_results(tau_list, S2_NMM, func_quadratic, 100)

    mc_results = "/Users/shaeermoeed/Github/DVRPairMC/Results/Entropy/Test/12_06_2024_07_16_13/Estimator Statistics.csv"
    g,T,P,tau_mc,S2_mc,S2_mc_err = np.loadtxt(mc_results, delimiter=',', skiprows=2, unpack=True)

    tau_mc, S2_mc = zip(*sorted(zip(tau_mc, S2_mc), reverse=True))
    tau_mc, S2_mc_err = zip(*sorted(zip(tau_mc, S2_mc_err), reverse=True))
    tau_mc = tau_mc
    S2_mc = np.array(S2_mc)
    S2_mc_err = S2_mc_err

    fit_tau_mc, fit_S2_mc, fit_S2_err_mc = extrapolate_results(tau_mc, S2_mc, func_cubic, 100)

    plt.figure()
    plt.scatter(tau_list, S2_NMM, label='NMM', color="C0")
    plt.scatter(tau_mc, S2_mc, label="MC", color='C3')
    #plt.errorbar(tau_mc, S2_mc, S2_mc_err, fmt='None', capsize=5, color='black')
    plt.plot(tau_list, S2_NMM, label="NMM", color='C0')
    #plt.plot(fit_tau_cubic, fit_S2_cubic, label="Cubic Fit")
    plt.plot(np.linspace(min(tau_list), max(tau_list), 50), [S2_ed]*len(np.linspace(min(tau_list), max(tau_list), 50)), ".", label="ED", color='C2')
    plt.plot(np.linspace(min(tau_list), max(tau_list), 50), [S2_dmrg]*len(np.linspace(min(tau_list), max(tau_list), 50)), label="DMRG", color='C4')
    #plt.errorbar([0.0], [fit_S2_cubic[-1]], S2_err_cubic, label="Cubic Fit Error", capsize=5)
    #plt.plot(fit_tau_mc, fit_S2_mc, label="MC Fit")
    #plt.errorbar([0.0], [fit_S2_mc[-1]], fit_S2_err_mc, label="MC Fit Error", capsize=5)
    plt.title("S2 Tau Convergence")
    plt.xlabel("Tau")
    plt.ylabel("S2")
    #plt.ylim(0.6050,0.6150)
    plt.legend()
    plt.savefig("/Users/shaeermoeed/Github/DVRPairMC/S2 NMM Results")
    
    