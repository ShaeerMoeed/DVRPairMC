import numpy as np
import os

def dvr_k(l):

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

def dvr_v_n3(l, g):
     
    num_grid_pts = 2 * l + 1
    precision = np.float64
    V_DVR = np.zeros((num_grid_pts**3, num_grid_pts**3), dtype=precision)
    for i in range(num_grid_pts):
        grid_pt_1 = 2.0 * np.pi * i/num_grid_pts
        s1 = np.sin(grid_pt_1)
        c1 = np.cos(grid_pt_1)
        for j in range(num_grid_pts):
            grid_pt_2 = 2.0 * np.pi * j/num_grid_pts
            s2 = np.sin(grid_pt_2)
            c2 = np.cos(grid_pt_2)
            for k in range(num_grid_pts):
                grid_pt_3 = 2.0 * np.pi * k/num_grid_pts
                s3 = np.sin(grid_pt_3)
                c3 = np.cos(grid_pt_3)
                index = (i * (num_grid_pts**2)) + (j * num_grid_pts) + k
                V_DVR[index, index] += g * ((s1 * s2) - (2 * c1 * c2)) + g * ((s2 * s3) - (2 * c2 * c3))

    return V_DVR

def get_corr_vec(num_grid_pts):

    corr_vec = np.zeros(num_grid_pts**3)
    for i in range(num_grid_pts):
        grid_pt_1 = 2 * np.pi * i/num_grid_pts
        for j in range(num_grid_pts):
            grid_pt_2 = 2 * np.pi * j/num_grid_pts
            corr_1 = np.cos(grid_pt_1 - grid_pt_2)
            for k in range(num_grid_pts):
                grid_pt_3 = 2 * np.pi * k/num_grid_pts
                corr_2 = np.cos(grid_pt_2 - grid_pt_3)
                corr_vec[i*(num_grid_pts**2) + j*num_grid_pts + k] += corr_1 + corr_2

    return corr_vec

def get_binder(num_grid_pts):

    binder_vec = np.zeros(num_grid_pts**3)
    for i in range(num_grid_pts):
        grid_pt_1 = 2 * np.pi * i/num_grid_pts
        pol_1 = np.cos(grid_pt_1)
        for j in range(num_grid_pts):
            grid_pt_2 = 2 * np.pi * j/num_grid_pts
            pol_2 = np.cos(grid_pt_2)
            for k in range(num_grid_pts):
                grid_pt_3 = 2 * np.pi * k/num_grid_pts
                pol_3 = np.cos(grid_pt_3)
                binder_vec[i*(num_grid_pts**2) + j*num_grid_pts + k] += 1.0 - 3.0*(pol_1**4 + pol_2**4 + pol_3**4)/(3.0 * ((pol_1**2 + pol_2**2 + pol_3**2)**2))

    return binder_vec

def prop_n3_SOS_preprocess(l, g):

    k_1_body = dvr_k(l)

    num_grid_pts = 2*l+1
    k_3_body = np.kron(k_1_body, np.eye(num_grid_pts**2)) + np.kron(np.kron(np.eye(num_grid_pts), k_1_body),np.eye(num_grid_pts)) + np.kron(np.eye(num_grid_pts**2), k_1_body) 
    #k_3_body = np.kron(k_1_body, np.kron(np.eye(num_grid_pts), np.eye(num_grid_pts))) + np.kron(np.eye(num_grid_pts), np.kron(k_1_body, np.eye(num_grid_pts))) + np.kron(np.eye(num_grid_pts), np.kron(np.eye(num_grid_pts), k_1_body))
    v_3_body = dvr_v_n3(l, g)

    h_3_body = k_3_body + v_3_body
    evals, evecs = np.linalg.eigh(h_3_body)
    print(evals[0])

    potential_vec = np.diag(v_3_body)
    corr_vec = get_corr_vec(num_grid_pts)
    binder_vec = get_binder(num_grid_pts)

    return evals, evecs, corr_vec, potential_vec, binder_vec

def prop_n3_SOS(eigvals, eigstates, corr_vec, pot_vec, binder_vec, T, num_grid_pts):

    beta = 1.0/T
    partition_func = 0.0
    v_exp = 0.0
    corr_exp = 0.0
    binder_exp = 0.0
    energy_exp = 0.0
    for i in range(num_grid_pts**3):
        prob = np.exp(-beta * eigvals[i])
        partition_func += prob
        state_vec = eigstates[:,i]
        v_exp_term = 0.0
        corr_exp_term = 0.0
        binder_exp_term = 0.0
        for j in range(num_grid_pts**3):
            v_exp_term += (state_vec[j]**2) * pot_vec[j]
            corr_exp_term += (state_vec[j]**2) * corr_vec[j]
            binder_exp_term += (state_vec[j]**2) * binder_vec[j]
        v_exp += v_exp_term * prob
        corr_exp += corr_exp_term * prob
        binder_exp += binder_exp_term * prob
        energy_exp += eigvals[i] * prob

    v_exp /= partition_func
    corr_exp /= partition_func
    binder_exp /= partition_func
    energy_exp /= partition_func

    ground_state = eigstates[:, 0]
    V_operator = np.diag(pot_vec)
    corr_operator = np.diag(corr_vec)
    binder_operator = np.diag(binder_vec)
    gs_exp_v = np.dot(np.dot(np.transpose(ground_state), V_operator), ground_state)
    gs_exp_corr = np.dot(np.dot(np.transpose(ground_state), corr_operator), ground_state)
    gs_exp_binder = np.dot(np.dot(np.transpose(ground_state), binder_operator), ground_state)

    print("GS <V> = ", gs_exp_v)
    print("GS <Corr> = ", gs_exp_corr)
    print("GS <Binder> = ", gs_exp_binder)

    print("Ground State Energy ED (N=3) = ", eigvals[0])

    return v_exp, corr_exp, binder_exp, energy_exp

def prop_n3_NMM(l, g, P, T, corr_vec, potential_vec, binder_vec):

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

    prop_pair = np.divide(prop_h_2, prop_2_k)

    '''
    effective_T = str(1.0/tau)[:3]
    pair_filename = "Pair_Propagator_DVR_l_{}_g_{}_T_{}.csv".format(str(l), str(g)[:3], effective_T)
    pair_rho_path = os.path.join("/Users/shaeermoeed/Github/DVRPairMC/ProbabilityDensities", pair_filename)
    prop_pair_file = np.loadtxt(pair_rho_path, delimiter=",", skiprows=2)
    print(np.array_equiv(prop_pair_file, prop_pair))
    '''

    rho_tau = np.ones((num_grid_pts**3, num_grid_pts**3))
    for i in range(num_grid_pts):
        for j in range(num_grid_pts):
            for k in range(num_grid_pts):
                for ip in range(num_grid_pts):
                    for jp in range(num_grid_pts):
                        for kp in range(num_grid_pts):
                            index = (i * (num_grid_pts**2)) + (j * num_grid_pts) + k
                            index_p = (ip * (num_grid_pts**2)) + (jp * num_grid_pts) + kp
                            rho_tau[index, index_p] *= prop_k[i,ip] * prop_k[j,jp] * prop_k[k,kp] * prop_pair[i * num_grid_pts + j, ip * num_grid_pts + jp] * prop_pair[j * num_grid_pts + k, jp * num_grid_pts + kp] 
    
    rho_beta = np.copy(rho_tau)
    delta_phi = 2.0 * np.pi/num_grid_pts
    for i in range(P-1):
        rho_beta = (delta_phi**3) * np.matmul(rho_beta, rho_tau)

    v_exp = 0.0
    corr_exp = 0.0
    binder_exp = 0.0
    partition_func = 0.0
    for i in range(num_grid_pts**3):
        v_exp += rho_beta[i,i]*potential_vec[i]
        corr_exp += rho_beta[i,i]*corr_vec[i]
        binder_exp += rho_beta[i,i]*binder_vec[i]
        partition_func += rho_beta[i,i]

    v_exp_check_nmm = np.trace(np.matmul(np.diag(potential_vec), rho_beta))
    Z_check = np.trace(rho_beta)
    v_exp_check_nmm /= Z_check
    print(v_exp_check_nmm)

    v_exp *= 1.0 / partition_func
    corr_exp *= 1.0 / partition_func
    binder_exp *= 1.0 / partition_func

    return v_exp, corr_exp, binder_exp

if __name__ == "__main__":

    l = 5
    g = 1.0
    P = 4
    T = 1.0

    evals, evecs, corr_vec, potential_vec, binder_vec = prop_n3_SOS_preprocess(l, g)
    v_exp_sos, corr_exp_sos, binder_exp_sos, energy_exp_sos = prop_n3_SOS(evals, evecs, corr_vec, potential_vec, binder_vec, T, 2*l+1)

    #v_3_body = dvr_v_n3(l, g)
    #potential_vec = np.diag(v_3_body)
    #corr_vec = get_corr_vec(2 * l + 1)
    #binder_vec = get_binder(2 * l + 1)
    v_exp_nmm, corr_exp_nmm, binder_exp_nmm = prop_n3_NMM(l, g, P, T, corr_vec, potential_vec, binder_vec)

    print("V (NMM) = ", v_exp_nmm)
    print("V (SOS) = ", v_exp_sos)
    print("Correlation (NMM) = ", corr_exp_nmm)
    print("Correlation (SOS) = ", corr_exp_sos)
    print("Binder (NMM) = ", binder_exp_nmm)
    print("Binder (SOS) = ", binder_exp_sos)
    print("Energy (SOS) = ", energy_exp_sos)
    







