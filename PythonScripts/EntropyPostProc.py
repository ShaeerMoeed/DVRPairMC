import numpy as np
import math
import os
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt

# Currently, this only supports single block MC (unclear whether we will need to break it into multiple blocks)

def process_mc_file(file):

    skip_rows = 0
    with open(file, 'r') as f:
        line_count = 0
        for line in f.readlines():
            data = line.strip().split(",")
            line_count += 1
            if line_count == 1:
                header = data
                sim_steps = int(header[-1].strip().split("=")[-1].strip())
                numerator_array = np.zeros(sim_steps)
                denominator_array = np.zeros(sim_steps)
                step_number_array = np.zeros(sim_steps)
            if line_count > 2:
                if ((float(data[1]) == 0.0) or (float(data[2]) == 0.0)):
                    skip_rows += 1
                numerator_array[line_count - 3] = float(data[1])
                denominator_array[line_count - 3] = float(data[2])
                step_number_array[line_count - 3] = float(data[0])

    return numerator_array, denominator_array, step_number_array, skip_rows

# When using this function, make sure there is a file for start position all zero and all pi and no random start position results in the folder
def process_mc_outputs_ordered_starts(dt_folder):

    parameter_sets = os.listdir(dt_folder)
    outFile = os.path.join(dt_folder, "Estimator Statistics.csv")
    
    str_data_out = ""
    out_header_2 = "g, T, L, P, tau, Purity, Purity Error\n"
    for i in range(len(parameter_sets)):
        parameter_file_id = parameter_sets[i]
        if parameter_file_id[-3:] == "csv" or parameter_file_id[-3:] == "png" or parameter_file_id[-3:] == "txt":
            continue
        parameter_list = parameter_file_id.split("_")
        if parameter_list[0] == '.DS':
            continue
        g = float(parameter_list[1])
        T = float(parameter_list[3])
        P = int(parameter_list[5])
        N_A = int(parameter_list[7])
        N_B = int(parameter_list[9])
        l = int(parameter_list[11])
        start_config = parameter_list[13]
        exchange_type = parameter_list[16]
        exchange_length = parameter_list[19]
        if start_config != "Zero":
            continue
        if exchange_type != "SwappedUnswapped":
            continue

        filepath_start_zero = os.path.join(dt_folder, parameter_file_id + "/MC Step Outputs.csv")
        num, denom, step_num, skip_rows_1 = process_mc_file(filepath_start_zero)
        h_0 = np.divide(num, step_num)
        h_0 = np.divide(h_0, 1-h_0)

        parameter_file_id_start_pi = "g_{}_T_{}_P_{}_NA_{}_NB_{}_l_{}_Start_{}_exchange_type_{}_exchange_length_{}".format(g, T, P, N_A, N_B, l, "Pi", exchange_type, exchange_length)
        filepath_start_pi = os.path.join(dt_folder, parameter_file_id_start_pi + "/MC Step Outputs.csv")
        num, denom, step_num, skip_rows_2 = process_mc_file(filepath_start_pi)
        
        h_pi = np.divide(num, step_num)
        h_pi = np.divide(h_pi, 1.0-h_pi)
        skip_rows = max((skip_rows_1, skip_rows_2, 1))

        h_0_list = [h_0]
        h_pi_list = [h_pi]

        N = N_B
        N_A_list = [0]
        for j in range(int(N/4)-1):

            N_A += 2
            N_B -= 2

            N_A_list.append(N_A)

            parameter_file_id_start_zero = "g_{}_T_{}_P_{}_NA_{}_NB_{}_l_{}_Start_{}_exchange_type_{}_exchange_length_{}".format(g, T, P, N_A, N_B, l, "Zero", "SwappedSectors", exchange_length)
            filepath_start_0 = os.path.join(dt_folder, parameter_file_id_start_zero + "/MC Step Outputs.csv")
            num, denom, step_num, skip_rows_1 = process_mc_file(filepath_start_0)
            
            h_0 = np.divide(num, step_num)
            h_0 = np.divide(h_0, 1-h_0)

            skip_rows = max((skip_rows, skip_rows_1, 1))

            parameter_file_id_start_pi = "g_{}_T_{}_P_{}_NA_{}_NB_{}_l_{}_Start_{}_exchange_type_{}_exchange_length_{}".format(g, T, P, N_A, N_B, l, "Pi", "SwappedSectors", exchange_length)
            filepath_start_pi = os.path.join(dt_folder, parameter_file_id_start_pi + "/MC Step Outputs.csv")
            num, denom, step_num, skip_rows_1 = process_mc_file(filepath_start_pi)

            h_pi = np.divide(num, step_num)
            h_pi = np.divide(h_pi, 1.0-h_pi)
            skip_rows = max((skip_rows, skip_rows_1, 1))

            h_0_list.append(np.multiply(h_0_list[-1], h_0))
            h_pi_list.append(np.multiply(h_pi_list[-1], h_pi))
            #h_0_list.append(h_0)
            #h_pi_list.append(h_pi)

        purity_mean_list = []
        purity_error_list = []
        for j in range(len(h_0_list)):

            Z_0 = h_0_list[j][skip_rows:]
            Z_pi = h_pi_list[j][skip_rows:]
            
            samples_1, num_blocks_1 = jack_knife_samples(100, Z_0)
            samples_2, num_blocks_2 = jack_knife_samples(100, Z_pi)
        
            num_blocks = num_blocks_1 + num_blocks_2
            samples = np.concatenate((samples_1, samples_2), axis=0)

            # TODO: Try just computing the purity here so that we can fit purity and not S2 by changing func_entropy to func_purity
            purity_mean, purity_err = calculate_jk_estimates(samples, func_identity, num_blocks)
            #samples = np.concatenate((Z_0, Z_pi), axis=0)
            #purity_mean, purity_err = calculateError_byBinning(samples[10000:])
            #purity_mean = func_purity(sample_mean)
            #purity_err = np.abs(sample_err/(1.0 - (sample_mean**2)))

            purity_mean_list.append(purity_mean)
            purity_error_list.append(purity_err)

            if T == 0.0:
                T = "0.05"
            print(T)
            tau = 1.0/(float(T) * float(P))
            str_data_out += str(g) + "," + str(T) + "," + str(N_A_list[j]) + "," + str(P) + "," + str(tau) + "," + str(purity_mean) + "," + str(purity_err) + "\n"

    out_header_1 = "N = {}, l = {}\n".format(N, l)
    out_header = out_header_1 + out_header_2
    with open(outFile, 'w') as f:
        f.write(out_header)
        f.write(str_data_out)

    return 0

def process_mc_outputs_random_start(dt_folder):

    parameter_sets = os.listdir(dt_folder)
    outFile = os.path.join(dt_folder, "Estimator Statistics.csv")
    
    str_data_out = ""
    out_header_2 = "g, T, L, P, tau, Purity, Purity Error\n"
    for i in range(len(parameter_sets)):
        parameter_file_id = parameter_sets[i]
        if parameter_file_id[-3:] == "csv" or parameter_file_id[-3:] == "png" or parameter_file_id[-3:] == "txt":
            continue
        parameter_list = parameter_file_id.split("_")
        if parameter_list[0] == '.DS':
            continue
        g = float(parameter_list[1])
        T = float(parameter_list[3])
        P = int(parameter_list[5])
        N_A = int(parameter_list[7])
        N_B = int(parameter_list[9])
        l = int(parameter_list[11])
        start_config = parameter_list[13]
        exchange_type = parameter_list[16]
        exchange_length = parameter_list[19]
        if start_config != "Random":
            continue
        if exchange_type != "SwappedUnswapped":
            continue

        filepath_start_zero = os.path.join(dt_folder, parameter_file_id + "/MC Step Outputs.csv")
        num, denom, step_num, skip_rows_1 = process_mc_file(filepath_start_zero)
        h_0 = np.divide(num, step_num)
        h_0 = np.divide(h_0, 1-h_0)

        skip_rows = skip_rows_1

        h_0_list = [h_0]

        N = N_B
        N_A_list = [0]
        for j in range(int(N/4)-1):

            N_A += 2
            N_B -= 2

            N_A_list.append(N_A)

            parameter_file_id_start_zero = "g_{}_T_{}_P_{}_NA_{}_NB_{}_l_{}_Start_{}_exchange_type_{}_exchange_length_{}".format(g, T, P, N_A, N_B, l, "Random", "SwappedSectors", exchange_length)
            filepath_start_0 = os.path.join(dt_folder, parameter_file_id_start_zero + "/MC Step Outputs.csv")
            num, denom, step_num, skip_rows_1 = process_mc_file(filepath_start_0)
            
            h_0 = np.divide(num, step_num)
            h_0 = np.divide(h_0, 1-h_0)

            skip_rows = max((skip_rows, skip_rows_1, 1))

            h_0_list.append(np.multiply(h_0_list[-1], h_0))

        purity_mean_list = []
        purity_error_list = []
        for j in range(len(h_0_list)):

            Z_0 = h_0_list[j][skip_rows+5000:]
            
            #samples_1, num_blocks_1 = jack_knife_samples(10000, Z_0)
        
            #num_blocks = num_blocks_1
            #samples = samples_1
            #samples = np.concatenate((Z_0, Z_pi), axis=0)

            # TODO: Try just computing the purity here so that we can fit purity and not S2 by changing func_entropy to func_purity
            #purity_mean, purity_err = calculate_jk_estimates(samples, func_identity, num_blocks)
            samples = Z_0
            purity_mean, purity_err = calculateError_byBinning(samples)
            #purity_mean = func_purity(sample_mean)
            #purity_err = np.abs(sample_err/(1.0 - (sample_mean**2)))

            purity_mean_list.append(purity_mean)
            purity_error_list.append(purity_err)

            tau = 1.0/(float(T) * float(P))
            str_data_out += str(g) + "," + str(T) + "," + str(N_A_list[j]) + "," + str(P) + "," + str(tau) + "," + str(purity_mean) + "," + str(purity_err) + "\n"

    out_header_1 = "N = {}, l = {}\n".format(N, l)
    out_header = out_header_1 + out_header_2
    with open(outFile, 'w') as f:
        f.write(out_header)
        f.write(str_data_out)

    return 0

def process_estimator_outputs(dt_folder):

    inFile = os.path.join(dt_folder, "Estimator Statistics.csv")
    g_T_L_pair_list = []
    purity_data = []
    P_data = []
    with open(inFile, 'r') as f:
        line_count = 0
        for line in f.readlines():
            line_count += 1
            if line_count == 1:
                out_header_line = line
            if line_count > 2: 
                data = line.strip().split(",")
                g_T_L_pair = (float(data[0]), float(data[1]), float(data[2]))
                P_data_g_T_L = (float(data[2]), float(data[3]))
                purity_data_g_T_L = (float(data[5]), float(data[6]))
                if g_T_L_pair not in g_T_L_pair_list:
                    g_T_L_pair_list.append(g_T_L_pair)
                    P_data.append([P_data_g_T_L])
                    purity_data.append([purity_data_g_T_L])
                else: 
                    for i in range(len(g_T_L_pair_list)):
                        if g_T_L_pair_list[i] == g_T_L_pair:
                            P_data[i].append(P_data_g_T_L)
                            purity_data[i].append(purity_data_g_T_L)
                            break

    outfile = os.path.join(dt_folder, "Parameter Sweep.csv")
    with open(outfile, 'w') as f:
        f.write(out_header_line)
        f.write("g,T,L,S2,S2 Error\n")
        for i in range(len(g_T_L_pair_list)):
            g_T_L_pair = g_T_L_pair_list[i]
            g = g_T_L_pair[0]
            T = g_T_L_pair[1]
            L = g_T_L_pair[2]
            P_data_g_T_L = P_data[i]
            purity_data_g_T_L = purity_data[i]
            tau_list = list(map(lambda x :float(x[1]), P_data_g_T_L))
            purity_mean_list = list(map(lambda x :float(x[0]), purity_data_g_T_L))
            purity_err_list = list(map(lambda x :float(x[1]), purity_data_g_T_L))
            filename_purity = 'Renyi Entropy Fit (g = {}, T = {}, L = {}).png'.format(g, T, L)
            S2_fit_results = extrapolate_results(tau_list, purity_mean_list, func_pair, 100, purity_err_list, 
                'S2', filename_purity, os.path.join(dt_folder, filename_purity), use_pts=True, exact_val=None, 
                use_last_pts=False)
            f.write(str(g) + "," + str(T) + "," + str(L) + "," + str(S2_fit_results[1][-1]) + "," + \
                str(S2_fit_results[2]) + "\n")

    return 0

# TODO: Need to update for entropy
def process_parameter_sweeps(dt_folder, process_dmrg=False):

    inFile = os.path.join(dt_folder, "Parameter Sweep.csv")
    T_list = []
    g_list = []
    L_list = []
    data_list = []
    with open(inFile, 'r') as f:
        line_count = 0
        for line in f.readlines():
            line_count += 1
            if line_count == 1:
                out_header_line = line
            if line_count > 2: 
                data = line.strip().split(",")
                T = float(data[1])
                g = float(data[0])
                L = float(data[2])
                line_data = (float(data[3]), float(data[4]))
                if T not in T_list:
                    T_list.append(T)
                    g_list.append([g])
                    L_list.append([[L]])
                    data_list.append([[line_data]])
                else:
                    for i in range(len(T_list)):
                        if T_list[i] == T:
                            if g not in g_list[i]:
                                g_list[i].append(g)
                                L_list[i].append([L])
                                data_list[i].append([line_data])
                            else:
                                for j in range(len(g_list[i])):
                                    if g_list[i][j] == g:
                                        L_list[i][j].append(L)
                                        data_list[i][j].append(line_data)
                                        break
                            break

    for i in range(len(T_list)):
        temperature=T_list[i]
        g_T_list = g_list[i]
        S2_g_means = []
        S2_g_err = []
        for j in range(len(g_T_list)):
            g = g_T_list[j]
            L_g_T_list = L_list[i][j]
            data_T_g_list = data_list[i][j]
            S2_mean_list = list(map(lambda x :float(x[0]), data_T_g_list))
            S2_err_list = list(map(lambda x :float(x[1]), data_T_g_list))
            L_g_T_list, S2_mean_list, S2_err_list = zip(*sorted(zip(L_g_T_list, S2_mean_list, S2_err_list)))

            S2_g_means.append(S2_mean_list[-1])
            S2_g_err.append(S2_err_list[-1])

            plt.figure()
            plt.scatter(L_g_T_list, S2_mean_list, label="g={}".format(g))
            plt.errorbar(L_g_T_list, S2_mean_list, S2_err_list, capsize=5, fmt='None', color='black')
            plt.title("S2 (g={}, T={})".format(g, temperature))
            plt.ylabel("S2")
            plt.xlabel("L")
            plt.legend()
            S2PlotName = "S2 (g = {}, T = {}).png".format(g, temperature)
            S2PlotPath = os.path.join(dt_folder, S2PlotName)
            plt.savefig(S2PlotPath)

        plt.figure()
        #plt.plot(g_T_list, S2_mean_list)
        plt.scatter(g_T_list, S2_g_means, label="MC", color='C0')
        #plt.plot(g_dmrg, S2_dmrg, label="DMRG", color="C3")
        plt.errorbar(g_T_list, S2_g_means, S2_g_err,capsize=5,fmt="none", color='black')
        plt.title("S2 (T = {})".format(temperature))
        plt.ylabel("S2")
        plt.xlabel("g")
        plt.legend()
        S2PlotName = "S2 (T = {}).png".format(temperature)
        S2PlotPath = os.path.join(dt_folder, S2PlotName)
        plt.savefig(S2PlotPath)
    
    return 0


def extrapolate_results(x_pts, y_pts, func, num_pts_out, error, ylabel, title, filepath, use_pts=False, exact_val = None, use_last_pts=False):

    if len(y_pts) > 2:

        if use_last_pts:
            x_pts_sorted, y_pts_sorted = zip(*sorted(zip(x_pts, y_pts)))
            x_pts_sorted, error_sorted = zip(*sorted(zip(x_pts, error)))
            x_pts = x_pts_sorted[:3]
            y_pts = y_pts_sorted[:3]
            error = error_sorted[:3]

        coefficients, cov = curve_fit(func, x_pts, y_pts, sigma=error)
        err = np.sqrt(np.diag(cov))[1]

        fit_x = []
        fit_y = []
        for i in range(num_pts_out):
            x = max(x_pts) * (1.0 - (float(i)/float(num_pts_out-1.0)))
            y = func(x, coefficients[0], coefficients[1])
            fit_x.append(x)
            fit_y.append(y)

        plt.figure()
        plt.plot(fit_x, fit_y, 'b', label='Fit')
        plt.scatter(x_pts, y_pts, color="red", label='MC Data')
        plt.errorbar(x_pts, y_pts, error, ecolor='black', fmt='none',capsize=5)
        plt.scatter([0.0], fit_y[-1], marker="*", color="red", label='Extrapolated')
        plt.errorbar([0.0], fit_y[-1], [err], ecolor='magenta', fmt='none',capsize=5)
        if (exact_val != None):
             plt.plot(fit_x, [exact_val] * len(fit_x), 'red', label='Exact')
        plt.legend()
        plt.xlabel('Tau')
        plt.ylabel(ylabel)
        plt.title(title)
        plt.savefig(filepath)
        plt.close()

        err /= fit_y[-1]
        fit_y = -np.log(fit_y)

    else:
        x_pts_sorted, y_pts_sorted = zip(*sorted(zip(x_pts, y_pts)))
        x_pts_sorted, error_sorted = zip(*sorted(zip(x_pts, error)))
        fit_x = [x_pts_sorted[0]]
        fit_y = [-np.log(y_pts_sorted[0])]
        err = error_sorted[0]/y_pts_sorted[0]

    if use_pts:
        x_pts_sorted, y_pts_sorted = zip(*sorted(zip(x_pts, y_pts)))
        x_pts_sorted, error_sorted = zip(*sorted(zip(x_pts, error)))
        fit_x = [x_pts_sorted[0]]
        fit_y = [-np.log(y_pts_sorted[0])]
        err = error_sorted[0]/y_pts_sorted[0]

    return fit_x, fit_y, err

def errorpropagation(data):
    ndim   = len(data)
    error = np.std(data,ddof=0)/np.sqrt(ndim)
    return error

def maxError_byBinning(data, workingNdim):
    if(workingNdim<=1):
        raise Exception('Not enough points MC steps were used for the binning method, please increase the number of MC steps')
    error = np.zeros(workingNdim)
    i = 0
    error[0] = errorpropagation(data)

    for i in range(1,workingNdim):
        ndim = int(len(data)/2)
        data1 = np.zeros(ndim)

        for j in range(ndim):
            data1[j] = 0.5*(data[2*j]+data[2*j+1])
        data = data1
        error[i] = errorpropagation(data)
    return np.max(error)

def calculateError_byBinning(arr):
    # Finding the average and standard error using the binning method
    # This method requires 2^n data points, this truncates the data to fit this
    workingNdim  = int(math.log(len(arr))/math.log(2))
    trunc = int(len(arr)-2**workingNdim)
    mean = np.mean(arr[trunc:])
    mean = arr[-1]
    standardError = maxError_byBinning(arr[trunc:], workingNdim-6)
    return mean, standardError

def jack_knife_samples(decorrelation_time, h_vec):

    block_size = decorrelation_time
    total_pts = np.size(h_vec)

    num_pts_to_use = int(np.floor(total_pts/block_size) * block_size)
    num_blocks = int(num_pts_to_use/block_size)

    h_vec = h_vec[(total_pts - num_pts_to_use):]
    block_data = h_vec.reshape((num_blocks, block_size))
    print(np.shape(block_data))

    samples = np.zeros((num_blocks))
    for i in range(num_blocks):
        samples[i] = block_data[i,block_size-1]

    print(len(samples))

    return samples, num_blocks

def calculate_jk_estimates(samples, func, num_blocks):

    j_mean_array = np.zeros((num_blocks))
    func_array = np.zeros((num_blocks))
    full_mean = np.mean(samples)
    print(full_mean)
    for i in range(num_blocks):
        j_mean_array[i] = ((full_mean * num_blocks) - samples[i])/(num_blocks-1)
        func_array[i] = func(j_mean_array[i])
    
    func_mean_jk = np.mean(func_array)
    func_std_dev_jk = np.std(func_array)
    func_std_error_jk = np.sqrt((num_blocks-1)) * func_std_dev_jk

    func_mean = num_blocks * func(full_mean) - (num_blocks - 1)*func_mean_jk

    return func_mean, func_std_error_jk

def func_entropy(x):

    return -np.log(x/(1.0 - x))

def func_purity(x):

    return x/(1.0 - x)

def func_identity(x):

    return x

def func_pair(x, a, b):
    return a + b*x*x*x

if __name__=="__main__":
    
    #mc_output_folder = "/Users/shaeermoeed/Github/DVRPairMC/Results/Entropy/Test/20_06_2024_17_57_02"
    mc_output_folder = "/Users/shaeermoeed/Github/DVRPairMC/Results/Entropy/Test/25_11_2025_16_42_08"
    process_mc_outputs_ordered_starts(mc_output_folder)

    #mc_output_folder = "/Users/shaeermoeed/Github/DVRPairMC/Results/Entropy/Test/25_11_2025_06_18_59"
    #process_mc_outputs_random_start(mc_output_folder)
    
    process_estimator_outputs(mc_output_folder)
    process_parameter_sweeps(mc_output_folder)

    '''
    file_1 = "/Users/shaeermoeed/Github/DVRPairMC/Results/Entropy/Test/20_06_2024_17_57_02/Parameter Sweep Binning.csv"
    file_2 = "/Users/shaeermoeed/Github/DVRPairMC/Results/Entropy/Test/20_06_2024_17_57_02/Parameter Sweep.csv"
    file_3 = "/Users/shaeermoeed/Github/DVRPairMC/Results/Entropy/Test/20_06_2024_17_57_02/entropy_Renyi.txt"
    g1, T1, S2_1, S2_err_1 = np.loadtxt(file_1, unpack=True, delimiter=",", skiprows=2)
    g2, T2, S2_2, S2_err_2 = np.loadtxt(file_2, unpack=True, delimiter=",", skiprows=2)
    #dmrg_g, dmrg_S2 = np.loadtxt(os.path.join(mc_output_folder, "entropy_Renyi.txt"), skiprows=3, unpack=True)

    log2 = [np.log(2.0)]*len(g1)

    plt.figure()
    plt.rcParams['mathtext.fontset']='stix'
    plt.scatter(g1, S2_1, label="MC", color='C0')
    plt.errorbar(g1, S2_1, S2_err_2, fmt='None', capsize=5)
    #plt.scatter(g2, S2_2, label="Jack-Knife", color='C0')
    #plt.errorbar(g2, S2_2, S2_err_2, fmt='None', capsize=5, color="black")
    #plt.plot(dmrg_g[:14], dmrg_S2[:14], label="DMRG", color='C3')
    plt.plot(g1, log2, label=r"ln(2)", color='black')
    plt.ylabel(r"$S_2$", fontsize=20)
    plt.xlabel(r"$g$", fontsize=20)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.legend()
    plt.savefig("Entropy.png", bbox_inches='tight', dpi=1200)
    plt.show()
    '''