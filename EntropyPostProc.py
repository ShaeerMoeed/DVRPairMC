import numpy as np
import math
import os
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt

# Currently, this only supports single block MC (unclear whether we will need to break it into multiple blocks)

# When using this function, make sure there is a file for start position all zero and all pi and no random start position results in the folder
def process_mc_outputs(dt_folder):

    parameter_sets = os.listdir(dt_folder)
    outFile = os.path.join(dt_folder, "Estimator Statistics.csv")
    
    str_data_out = ""
    out_header_2 = "g,T,P,tau, Purity, Purity Error\n"
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
        if start_config != "Zero":
            continue
        filepath_start_zero = os.path.join(dt_folder, parameter_file_id + "/MC Step Outputs.csv")

        skip_rows_1 = 0
        with open(filepath_start_zero, 'r') as f:
            line_count = 0
            for line in f.readlines():
                data = line.strip().split(",")
                line_count += 1
                if line_count == 1:
                    header = data
                    sim_steps = int(header[-1].strip().split("=")[-1].strip())
                    numerator_array_start_zero = np.zeros(sim_steps)
                    denominator_array_start_zero = np.zeros(sim_steps)
                    step_number_array_start_zero = np.zeros(sim_steps)
                if line_count > 2:
                    if ((float(data[1]) == 0.0) or (float(data[2]) == 0.0)):
                        skip_rows_1 += 1
                    step_number_array_start_zero[line_count - 3] = float(data[0])
                    numerator_array_start_zero[line_count - 3] = float(data[1])
                    denominator_array_start_zero[line_count - 3] = float(data[2])

        # TODO: Generalize so that all parameter sets with the same parameters and different start positions are automatically processed in a loop  
        skip_rows_2 = 0
        parameter_file_id_start_pi = "g_{}_T_{}_P_{}_NA_{}_NB_{}_l_{}_Start_{}".format(g, T, P, N_A, N_B, l, "Pi")
        filepath_start_pi = os.path.join(dt_folder, parameter_file_id_start_pi + "/MC Step Outputs.csv")
        with open(filepath_start_pi, 'r') as f:
            line_count = 0
            for line in f.readlines():
                data = line.strip().split(",")
                line_count += 1
                if line_count == 1:
                    header = data
                    sim_steps = int(header[-1].strip().split("=")[-1].strip())
                    numerator_array_start_pi = np.zeros(sim_steps)
                    denominator_array_start_pi = np.zeros(sim_steps)
                    step_number_array_start_pi = np.zeros(sim_steps)
                if line_count > 2:
                    if ((float(data[1]) == 0.0) or (float(data[2]) == 0.0)):
                        skip_rows_2 += 1
                    numerator_array_start_pi[line_count - 3] = float(data[1])
                    denominator_array_start_pi[line_count - 3] = float(data[2])
                    step_number_array_start_pi[line_count - 3] = float(data[0])
                    #numerator_array_start_pi[line_count - 3] = float(data[1])
                    #denominator_array_start_pi[line_count - 3] = float(data[2])

        skip_rows = max((skip_rows_1, skip_rows_2, 1))

        h_1 = np.divide(numerator_array_start_zero, step_number_array_start_zero)
        h_2 = np.divide(numerator_array_start_pi, step_number_array_start_pi)

        h_1 = h_1[skip_rows:]
        h_2 = h_2[skip_rows:]

        samples_1, num_blocks_1 = jack_knife_samples(1, h_1)
        samples_2, num_blocks_2 = jack_knife_samples(1, h_2)

        num_blocks = num_blocks_1 + num_blocks_2
        samples = np.concatenate((samples_1, samples_2), axis=0)
        # TODO: Try just computing the purity here so that we can fit purity and not S2 by changing func_entropy to func_purity
        purity_mean, purity_err = calculate_jk_estimates(samples, func_purity, num_blocks)
        #sample_mean, sample_err = calculateError_byBinning(samples)
        #purity_mean = func_purity(sample_mean)
        #purity_err = np.abs(sample_err/(1.0 - (sample_mean**2)))

        tau = 1.0/(float(T) * float(P))
        str_data_out += str(g) + "," + str(T) + "," + str(P) + "," + str(tau) + "," + str(purity_mean) + "," + str(purity_err) + "\n"

    out_header_1 = "NA = {}, NB = {}, l = {}, Steps = {}\n".format(N_A, N_B, l, sim_steps)
    out_header = out_header_1 + out_header_2
    with open(outFile, 'w') as f:
        f.write(out_header)
        f.write(str_data_out)

    return 0

def process_estimator_outputs(dt_folder):

    inFile = os.path.join(dt_folder, "Estimator Statistics.csv")
    g_T_pair_list = []
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
                g_T_pair = (float(data[0]), float(data[1]))
                P_data_g_T = (float(data[2]), float(data[3]))
                purity_data_g_T = (float(data[4]), float(data[5]))
                if g_T_pair not in g_T_pair_list:
                    g_T_pair_list.append(g_T_pair)
                    P_data.append([P_data_g_T])
                    purity_data.append([purity_data_g_T])
                else: 
                    for i in range(len(g_T_pair_list)):
                        if g_T_pair_list[i] == g_T_pair:
                            P_data[i].append(P_data_g_T)
                            purity_data[i].append(purity_data_g_T)
                            break

    outfile = os.path.join(dt_folder, "Parameter Sweep.csv")
    with open(outfile, 'w') as f:
        f.write(out_header_line)
        f.write("g,T,S2, S2 Error\n")
        for i in range(len(g_T_pair_list)):
            g_T_pair = g_T_pair_list[i]
            g = g_T_pair[0]
            T = g_T_pair[1]
            P_data_g_T = P_data[i]
            purity_data_g_T = purity_data[i]
            tau_list = list(map(lambda x :float(x[1]), P_data_g_T))
            purity_mean_list = list(map(lambda x :float(x[0]), purity_data_g_T))
            purity_err_list = list(map(lambda x :float(x[1]), purity_data_g_T))
            filename_purity = 'Renyi Entropy Fit (g = {}, T = {}).png'.format(g, T)
            S2_fit_results = extrapolate_results(tau_list, purity_mean_list, func_pair, 100, purity_err_list, 
                'S2', filename_purity, os.path.join(dt_folder, filename_purity), use_pts=True, exact_val=None, use_last_pts=False)
            f.write(str(g) + "," + str(T) + "," + str(S2_fit_results[1][-1]) + "," + \
                str(S2_fit_results[2]) + "\n")

    return 0

# TODO: Need to update for entropy
def process_parameter_sweeps(dt_folder):

    inFile = os.path.join(dt_folder, "Parameter Sweep.csv")
    T_list = []
    g_list = []
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
                line_data = (float(data[2]), float(data[3]))
                if T not in T_list:
                    T_list.append(T)
                    g_list.append([g])
                    data_list.append([line_data])
                else:
                    for i in range(len(T_list)):
                        if T_list[i] == T:
                            g_list[i].append(g)
                            data_list[i].append(line_data)
                            break
    
    g_dmrg, S2_dmrg = np.loadtxt(os.path.join(dt_folder, "entropy_Renyi.txt"), skiprows=3, unpack=True)

    for i in range(len(T_list)):
        temperature=T_list[i]
        g_T_list = g_list[i]
        data_T_list = data_list[i]
        S2_mean_list = list(map(lambda x :float(x[0]), data_T_list))
        S2_err_list = list(map(lambda x :float(x[1]), data_T_list))
        g_T_list, S2_mean_list, S2_err_list = zip(*sorted(zip(g_T_list, S2_mean_list, S2_err_list)))

        if max(g_dmrg) > max(g_T_list):
            for i in range(len(g_dmrg)):
                if g_dmrg[i] > max(g_T_list):
                    g_dmrg = g_dmrg[:i]
                    S2_dmrg = S2_dmrg[:i]
                    break

        plt.figure()
        #plt.plot(g_T_list, S2_mean_list)
        plt.scatter(g_T_list, S2_mean_list, label="MC", color='C0')
        plt.plot(g_dmrg, S2_dmrg, label="DMRG", color="C3")
        plt.errorbar(g_T_list, S2_mean_list, S2_err_list,capsize=5,fmt="none", color='black')
        plt.title("Effect Of Coupling Strength On S2 (T = {})".format(temperature))
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

def func_pair(x, a, b):
    return a + b*x*x*x

if __name__=="__main__":
    
    dt_folder = "/Users/shaeermoeed/Github/DVRPairMC/Results/Entropy/Test/20_06_2024_17_57_02"
    #process_mc_outputs(dt_folder)
    
    process_estimator_outputs(dt_folder)
    process_parameter_sweeps(dt_folder)
    
    file_1 = "/Users/shaeermoeed/Github/DVRPairMC/Results/Entropy/Test/20_06_2024_17_57_02/Parameter Sweep Binning.csv"
    file_2 = "/Users/shaeermoeed/Github/DVRPairMC/Results/Entropy/Test/20_06_2024_17_57_02/Parameter Sweep.csv"
    file_3 = "/Users/shaeermoeed/Github/DVRPairMC/Results/Entropy/Test/20_06_2024_17_57_02/entropy_Renyi.txt"
    g1, T1, S2_1, S2_err_1 = np.loadtxt(file_1, unpack=True, delimiter=",", skiprows=2)
    g2, T2, S2_2, S2_err_2 = np.loadtxt(file_2, unpack=True, delimiter=",", skiprows=2)
    g_dmrg, S2_dmrg = np.loadtxt(os.path.join(dt_folder, "entropy_Renyi.txt"), skiprows=3, unpack=True)

    plt.figure()
    plt.scatter(g1, S2_1, label="MC", color='C0')
    plt.errorbar(g1, S2_1, S2_err_2, fmt='None', capsize=5)
    #plt.scatter(g2, S2_2, label="Jack-Knife", color='C0')
    #plt.errorbar(g2, S2_2, S2_err_2, fmt='None', capsize=5, color="black")
    plt.plot(g_dmrg[:14], S2_dmrg[:14], label="DMRG", color='C3')
    plt.title("Renyi Entropy Convergence")
    plt.ylabel("S2")
    plt.xlabel("g")
    plt.legend()
    plt.savefig("Entropy.png", bbox_inches='tight')
    plt.show()
    


