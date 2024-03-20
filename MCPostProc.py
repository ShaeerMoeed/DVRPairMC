import numpy as np
import math
import os
from scipy.optimize import curve_fit
from matplotlib import pyplot as plt

# Currently, this only supports single block MC (unclear whether we will need to break it into multiple blocks)

def process_mc_outputs(dt_folder):

    parameter_sets = os.listdir(dt_folder)
    outFile = os.path.join(dt_folder, "Estimator Statistics.csv")
    
    str_data_out = ""
    out_header_2 = "g,T,P,tau,Potential Mean,Potential Error,Correlation Mean,Correlation Error\n"
    for i in range(len(parameter_sets)):
        parameter_file_id = parameter_sets[i]
        if parameter_file_id[-3:] == "csv" or parameter_file_id[-3:] == "png":
            continue
        parameter_list = parameter_file_id.split("_")
        g = float(parameter_list[1])
        T = float(parameter_list[3])
        P = int(parameter_list[5])
        N = int(parameter_list[7])
        l = int(parameter_list[9])
        filepath = os.path.join(dt_folder, parameter_file_id + "/Block_1/MC Step Outputs.csv")

        with open(filepath, 'r') as f:
            line_count = 0
            for line in f.readlines():
                data = line.strip().split(",")
                line_count += 1
                if line_count == 1:
                    header = data
                    sim_steps = int(header[-1].strip().split("=")[-1].strip())
                    potential_array = np.zeros(sim_steps)
                    corr_array = np.zeros(sim_steps)
                if line_count > 2:
                    potential_array[line_count - 3] = float(data[1])
                    corr_array[line_count - 3] = float(data[2])

        potential_mean, potential_se = calculateError_byBinning(potential_array)
        corr_mean, corr_se = calculateError_byBinning(corr_array)
        tau = 1.0/(float(T) * float(P))
        str_data_out += str(g) + "," + str(T) + "," + str(P) + "," + str(tau) + "," + str(potential_mean) + \
        "," + str(potential_se) + "," + str(corr_mean) + "," + str(corr_se) + "\n"

    out_header_1 = "N = {}, l = {}, Steps = {}\n".format(N, l, sim_steps)
    out_header = out_header_1 + out_header_2
    with open(outFile, 'w') as f:
        f.write(out_header)
        f.write(str_data_out)

    return 0

def process_estimator_outputs(dt_folder):

    inFile = os.path.join(dt_folder, "Estimator Statistics.csv")
    g_T_pair_list = []
    potential_data = []
    corr_data = []
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
                potential_data_g_T = (float(data[4]), float(data[5]))
                corr_data_g_T = (float(data[6]), float(data[7]))
                if g_T_pair not in g_T_pair_list:
                    g_T_pair_list.append(g_T_pair)
                    P_data.append([P_data_g_T])
                    potential_data.append([potential_data_g_T])
                    corr_data.append([corr_data_g_T])
                else: 
                    for i in range(len(g_T_pair_list)):
                        if g_T_pair_list[i] == g_T_pair:
                            P_data[i].append(P_data_g_T)
                            potential_data[i].append(potential_data_g_T)
                            corr_data[i].append(corr_data_g_T)
                            break

    outfile = os.path.join(dt_folder, "Parameter Sweep.csv")
    with open(outfile, 'w') as f:
        f.write(out_header_line)
        f.write("g,T,Potential Mean,Potential Error,Correlation Mean,Correlation Error\n")
        for i in range(len(g_T_pair_list)):
            g_T_pair = g_T_pair_list[i]
            g = g_T_pair[0]
            T = g_T_pair[1]
            P_data_g_T = P_data[i]
            potential_data_g_T = potential_data[i]
            corr_data_g_T = corr_data[i]
            tau_list = list(map(lambda x :float(x[1]), P_data_g_T))
            potential_mean_list = list(map(lambda x :float(x[0]), potential_data_g_T))
            potential_err_list = list(map(lambda x :float(x[1]), potential_data_g_T))
            corr_mean_list = list(map(lambda x :float(x[0]), corr_data_g_T))
            corr_err_list = list(map(lambda x :float(x[1]), corr_data_g_T))
            filename_potential = 'Potential Energy Fit (g = {}, T = {}).png'.format(g, T)
            filename_correlation = 'Correlation Fit (g = {}, T = {}).png'.format(g, T)
            potential_fit_results = extrapolate_results(tau_list, potential_mean_list, func_pair, 100, potential_err_list, 'V', filename_potential, os.path.join(dt_folder, filename_potential))
            corr_fit_results = extrapolate_results(tau_list, corr_mean_list, func_pair, 100, corr_err_list, 'eiej', filename_correlation, os.path.join(dt_folder, filename_correlation))
            f.write(str(g) + "," + str(T) + "," + str(potential_fit_results[1][-1]) + "," + str(potential_fit_results[2]) + "," + str(corr_fit_results[1][-1]) + "," + str(corr_fit_results[2]) + "\n")

    return 0

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
                line_data = (float(data[2]), float(data[3]), float(data[4]), float(data[5]))
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
    
    for i in range(len(T_list)):
        temperature=T_list[i]
        g_T_list = g_list[i]
        data_T_list = data_list[i]
        potential_mean_list = list(map(lambda x :float(x[0]), data_T_list))
        potential_error_list = list(map(lambda x :float(x[1]), data_T_list))
        corr_mean_list = list(map(lambda x :float(x[2]), data_T_list))
        corr_error_list = list(map(lambda x :float(x[3]), data_T_list))
        g_T_list, potential_mean_list, potential_error_list, corr_mean_list, corr_error_list = \
            zip(*sorted(zip(g_T_list, potential_mean_list, potential_error_list, corr_mean_list, corr_error_list)))
        
        plt.figure()
        plt.plot(g_T_list, potential_mean_list)
        plt.scatter(g_T_list, potential_mean_list)
        plt.errorbar(g_T_list, potential_mean_list, potential_error_list,capsize=5,fmt="none")
        plt.title("Effect Of Coupling Strength On Potential (T = {})".format(temperature))
        plt.ylabel("V")
        plt.xlabel("g")
        potentialPlotName = "Potential (T = {}).png".format(temperature)
        potentialPlotPath = os.path.join(dt_folder, potentialPlotName)
        plt.savefig(potentialPlotPath)

        plt.figure()
        plt.plot(g_T_list, corr_mean_list)
        plt.scatter(g_T_list, corr_mean_list)
        plt.errorbar(g_T_list,corr_mean_list,corr_error_list,capsize=5,fmt="none")
        plt.title("Effect Of Coupling Strength On Correlation (T = {})".format(temperature))
        plt.ylabel("eiej")
        plt.xlabel("g")
        corrPlotName = "EiEj (T = {}).png".format(temperature)
        corrPlotPath = os.path.join(dt_folder, corrPlotName)
        plt.savefig(corrPlotPath)
    
    return 0


def extrapolate_results(x_pts, y_pts, func, num_pts_out, error, ylabel, title, filepath):

    if len(y_pts) > 2:

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
        plt.legend()
        plt.xlabel('Tau')
        plt.ylabel(ylabel)
        plt.title(title)
        plt.savefig(filepath)
        plt.close()

    else:
        fit_x = [x_pts[-1]]
        fit_y = [y_pts[-1]]
        err = error[-1]

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
    standardError = maxError_byBinning(arr[trunc:], workingNdim-6)
    return mean, standardError

def func_pair(x, a, b):
    return a*x*x*x + b

if __name__=="__main__":
    
    dt_folder = "/Users/shaeermoeed/Github/DVRPairMC/Results/20_03_2024_10_02_53"
    process_mc_outputs(dt_folder)
    process_estimator_outputs(dt_folder)
    process_parameter_sweeps(dt_folder)
