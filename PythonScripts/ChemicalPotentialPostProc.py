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
    out_header_2 = "N,T,P,tau,Potential Mean,Potential Error,Correlation Mean,Correlation Error, Binder Mean, Binder Error, Energy Mean, Energy Error\n"
    for i in range(len(parameter_sets)):
        parameter_file_id = parameter_sets[i]
        if parameter_file_id[-3:] == "csv" or parameter_file_id[-3:] == "png" or parameter_file_id[-3:] == "txt":
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
                    binder_array = np.zeros(sim_steps)
                    energy_array = np.zeros(sim_steps)
                if line_count > 2:
                    potential_array[line_count - 3] = float(data[1])
                    corr_array[line_count - 3] = float(data[2])
                    binder_array[line_count - 3] = float(data[3])
                    energy_array[line_count - 3] = float(data[4])

        potential_mean, potential_se = calculateError_byBinning(potential_array)

        corr_mean, corr_se = calculateError_byBinning(corr_array)
        corr_jk_samples, num_blocks = jack_knife_samples(5000, corr_array, 20)
        corr_mean, corr_se = calculate_corr_jk_estimates(corr_jk_samples, num_blocks, N) 

        #binder_mean, binder_se = calculateError_byBinning(binder_array)
        binder_denom_array = np.copy(binder_array)
        binder_num_array = list(map(lambda x :x**2, binder_array))

        if g <= 0.61:
            binder_num_jk_samples, num_blocks = jack_knife_samples(20000, binder_num_array, 5)
            binder_denom_jk_samples, num_blocks = jack_knife_samples(20000, binder_denom_array, 5)
            binder_mean, binder_error = calculate_binder_jk_estimates(binder_num_jk_samples, binder_denom_jk_samples, num_blocks)
        if g > 0.61:
            binder_num_jk_samples, num_blocks = jack_knife_samples(1, binder_num_array, 5)
            binder_denom_jk_samples, num_blocks = jack_knife_samples(1, binder_denom_array, 5)
            binder_mean, binder_error = calculate_binder_jk_estimates(binder_num_jk_samples, binder_denom_jk_samples, num_blocks) 

        '''
        binder_denom_mean, binder_denom_se = calculateError_byBinning(binder_denom_array)
        binder_num_mean, binder_num_se = calculateError_byBinning(binder_num_array)
        '''

        #binder_mean = 1.0 - binder_num_mean/(3.0 * (binder_denom_mean**2))
        # TODO: Figure out correct statistical error for binder
        #binder_se = binder_denom_se/(N**2)

        #energy_array = energy_array[5000:]
        energy_mean, energy_se = calculateError_byBinning(energy_array[5000:])
        #e0_jk_samples, num_blocks = jack_knife_samples(int(np.shape(energy_array)[0]/10), energy_array, 5)
        #energy_mean, energy_se = calculate_energy_jk_estimates(e0_jk_samples, num_blocks) 
        tau = 1.0/(float(T) * float(P))
        str_data_out += str(N) + "," + str(T) + "," + str(P) + "," + str(tau) + "," + str(potential_mean) + \
        "," + str(potential_se) + "," + str(corr_mean) + "," + str(corr_se) + "," + str(binder_mean) + "," + str(binder_error) + \
        "," + str(energy_mean) + "," + str(energy_se) + "\n"

    out_header_1 = "g = {}, l = {}, Steps = {}\n".format(g, l, sim_steps)
    out_header = out_header_1 + out_header_2
    with open(outFile, 'w') as f:
        f.write(out_header)
        f.write(str_data_out)

    return 0

def process_estimator_outputs(dt_folder):

    inFile = os.path.join(dt_folder, "Estimator Statistics.csv")
    N_T_pair_list = []
    potential_data = []
    corr_data = []
    binder_data = []
    energy_data = []
    P_data = []
    with open(inFile, 'r') as f:
        line_count = 0
        for line in f.readlines():
            line_count += 1
            if line_count == 1:
                out_header_line = line
            if line_count > 2: 
                data = line.strip().split(",")
                N_T_pair = (float(data[0]), float(data[1]))
                P_data_N_T = (float(data[2]), float(data[3]))
                potential_data_N_T = (float(data[4]), float(data[5]))
                corr_data_N_T = (float(data[6]), float(data[7]))
                binder_data_N_T = (float(data[8]), float(data[9]))
                energy_data_N_T = (float(data[10]), float(data[11]))
                if N_T_pair not in N_T_pair_list:
                    N_T_pair_list.append(N_T_pair)
                    P_data.append([P_data_N_T])
                    potential_data.append([potential_data_N_T])
                    corr_data.append([corr_data_N_T])
                    binder_data.append([binder_data_N_T])
                    energy_data.append([energy_data_N_T])
                else: 
                    for i in range(len(N_T_pair_list)):
                        if N_T_pair_list[i] == N_T_pair:
                            P_data[i].append(P_data_N_T)
                            potential_data[i].append(potential_data_N_T)
                            corr_data[i].append(corr_data_N_T)
                            binder_data[i].append(binder_data_N_T)
                            energy_data[i].append(energy_data_N_T)
                            break
    
    g = float(out_header_line.split(",")[0].split("=")[1].strip())
    outfile = os.path.join(dt_folder, "Parameter Sweep.csv")
    with open(outfile, 'w') as f:
        f.write(out_header_line)
        f.write("N,T,Potential Mean,Potential Error,Correlation Mean,Correlation Error, Binder Mean, Binder Error, Energy Mean, Energy Error\n")
        for i in range(len(N_T_pair_list)):
            N_T_pair = N_T_pair_list[i]
            N = N_T_pair[0]
            T = N_T_pair[1]
            P_data_N_T = P_data[i]
            potential_data_N_T = potential_data[i]
            corr_data_N_T = corr_data[i]
            binder_data_N_T = binder_data[i]
            energy_data_N_T = energy_data[i]
            tau_list = list(map(lambda x :float(x[1]), P_data_N_T))
            potential_mean_list = list(map(lambda x :float(x[0]), potential_data_N_T))
            potential_err_list = list(map(lambda x :float(x[1]), potential_data_N_T))
            corr_mean_list = list(map(lambda x :float(x[0]), corr_data_N_T))
            corr_err_list = list(map(lambda x :float(x[1]), corr_data_N_T))
            binder_mean_list = list(map(lambda x :float(x[0]), binder_data_N_T))
            binder_err_list = list(map(lambda x :float(x[1]), binder_data_N_T))
            energy_mean_list = list(map(lambda x :float(x[0]), energy_data_N_T))
            energy_err_list = list(map(lambda x :float(x[1]), energy_data_N_T))
            filename_potential = 'Potential Energy Fit (N = {}, T = {}, g = {}).png'.format(N, T, g)
            filename_correlation = 'Correlation Fit (N = {}, T = {}, g = {}).png'.format(N, T, g)
            filename_binder = 'Binder Ratio Fit (N = {}, T = {}, g = {}).png'.format(N, T, g)
            filename_energy = 'Energy Fit (N = {}, T = {}, g = {}).png'.format(N, T, g)
            potential_fit_results = extrapolate_results(tau_list, potential_mean_list, func_quadratic, 100, potential_err_list, 'V', filename_potential, os.path.join(dt_folder, filename_potential))
            
            if g == 0.5 or g == 0.6:
                binder_fit_results = extrapolate_results(tau_list, binder_mean_list, func_quadratic, 100, binder_err_list, 'Binder Ratio', filename_binder, os.path.join(dt_folder, filename_binder), use_pts=False)
                corr_fit_results = extrapolate_results(tau_list, corr_mean_list, func_linear, 100, corr_err_list, 'eiej', filename_correlation, os.path.join(dt_folder, filename_correlation), use_pts=False)
            else:
                corr_fit_results = extrapolate_results(tau_list, corr_mean_list, func_quadratic, 100, corr_err_list, 'eiej', filename_correlation, os.path.join(dt_folder, filename_correlation), use_pts=True)
                binder_fit_results = extrapolate_results(tau_list, binder_mean_list, func_quadratic, 100, binder_err_list, 'Binder Ratio', filename_binder, os.path.join(dt_folder, filename_binder), use_pts=True)
            '''
            if g == 0.5 or g == 0.6:
                binder_fit_results = extrapolate_results(tau_list, binder_mean_list, func_linear, 100, binder_err_list, 'Binder Ratio', filename_binder, os.path.join(dt_folder, filename_binder), use_pts=False)
                corr_fit_results = extrapolate_results(tau_list, corr_mean_list, func_linear, 100, corr_err_list, 'eiej', filename_correlation, os.path.join(dt_folder, filename_correlation), use_pts=False)
            else:
                corr_fit_results = extrapolate_results(tau_list, corr_mean_list, func_cubic, 100, corr_err_list, 'eiej', filename_correlation, os.path.join(dt_folder, filename_correlation), use_pts=True)
                binder_fit_results = extrapolate_results(tau_list, binder_mean_list, func_cubic, 100, binder_err_list, 'Binder Ratio', filename_binder, os.path.join(dt_folder, filename_binder), use_pts=True)
            '''
            energy_fit_results = extrapolate_results(tau_list, energy_mean_list, func_quadratic, 100, energy_err_list, 'Energy', filename_energy, os.path.join(dt_folder, filename_energy), use_pts=False)
            f.write(str(N) + "," + str(T) + "," + str(potential_fit_results[1][-1]) + "," + \
                    str(potential_fit_results[2]) + "," + str(corr_fit_results[1][-1]) + "," + \
                    str(corr_fit_results[2]) + "," + str(binder_fit_results[1][-1]) + "," + \
                    str(binder_fit_results[2]) + "," + str(energy_fit_results[1][-1]) + "," + \
                    str(energy_fit_results[2]) + "\n")

    return 0

def process_parameter_sweeps(dt_folder, cumulative):
    
    if cumulative: 
        file_name = "Parameter Sweep Cumulative.csv"
    else: 
        file_name = "Parameter Sweep.csv"
    inFile = os.path.join(dt_folder, file_name)
    T_list = []
    N_list = []
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
                N = float(data[0])
                line_data = (float(data[2]), float(data[3]), float(data[4]), float(data[5]), float(data[6]), float(data[7]), float(data[8]), float(data[9]))
                if T not in T_list:
                    T_list.append(T)
                    N_list.append([N])
                    data_list.append([line_data])
                else:
                    for i in range(len(T_list)):
                        if T_list[i] == T:
                            N_list[i].append(N)
                            data_list[i].append(line_data)
                            break
    
    for i in range(len(T_list)):
        temperature=T_list[i]
        N_T_list = N_list[i]
        data_T_list = data_list[i]
        potential_mean_list = list(map(lambda x :float(x[0]), data_T_list))
        potential_error_list = list(map(lambda x :float(x[1]), data_T_list))
        corr_mean_list = list(map(lambda x :float(x[2]), data_T_list))
        corr_error_list = list(map(lambda x :float(x[3]), data_T_list))
        binder_mean_list = list(map(lambda x :float(x[4]), data_T_list))
        binder_error_list = list(map(lambda x :float(x[5]), data_T_list))
        energy_mean_list = list(map(lambda x :float(x[6]), data_T_list))
        energy_error_list = list(map(lambda x :float(x[7]), data_T_list))
        N_T_list, potential_mean_list, potential_error_list, corr_mean_list, corr_error_list, binder_mean_list, binder_error_list, energy_mean_list, energy_error_list = \
            zip(*sorted(zip(N_T_list, potential_mean_list, potential_error_list, corr_mean_list, corr_error_list, binder_mean_list, binder_error_list,
                            energy_mean_list, energy_error_list)))
        
        if cumulative:
            potentialPlotName = "Potential Cumulative (T = {}).png".format(temperature)
            corrPlotName = "EiEj Cumulative (T = {}).png".format(temperature)
            binderPlotName = "Binder Ratio Cumulative (T = {}).png".format(temperature)
            energyPlotName = "Energy Cumulative (T = {}).png".format(temperature)
        else:
            potentialPlotName = "Potential (T = {}).png".format(temperature)
            corrPlotName = "EiEj (T = {}).png".format(temperature)
            binderPlotName = "Binder Ratio (T = {}).png".format(temperature)
            energyPlotName = "Energy (T = {}).png".format(temperature)
        
        plt.figure()
        plt.plot(N_T_list, potential_mean_list)
        plt.scatter(N_T_list, potential_mean_list)
        plt.errorbar(N_T_list, potential_mean_list, potential_error_list,capsize=5,fmt="none")
        plt.title("Effect Of Particle Number On Potential (T = {})".format(temperature))
        plt.ylabel("V")
        plt.xlabel("N")
        potentialPlotPath = os.path.join(dt_folder, potentialPlotName)
        plt.savefig(potentialPlotPath)
        plt.close()

        plt.figure()
        plt.plot(N_T_list, corr_mean_list)
        plt.scatter(N_T_list, corr_mean_list)
        plt.errorbar(N_T_list,corr_mean_list,corr_error_list,capsize=5,fmt="none")
        plt.title("Effect Of Particle Number On Correlation (T = {})".format(temperature))
        plt.ylabel("eiej")
        plt.xlabel("N")
        corrPlotPath = os.path.join(dt_folder, corrPlotName)
        plt.savefig(corrPlotPath)
        plt.close()

        plt.figure()
        plt.plot(N_T_list, binder_mean_list)
        plt.scatter(N_T_list, binder_mean_list)
        plt.errorbar(N_T_list,binder_mean_list,binder_error_list,capsize=5,fmt="none")
        plt.title("Effect Of Particle Number On Binder Ratio (T = {})".format(temperature))
        plt.ylabel("Binder")
        plt.xlabel("N")
        binderPlotPath = os.path.join(dt_folder, binderPlotName)
        plt.savefig(binderPlotPath)
        plt.close()

        plt.figure()
        plt.plot(N_T_list, energy_mean_list)
        plt.scatter(N_T_list, energy_mean_list)
        plt.errorbar(N_T_list,energy_mean_list,energy_error_list,capsize=5,fmt="none")
        plt.title("Effect Of Particle Number On Energy (T = {})".format(temperature))
        plt.ylabel("Energy")
        plt.xlabel("N")
        energyPlotPath = os.path.join(dt_folder, energyPlotName)
        plt.savefig(energyPlotPath)

        plt.close()
    
    return 0


def extrapolate_results(x_pts, y_pts, func, num_pts_out, error, ylabel, title, filepath, use_pts=False):

    if len(y_pts) > 2:

        coefficients, cov = curve_fit(func, x_pts, y_pts, sigma=error)
        err = np.sqrt(np.diag(cov))[1]

        fit_x = []
        fit_y = []
        for i in range(num_pts_out):
            x = max(x_pts) * (1.0 - (float(i)/float(num_pts_out-1.0)))
            if len(coefficients) == 2:
                y = func(x, coefficients[0], coefficients[1])
            elif len(coefficients) == 3:
                y = func(x, coefficients[0], coefficients[1], coefficients[2])
            else:
                raise Exception("Length of coefficients should be either 2 or 3\n")
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
        x_pts_sorted, y_pts_sorted = zip(*sorted(zip(x_pts, y_pts)))
        x_pts_sorted, error_sorted = zip(*sorted(zip(x_pts, error)))
        fit_x = [x_pts_sorted[0]]
        fit_y = [y_pts_sorted[0]]
        err = error_sorted[0]

    if use_pts:
        x_pts_sorted, y_pts_sorted = zip(*sorted(zip(x_pts, y_pts)))
        x_pts_sorted, error_sorted = zip(*sorted(zip(x_pts, error)))
        fit_x = [x_pts_sorted[0]]
        fit_y = [y_pts_sorted[0]]
        err = error_sorted[0]

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

def func_linear(x, a, b):
    return a*x + b

def func_cubic(x, a, b):
    return a*x*x*x + b

def func_quadratic(x,a, b):
    return a*x*x + b

def jack_knife_samples(decorrelation_time, h_vec, block_size=50):

    # block_size-1 is the number of points that will be dropped between successive samples
    h_vec = np.array(h_vec[decorrelation_time:])
    #block_size = decorrelation_time
    total_pts = np.size(h_vec)

    num_pts_to_use = int(np.floor(total_pts/block_size) * block_size)
    num_blocks = int(num_pts_to_use/block_size)

    h_vec = h_vec[(total_pts - num_pts_to_use):]
    block_data = h_vec.reshape((num_blocks, block_size))

    samples = np.zeros((num_blocks))
    for i in range(num_blocks):
        samples[i] = block_data[i,block_size-1]

    return samples, num_blocks

def calculate_binder_jk_estimates(samples_num, samples_denom, num_blocks):

    j1_mean_array = np.zeros((num_blocks))
    j2_mean_array = np.zeros((num_blocks))
    func_array = np.zeros((num_blocks))
    full_mean_num = np.mean(samples_num)
    full_mean_denom = np.mean(samples_denom)

    for i in range(num_blocks):
        j1_mean_array[i] = ((full_mean_num * num_blocks) - samples_num[i])/(num_blocks-1)
        j2_mean_array[i] = ((full_mean_denom * num_blocks) - samples_denom[i])/(num_blocks-1)
        func_array[i] = func_binder_ratio(j1_mean_array[i], j2_mean_array[i])
    
    func_mean_jk = np.mean(func_array)
    func_std_dev_jk = np.std(func_array)
    func_std_error_jk = np.sqrt((num_blocks-1)) * func_std_dev_jk

    func_mean = num_blocks * func_binder_ratio(full_mean_num, full_mean_denom) - (num_blocks - 1)*func_mean_jk

    return func_mean, func_std_error_jk

def calculate_corr_jk_estimates(jk_samples, num_blocks, num_rotors):

    jk_mean_array = np.zeros((num_blocks))
    func_array = np.zeros((num_blocks))
    full_mean = np.mean(jk_samples)

    for i in range(num_blocks):
        jk_mean_array[i] = ((full_mean * num_blocks) - jk_samples[i])/(num_blocks-1)
        func_array[i] = jk_mean_array[i]/num_rotors
    
    func_mean_jk = np.mean(func_array)
    func_std_dev_jk = np.std(func_array)
    func_std_error_jk = np.sqrt((num_blocks-1)) * func_std_dev_jk

    func_mean = num_blocks * (full_mean/num_rotors) - (num_blocks - 1)*func_mean_jk

    return func_mean, func_std_error_jk

def calculate_energy_jk_estimates(jk_samples, num_blocks):

    jk_mean_array = np.zeros((num_blocks))
    func_array = np.zeros((num_blocks))
    full_mean = np.mean(jk_samples)

    for i in range(num_blocks):
        jk_mean_array[i] = ((full_mean * num_blocks) - jk_samples[i])/(num_blocks-1)
        func_array[i] = jk_mean_array[i]
    
    func_mean_jk = np.mean(func_array)
    func_std_dev_jk = np.std(func_array)
    func_std_error_jk = np.sqrt((num_blocks-1)) * func_std_dev_jk

    func_mean = num_blocks * (full_mean) - (num_blocks - 1)*func_mean_jk

    return func_mean, func_std_error_jk

def func_binder_ratio(num, denom):

    return 1.0 - (num/(3.0*(denom**2)))

def combine_parameter_sweeps(dt_folder1, dt_folder2, filename_1, filename_2):

    inFile1 = os.path.join(dt_folder1, filename_1)
    inFile2 = os.path.join(dt_folder2, filename_2)

    data_1 = np.loadtxt(inFile1, delimiter=",", skiprows=2)
    data_2 = np.loadtxt(inFile2, delimiter=",", skiprows=2)

    data = np.zeros(np.shape(data_1))
    for i in range(np.shape(data_1)[0]):
        data[i, 0] = data_1[i, 0]
        data[i, 1] = data_1[i, 1]

    for i in range(4):
        for j in range(np.shape(data_1)[0]):
            for k in range(np.shape(data_2)[0]):
                if ((data_1[j, 0] == data_2[k,0]) and (data_1[j,1] == data_2[k,1])):
                    data_1[j, 2+2*i] = 0.5*(data_1[j, 2+2*i] + data_2[k, 2+2*i])
                    data_1[j, 2+2*i + 1] = 0.5*np.sqrt((data_1[j,2+2*i+1])**2 + (data_2[k,2+2*i+1])**2)
            data[j, 2+2*i] = data_1[j, 2+2*i]
            data[j, 2+2*i + 1] = data_1[j,2+2*i+1]

    with open(inFile1, "r") as f:
        header_txt = f.readline()
        header_txt += f.readline().strip("\n")

    np.savetxt(os.path.join(dt_folder2, "Parameter Sweep Cumulative.csv"), data, delimiter=",", header=header_txt)

    return 0

def calculate_chemical_potential(folder, skip_pts_plot=0, cumulative=False):

    if cumulative:
        file_name = "Parameter Sweep Cumulative.csv"
    else:
        file_name = "Parameter Sweep.csv"
    param_sweep_file = os.path.join(folder, file_name)
    param_sweep_data = np.loadtxt(param_sweep_file, delimiter=",", skiprows=2, unpack=False)
    energy_vals = param_sweep_data[:,-2]
    energy_err = param_sweep_data[:, -1]
    N_vals = param_sweep_data[:, 0]
    N_vals_sorted, energy_vals_sorted, energy_err_sorted = zip(*sorted(zip(N_vals, energy_vals, energy_err)))
    N_vals = np.array(N_vals_sorted)
    energy_vals = np.array(energy_vals_sorted)
    energy_err = np.array(energy_err_sorted)
    chem_pot_results = np.zeros((len(N_vals), 3))
    for j in range(len(N_vals)):
        
        '''
        chem_pot_results[j, 1] = (energy_vals[j])/N_vals[j]
        chem_pot_results[j, 2] = energy_err[j]/N_vals[j]
        chem_pot_results[j, 0] = 1.0/(N_vals[j])
        '''
        '''
        chem_pot_results[j, 1] = (energy_vals[j+5] - energy_vals[0])/(N_vals[j+5] - N_vals[0])
        chem_pot_results[j, 2] = np.sqrt(energy_err[j+5]**2 + energy_err[0]**2)/(N_vals[j+5]- N_vals[0])
        chem_pot_results[j, 0] = N_vals[j+5]
        '''
        chem_pot_results[j, 1] = (energy_vals[j])/N_vals[j]
        chem_pot_results[j, 2] = energy_err[j]/N_vals[j]
        chem_pot_results[j, 0] = N_vals[j]
    
    with open(param_sweep_file, "r") as f:
        out_header = f.readline()
    
    print(chem_pot_results[:, 0])
    out_header += "N, Chemical Potential, Chemical Potential Error"
    out_file = os.path.join(folder, "Chemical Potential.csv")
    np.savetxt(out_file, chem_pot_results, delimiter=",", header=out_header)

    plt.figure()
    plt.scatter(chem_pot_results[skip_pts_plot:-1, 0], chem_pot_results[skip_pts_plot:-1, 1])
    plt.errorbar(chem_pot_results[skip_pts_plot:-1, 0], chem_pot_results[skip_pts_plot:-1, 1], chem_pot_results[skip_pts_plot:-1, 2], fmt="None", capsize=5)
    plt.plot(chem_pot_results[skip_pts_plot:-1, 0], chem_pot_results[skip_pts_plot:-1, 1])
    plt.title("Energy dependance on N")
    plt.xlabel(r"$N$")
    plt.ylabel(r"$E$")
    plt.savefig(os.path.join(folder, "Chemical Potential.png"))

    return 0

if __name__=="__main__":

    
    # g=1.0
    dt_folder = "/Users/shaeermoeed/Github/DVRPairMC/Results/PIGS/24_08_2024_13_47_08"
    process_mc_outputs(dt_folder)
    process_estimator_outputs(dt_folder)
    process_parameter_sweeps(dt_folder, cumulative=False)
    calculate_chemical_potential(dt_folder)

    dt_folder = "/Users/shaeermoeed/Github/DVRPairMC/Results/PIGS/24_08_2024_16_00_14"
    #process_mc_outputs(dt_folder)
    #process_estimator_outputs(dt_folder)
    #process_parameter_sweeps(dt_folder, cumulative=False)
    calculate_chemical_potential(dt_folder)
    

    
    g = 0.1
    dt_folder = "/Users/shaeermoeed/Github/DVRPairMC/Results/PIGS/24_08_2024_19_00_40"
    dt_folder_2 = "/Users/shaeermoeed/Github/DVRPairMC/Results/PIGS/20_09_2024_08_54_57" #started with -1 config
    dt_folder_3 = "/Users/shaeermoeed/Github/DVRPairMC/Results/PIGS/06_10_2024_08_39_19" 
    process_mc_outputs(dt_folder)
    process_estimator_outputs(dt_folder)
    process_mc_outputs(dt_folder_2)
    process_estimator_outputs(dt_folder_2)
    combine_parameter_sweeps(dt_folder, dt_folder_2, "Parameter Sweep.csv", "Parameter Sweep.csv")
    process_parameter_sweeps(dt_folder_2, cumulative=False)
    calculate_chemical_potential(dt_folder_2, cumulative=False)
    process_mc_outputs(dt_folder_3)
    process_estimator_outputs(dt_folder_3)
    process_parameter_sweeps(dt_folder_3, cumulative=False)
    calculate_chemical_potential(dt_folder_3, cumulative=False)
    
    
    '''
    # g = 0.3
    dt_folder = "/Users/shaeermoeed/Github/DVRPairMC/Results/PIGS/24_08_2024_15_00_12"
    process_mc_outputs(dt_folder)
    process_estimator_outputs(dt_folder)
    process_parameter_sweeps(dt_folder, cumulative=False)
    calculate_chemical_potential(dt_folder)

    # g = 0.2
    dt_folder = "/Users/shaeermoeed/Github/DVRPairMC/Results/PIGS/24_08_2024_19_35_29"
    process_mc_outputs(dt_folder)
    process_estimator_outputs(dt_folder)
    process_parameter_sweeps(dt_folder, cumulative=False)
    calculate_chemical_potential(dt_folder)
    
    # g = 0.4
    dt_folder = "/Users/shaeermoeed/Github/DVRPairMC/Results/PIGS/25_08_2024_05_06_14"
    process_mc_outputs(dt_folder)
    process_estimator_outputs(dt_folder)
    process_parameter_sweeps(dt_folder, cumulative=False)
    calculate_chemical_potential(dt_folder)
    '''

    g = 0.5
    dt_folder = "/Users/shaeermoeed/Github/DVRPairMC/Results/PIGS/25_08_2024_05_58_57"
    process_mc_outputs(dt_folder)
    process_estimator_outputs(dt_folder)
    process_parameter_sweeps(dt_folder, cumulative=False)
    calculate_chemical_potential(dt_folder)
    