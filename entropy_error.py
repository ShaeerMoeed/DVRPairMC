import numpy as np
import matplotlib.pyplot as plt

def read_data(file, skip_min=None):

    step_number, N, D = np.loadtxt(file, skiprows=2, unpack=True, delimiter=",")

    skip = 0
    for i in range(len(N)):
        if N[i] == 0:
            skip += 1
        elif D[i] == 0:
            skip += 1

    if skip_min == None:
        skip_min = 1
    skip = max(skip, skip_min)

    N = N[skip:]
    D = D[skip:]
    step_number = step_number[skip:]

    '''
    h = np.divide(N,step_number)
    std_dev_old = np.std(h[0:10000])
    std_dev_vec = np.zeros(len(h))
    std_dev_vec[0] = std_dev_old
    for i in range(1,len(h)):
        std_dev_new = np.std(h[i:i+10000])
        if (np.abs(std_dev_new-std_dev_old) < 1e-6):
            skip = i+10000
            break
        std_dev_old = std_dev_new
    '''

    '''
    for i in range(1,len(h)-10000):
        std_dev_vec[i] = np.std(h[i:i+10000])

    plt.figure()
    plt.plot(np.diff(std_dev_vec[skip:]))

    plt.figure()
    plt.plot(h)
    plt.plot([skip] * len(h), np.linspace(min(h), max(h), len(h)))
    plt.show()
        
    N = N[skip:]
    D = D[skip:]
    step_number = step_number[skip:]
    '''

    return N,D,step_number,skip

file1 = "/Users/shaeermoeed/Github/DVRPairMC/Results/Entropy/Test/12_06_2024_04_49_29/g_1.0_T_0.1_P_20_NA_2_NB_2_l_10_Start_Zero/MC Step Outputs.csv"
N1, D1, step_number_1, skip1 = read_data(file1, 10000)
file2 = "/Users/shaeermoeed/Github/DVRPairMC/Results/Entropy/Test/12_06_2024_04_49_29/g_1.0_T_0.1_P_20_NA_2_NB_2_l_10_Start_Pi/MC Step Outputs.csv"
N2, D2, step_number_2, skip2 = read_data(file2, 10000)
h1 = np.divide(N1, step_number_1)
h2 = np.divide(N2, step_number_2)
print(np.shape(h1))
print(np.shape(h2))
print(skip1)
print(skip2)
h = np.concatenate((h1, h2), axis=0)

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

    return func_mean, func_std_error_jk, func_array, j_mean_array

def func_entropy(x):

    return -np.log(x/(1.0 - x))

print(len(h))
samples_1, num_blocks_1 = jack_knife_samples(20, h1)
samples_2, num_blocks_2 = jack_knife_samples(20, h2)

num_blocks = num_blocks_1 + num_blocks_2
samples = np.concatenate((samples_1, samples_2), axis=0)
S2_jk_mean, S2_jk_err, S2_arr, j_mean_array = calculate_jk_estimates(samples, func_entropy, num_blocks)

plt.figure()
plt.hist(j_mean_array)
plt.title("Jack-Knife Means (<h>)")
plt.savefig("/Users/shaeermoeed/Github/DVRPairMC/Jack-Knife Means.png")

print(S2_jk_mean - S2_jk_err, S2_jk_mean + S2_jk_err)
print(S2_jk_mean)

plt.figure()
plt.plot(np.arange(len(S2_arr)), S2_arr, label="S2")
plt.plot(np.arange(len(S2_arr)), [S2_jk_mean - S2_jk_err] * len(S2_arr), label="- Error")
plt.plot(np.arange(len(S2_arr)), [S2_jk_mean + S2_jk_err] * len(S2_arr), label = "+ Error")
plt.legend()
plt.title("Jack-Knife S2 (<h>/(1-<h>))")
plt.savefig("/Users/shaeermoeed/Github/DVRPairMC/Jack-Knife S2.png")

plt.figure()
plt.plot(np.arange(len(samples)), samples, label="h")
#plt.legend()
plt.title("h")
plt.savefig("/Users/shaeermoeed/Github/DVRPairMC/h MC.png")
plt.show()