import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib import animation
import os

def plotSpecificTraj(in_file, bead_num, rotor_num):

    diag_data = np.loadtxt(in_file, skiprows=2, delimiter=",")
    numSteps = diag_data[-1, 0]
    numBeads = diag_data[-1, 1]
    numRotors = diag_data[-1, 2]
    if (bead_num < 1) or (bead_num > numBeads):
        raise Exception("Incorrect bead number") 
    if (rotor_num < 1) or (rotor_num > numRotors):
        raise Exception("Incorrect rotor number") 
    arr_size = np.shape(diag_data)[0]
    traj = np.zeros(int(numSteps))
    count = 0
    for i in range(arr_size):
        if diag_data[i,1] == bead_num:
            if diag_data[i,2] == rotor_num:
                traj[count] = diag_data[i,3]
                count += 1

    plt.figure()
    plt.plot(traj)
    plt.show()

    return 0

def plotSpecificBeadAndTime(in_file, sim_time, bead_num):

    diag_data = np.loadtxt(in_file, skiprows=2, delimiter=",")
    numSteps = diag_data[-1, 0]
    numBeads = diag_data[-1, 1]
    numRotors = diag_data[-1, 2]
    if (bead_num < 1) or (bead_num > numBeads):
        raise Exception("Incorrect bead number") 
    if (sim_time < 1) or (sim_time > numSteps):
        raise Exception("Incorrect time step number") 
    arr_size = np.shape(diag_data)[0]
    beadConfigT = np.zeros(int(numRotors))
    count = 0
    for i in range(arr_size):
        if diag_data[i,0] == sim_time:
            if diag_data[i,1] == bead_num:
                beadConfigT[count] = diag_data[i,3]
                count += 1

    plt.figure()
    plt.plot(beadConfigT)
    plt.show()

    return 0

def buildConfigTensor(in_file, numStepsToUse=100):

    diag_data = np.loadtxt(in_file, skiprows=2, delimiter=",", dtype=int)
    #diag_data = diag_data[:, :]
    numSteps = diag_data[-1, 0]
    numBeads = diag_data[-1, 1]
    numRotors = diag_data[-1, 2]
    configTensor = np.zeros((numStepsToUse, numBeads, numRotors), dtype=np.float64)
    #arr_size = np.shape(diag_data)[0]
    conversion_factor = ((2 * np.pi)/(21.0))
    for i in range(numStepsToUse):
        #step_number = i+1
        for j in range(numBeads):
            #bead_number = j+1
            for k in range(numRotors):
                #rotor_number = k+1
                configTensor[i, j, k] = np.cos(conversion_factor * diag_data[(i + numSteps - numStepsToUse) * (numBeads * numRotors) + (j * numRotors) + k, 3])

    return configTensor, numBeads, numRotors, numSteps

def plot_animation(configs_tensor, numBeads, numRotors, out_folder, num_frames, factor=50):

    fig, ax = plt.subplots()
    def init_heatmap():
        ax.clear()
        sns.heatmap(np.zeros((numBeads, numRotors)), ax=ax, vmax=1.0, vmin=-1.0, cbar=True, annot=False, fmt=".2f")
        return fig, ax

    def update_heatmap(i):
        ax.clear()
        sns.heatmap(configs_tensor[i*factor, :, :], ax=ax, vmax=1.0, vmin=-1.0, cbar=False, annot=False, fmt=".2f")

    anim = animation.FuncAnimation(fig, update_heatmap, init_func=init_heatmap, frames=int(num_frames/factor), repeat=False)
    anim.save(os.path.join(out_folder, 'animated_heatmap.mp4'), writer='ffmpeg')
    #plt.show()
    plt.close()

    return 0

g = 0.5
#diagFolder = "/Users/shaeermoeed/Github/DVRPairMC/Results/PIGS/28_08_2024_17_14_21/g_{}_T_0.1_P_60_N_10_l_10/Block_1".format(g)
diagFolder = "/Users/shaeermoeed/Github/DVRPairMC/Results/PIGS/28_08_2024_17_59_09/g_{}_T_0.1_P_60_N_150_l_10/Block_1".format(g)
#diagFolder = "/Users/shaeermoeed/Github/DVRPairMC/Results/PIGS/29_08_2024_08_47_51/g_{}_T_0.1_P_60_N_150_l_10/Block_1".format(g)
#diagFolder = "/Users/shaeermoeed/Github/DVRPairMC/Results/PIGS/29_08_2024_11_06_23/g_{}_T_0.1_P_60_N_150_l_10/Block_1".format(g)
final_config_file = os.path.join(diagFolder,"Final Config.csv")

numStepsToUse = 10000
configs, numBeads, numRotors, numSteps = buildConfigTensor(os.path.join(diagFolder, "MC Diagnostic Outputs.csv"), numStepsToUse)
#plotSpecificTraj(diagFile, 27, 10)
#plotSpecificBeadAndTime(diagFile, 4998, 31)
plot_animation(configs, numBeads, numRotors, diagFolder, numStepsToUse)
plt.figure()
final_config_file = os.path.join(diagFolder,"Final Config.csv")
np.savetxt(final_config_file, configs[-1, :, :], delimiter=",")

numRotors = 150
data_to_plot = np.loadtxt(final_config_file, delimiter=",", skiprows=0)
plt.bar(np.arange(1, numRotors+1), (data_to_plot[30, :]), width=0.5)
plt.ylabel(r"cos($\phi$)")
plt.xlabel("n")
plt.title("g = {}".format(g))
plt.savefig(os.path.join(diagFolder, "Final Configuration Middle Bead (g = {}).png".format(g)))
plt.show()

