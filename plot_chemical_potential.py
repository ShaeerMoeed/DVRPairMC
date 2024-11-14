import numpy as np
import matplotlib.pyplot as plt

g_0_file = "/Users/shaeermoeed/Github/DVRPairMC/Results/PIGS/06_10_2024_08_39_19/Chemical Potential.csv"
g_1_file = "/Users/shaeermoeed/Github/DVRPairMC/Results/PIGS/25_08_2024_05_58_57/Chemical Potential.csv"
g_2_file = "/Users/shaeermoeed/Github/DVRPairMC/Results/PIGS/24_08_2024_13_47_08/Chemical Potential.csv"

N_2_results_file = "/Users/shaeermoeed/Github/DVRPairMC/E0_g_N2.csv"
N_2_g, N_2_E = np.loadtxt(N_2_results_file, unpack=True, delimiter=",")
N_2_E = 0.5 * N_2_E

N_0,mu_0,mu_0_err = np.loadtxt(g_0_file, unpack=True, delimiter=",")
N_1,mu_1,mu_1_err = np.loadtxt(g_1_file, unpack=True, delimiter=",")
N_2,mu_2,mu_2_err = np.loadtxt(g_2_file, unpack=True, delimiter=",")

N_0_list = [2]
mu_0_list = [N_2_E[0]]
mu_0_err_list = [0.0]
for i in range(len(N_0)):
    N = N_0[i]
    if N > 35:
        break
    N_0_list.append(N)
    mu_0_list.append(mu_0[i])
    mu_0_err_list.append(mu_0_err[i])

N_1_list = [2]
mu_1_list = [N_2_E[1]]
mu_1_err_list = [0.0]
for i in range(len(N_1)):
    N = N_1[i]
    if N > 35:
        break
    N_1_list.append(N)
    mu_1_list.append(mu_1[i])
    mu_1_err_list.append(mu_1_err[i])

N_2_list = [2]
mu_2_list = [N_2_E[2]]
mu_2_err_list = [0.0]
for i in range(len(N_2)):
    N = N_2[i]
    if N > 35:
        break
    N_2_list.append(N)
    mu_2_list.append(mu_2[i])
    mu_2_err_list.append(mu_2_err[i])

plt.figure()
plt.rcParams['mathtext.fontset']='stix'
plt.scatter(N_0_list, mu_0_list, label="g=0.1", color="C0")
plt.errorbar(N_0_list, mu_0_list, mu_0_err_list, fmt="None", capsize=5, color="C0")
plt.scatter(N_1_list, mu_1_list, label="g=0.5", color="C3")
plt.errorbar(N_1_list, mu_1_list, mu_1_err_list, fmt="None", capsize=5, color="C3")
plt.scatter(N_2_list, mu_2_list, label="g=1.0", color='black')
plt.errorbar(N_2_list, mu_2_list, mu_2_err_list, fmt="None", capsize=5, color='black')
plt.legend(prop={'size': 12})
plt.ylabel(r"$E/N$", fontsize=20)
plt.xlabel(r"$N$", fontsize=20)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.savefig("chemical_potential.png", bbox_inches='tight', dpi=1200)
plt.show()
