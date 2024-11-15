import numpy as np
import matplotlib.pyplot as plt

dmrg_g, dmrg_corr = np.loadtxt("/Users/shaeermoeed/Github/DVRPairMC/Results/PIGS/20_07_2024_16_15_52/corr.txt", skiprows=4, unpack=True)
dmrg_g_2, dmrg_corr_2 = np.loadtxt("/Users/shaeermoeed/Github/DVRPairMC/Results/PIGS/22_03_2024_22_08_16/dmrg_corr.txt", skiprows=4, unpack=True)

mc_results = "/Users/shaeermoeed/Github/DVRPairMC/Results/PIGS/19_07_2024_17_48_46/Parameter Sweep Cumulative.csv" #300000 steps total
#mc_results = "/Users/shaeermoeed/Github/DVRPairMC/Results/PIGS/22_03_2024_22_08_16/Parameter Sweep.csv"
#mc_results = "/Users/shaeermoeed/Github/DVRPairMC/Results/PIGS/15_07_2024_05_49_33/Parameter Sweep Cumulative.csv"
mc_data = np.loadtxt(mc_results, delimiter=",", skiprows=2)
mc_g = mc_data[:,0]
mc_corr_mean = mc_data[:,4]
mc_corr_se = mc_data[:,5]

mc_resullts_2 = "/Users/shaeermoeed/Github/DVRPairMC/Results/PIGS/20_07_2024_23_53_49/Parameter Sweep.csv" # 200000 steps
mc_data_2 = np.loadtxt(mc_resullts_2, delimiter=",", skiprows=2)
mc_g_2 = mc_data_2[:,0]
mc_corr_mean_2 = mc_data_2[:,4]
mc_corr_se_2 = mc_data_2[:,5]

mc_g, mc_corr_mean, mc_corr_se = zip(*sorted(zip(mc_g, mc_corr_mean, mc_corr_se)))

#mc_corr_mean = list(map(lambda x: x/149, mc_corr_mean))
#mc_corr_se = list(map(lambda x: x/149, mc_corr_se))

mc_corr_mean = list(mc_corr_mean)
mc_corr_se = list(mc_corr_se)

mc_corr_mean_copy = np.copy(mc_corr_mean)
mc_corr_se_copy = np.copy(mc_corr_se)

'''
for i in range(len(mc_g)):
    if mc_g[i] == 0.5:
        for j in range(len(mc_g_2)):
            if mc_g_2[j] == 0.5:
                mc_corr_mean[i] = mc_corr_mean_2[j]
                mc_corr_se[i] = mc_corr_se_2[j]
                print("Updated Corr(0.5)")
'''
for i in range(len(mc_g)):
    if mc_g[i] == 0.6:
        for j in range(len(mc_g_2)):
            if mc_g_2[j] == 0.6:
                mc_corr_mean[i] = mc_corr_mean_2[j]
                mc_corr_se[i] = mc_corr_se_2[j]
                print("Updated Corr(0.6)")

plt.figure()
plt.rcParams['mathtext.fontset']='stix'

# These are in unitless percentages of the figure size. (0,0 is bottom left)
#left, bottom, width, height = [0.6, 0.25, 0.2, 0.2]
#ax2 = fig.add_axes([left, bottom, width, height])

mc_resullts_3 = "/Users/shaeermoeed/Github/DVRPairMC/Results/PIGS/01_08_2024_18_11_46/Parameter Sweep.csv" # 300000
mc_data_3 = np.loadtxt(mc_resullts_3, delimiter=",", skiprows=2)
mc_g_3 = mc_data_3[:,0]
mc_corr_mean_3 = mc_data_3[:,4]
mc_corr_se_3 = mc_data_3[:,5]

for i in range(len(mc_g)):
    if mc_g[i] == 0.5:
        for j in range(len(mc_g_3)):
            if mc_g_3[j] == 0.5:
                mc_corr_mean[i] = mc_corr_mean_3[j]
                mc_corr_se[i] = mc_corr_se_3[j]
                print("Updated Corr(0.5)")

'''
for i in range(len(mc_g)):
    if mc_g[i] == 0.6:
        for j in range(len(mc_g_3)):
            if mc_g_3[j] == 0.6:
                mc_corr_mean[i] = mc_corr_mean_3[j]
                mc_corr_se[i] = mc_corr_se_3[j]
                print("Updated Corr(0.6)")
'''

mc_corr_mean = list(map(lambda x: x/149, mc_corr_mean))
mc_corr_se = list(map(lambda x: x/149, mc_corr_se))
#ax1.plot(dmrg_g, dmrg_corr, label="DMRG (Nbond=200)", color="C3")
plt.plot(dmrg_g_2, dmrg_corr_2, label="DMRG", color="C3")
plt.scatter(mc_g, mc_corr_mean, label="MC", color="C0")
plt.errorbar(mc_g, mc_corr_mean, mc_corr_se, capsize=5, fmt="None", color="C0")
#ax1.scatter(mc_g_2, mc_corr_mean_3, label="P=80", color="C5")
#ax1.errorbar(mc_g_2, mc_corr_mean_3, mc_corr_se_3, capsize=5, color="C5")
plt.ylabel(r"$C$", fontsize=20)
plt.xlabel(r"$g$", fontsize=20)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(prop={'size': 12})
plt.savefig("Correlation DMRG Comparison Fixed.png", bbox_inches='tight', dpi=1200)

plt.figure()
mc_corr_mean_copy = list(map(lambda x: x/149, mc_corr_mean_copy))
mc_corr_se_copy = list(map(lambda x: x/149, mc_corr_se_copy))
#ax1.plot(dmrg_g, dmrg_corr, label="DMRG (Nbond=200)", color="C3")
plt.plot(dmrg_g_2, dmrg_corr_2, label="DMRG", color="C3")
plt.scatter(mc_g, mc_corr_mean_copy, label="MC", color="C0")
plt.errorbar(mc_g, mc_corr_mean_copy, mc_corr_se_copy, capsize=5, fmt="None", color="C0")
#ax1.scatter(mc_g_2, mc_corr_mean_3, label="P=80", color="C5")
#ax1.errorbar(mc_g_2, mc_corr_mean_3, mc_corr_se_3, capsize=5, color="C5")
plt.ylabel(r"$C$", fontsize=20)
plt.xlabel(r"$g$", fontsize=20)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(prop={'size': 12})
plt.savefig("Correlation DMRG Comparison Unfixed.png", bbox_inches='tight', dpi=1200)

plt.close()

'''
ax2.plot(dmrg_g_2, dmrg_corr_2, label="DMRG", color="C3")
ax2.scatter(mc_g[3:6], mc_corr_mean_copy[3:6], label="MC", color="C0")
ax2.errorbar(mc_g[3:6], mc_corr_mean_copy[3:6], mc_corr_se_copy[3:6], capsize=5, fmt="None")
ax2.set_ylabel(r"$C$", fontsize=20)
ax2.set_xlabel(r"$g$", fontsize=20)
ax2.set_yticks([0.0, 0.3, 0.6])
ax2.set_xticks([0.0, 0.5, 1,0])
'''

#ax1.set_title("Effect of Coupling Strength On Structural Properties")



