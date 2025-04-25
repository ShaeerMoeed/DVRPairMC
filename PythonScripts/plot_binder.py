import numpy as np
import matplotlib.pyplot as plt

mc_results = "/Users/shaeermoeed/Github/DVRPairMC/Results/PIGS/22_03_2024_22_08_16/Parameter Sweep No Fit.csv"
#mc_results = "/Users/shaeermoeed/Github/DVRPairMC/Results/PIGS/15_07_2024_05_49_33/Parameter Sweep Cumulative.csv"

mc_data = np.loadtxt(mc_results, delimiter=",", skiprows=2)
mc_g = mc_data[:,0]

mc_binder_mean = mc_data[:,6]
mc_binder_se = mc_data[:,7]

'''
mc_resullts_2 = "/Users/shaeermoeed/Github/DVRPairMC/Results/PIGS/20_07_2024_23_53_49/Parameter Sweep.csv"
mc_data_2 = np.loadtxt(mc_resullts_2, delimiter=",", skiprows=2)
mc_g_2 = mc_data_2[:,0]
mc_binder_mean_2 = mc_data_2[:,6]
mc_binder_se_2 = mc_data_2[:,7]

mc_g, mc_binder_mean, mc_binder_se = zip(*sorted(zip(mc_g, mc_binder_mean, mc_binder_se)))

mc_binder_mean = list(mc_binder_mean)
mc_binder_se = list(mc_binder_se)

mc_resullts_3 = "/Users/shaeermoeed/Github/DVRPairMC/Results/PIGS/01_08_2024_18_11_46/Parameter Sweep.csv"
mc_data_3 = np.loadtxt(mc_resullts_3, delimiter=",", skiprows=2)
mc_g_3 = mc_data_3[:,0]
mc_binder_mean_3 = mc_data_3[:,6]
mc_binder_se_3 = mc_data_3[:,7]
'''

'''
for i in range(len(mc_g)):
    if (mc_g[i] == 0.5):
        mc_binder_mean[i] = mc_binder_mean_2[0]
        mc_binder_se[i] = mc_binder_se_2[0]
    #if (mc_g[i] == 0.6):
    #    mc_binder_mean[i] = mc_binder_mean_2[1]
    #    mc_binder_se[i] = mc_binder_se_2[1]
'''

'''
for i in range(len(mc_g)):
    if (mc_g[i] == 0.5):
        mc_binder_mean[i] = mc_binder_mean_3[1]
        mc_binder_se[i] = mc_binder_se_3[1]
    if (mc_g[i] == 0.6):
        mc_binder_mean[i] = mc_binder_mean_3[0]
        mc_binder_se[i] = mc_binder_se_3[0]
'''

dmrg_g, dmrg_pol, dmrg_binder = np.loadtxt("/Users/shaeermoeed/Github/DVRPairMC/PythonScripts/binder_x.txt", skiprows=4, unpack=True)

plt.figure()
#ax2.plot(list(mc_g), list(mc_binder_mean), color="C0")
plt.rcParams['mathtext.fontset']='stix'
plt.scatter(mc_g, mc_binder_mean, color="C0", label="MC")
plt.errorbar(mc_g, mc_binder_mean, mc_binder_se, capsize=5, fmt="None")
plt.plot(dmrg_g, dmrg_binder, color="C3", label="DMRG")
plt.ylabel(r"$U$", fontsize=20)
plt.xlabel(r"$g$", fontsize=20)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(prop={'size': 12})
#ax2.set_yticks([0.0, 0.3, 0.6])
#ax2.set_xticks([0.0, 0.5, 1,0])
plt.savefig("Binder DMRG Comparison.png", bbox_inches='tight', dpi=1200)





