import numpy as np
from matplotlib import pyplot as plt
import os

dmrg_file = "/Users/shaeermoeed/Github/DVRPairMC/energy_dmrg.txt"
dmrg_g, dmrg_e0 = np.loadtxt(dmrg_file, skiprows=4, unpack=True)

mc_file = "/Users/shaeermoeed/Github/DVRPairMC/Results/PIGS/22_03_2024_22_08_16/Parameter Sweep.csv"
mc_data = np.loadtxt(mc_file, skiprows=2, delimiter=",")
mc_g = mc_data[:, 0]
mc_e0 = mc_data[:, -2]
mc_e0_err = mc_data[:, -1]

mc_e0 = list(map(lambda x,y :x*y, mc_e0, mc_g))
mc_e0_err = list(map(lambda x,y :x*y, mc_e0_err, mc_g))

plt.figure()
plt.rcParams['mathtext.fontset']='stix'
plt.scatter(mc_g, mc_e0, label="Pair", color="C0")
plt.errorbar(mc_g, mc_e0, mc_e0_err, fmt="None", capsize=5, color="C0")
plt.plot(dmrg_g[:-1], dmrg_e0[:-1], label="DMRG", color="C3")
plt.xlabel(r"$g$", fontsize=20)
plt.ylabel(r"$E_0$", fontsize=20)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(prop={'size': 12})
plt.savefig("Pair_Energy_g.png", bbox_inches='tight', dpi=1200)
plt.show()