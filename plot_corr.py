import numpy as np
import matplotlib.pyplot as plt

dmrg_g, dmrg_corr = np.loadtxt("/Users/shaeermoeed/Github/DVRPairMC/Results/PIGS/22_03_2024_21_58_01/corr.txt", skiprows=4, unpack=True)
mc_results = "/Users/shaeermoeed/Github/DVRPairMC/Results/PIGS/22_03_2024_21_58_01/Parameter Sweep.csv"

mc_data = np.loadtxt(mc_results, delimiter=",", skiprows=2)
mc_g = mc_data[:,0]
mc_corr_mean = mc_data[:,4]
mc_corr_se = mc_data[:,5]

mc_binder_mean = mc_data[:,6]
mc_binder_se = mc_data[:,7]

mc_g, mc_corr_mean, mc_binder_mean, mc_corr_se, mc_binder_se = zip(*sorted(zip(mc_g, mc_corr_mean, mc_binder_mean, mc_corr_se, mc_binder_se)))

mc_corr_mean = list(map(lambda x: x/100.0, mc_corr_mean))
mc_corr_se = list(map(lambda x: x/100.0, mc_corr_se))

fig, ax1 = plt.subplots()

# These are in unitless percentages of the figure size. (0,0 is bottom left)
left, bottom, width, height = [0.6, 0.25, 0.2, 0.2]
ax2 = fig.add_axes([left, bottom, width, height])

ax1.plot(dmrg_g[:10], dmrg_corr[:10], label="DMRG", color="C3")
ax1.scatter(mc_g, mc_corr_mean, label="MC", color="C0")
ax1.errorbar(mc_g, mc_corr_mean, mc_corr_se, capsize=5)
ax1.legend()
ax1.set_ylabel("C")
ax1.set_xlabel("g")

ax2.plot(list(mc_g), list(mc_binder_mean), color="C0")
ax2.scatter(mc_g, mc_binder_mean, color="C0")
ax2.errorbar(mc_g, mc_binder_mean, mc_binder_se, capsize=5)
ax2.set_ylabel("U")
ax2.set_xlabel("g")
ax2.set_yticks([0.0, 0.3, 0.6])
ax2.set_xticks([0.0, 0.5, 1,0])

ax1.set_title("Effect of Coupling Strength On Structural Properties")

plt.savefig("Correlation DMRG Comparison.png", bbox_inches='tight')



