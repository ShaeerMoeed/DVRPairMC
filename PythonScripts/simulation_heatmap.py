import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
from matplotlib import gridspec

def buildConfigMatrix(diag_data, mc_step, numStepsToUse=100):

    numSteps = diag_data[-1, 0]
    numBeads = diag_data[-1, 1]
    numRotors = diag_data[-1, 2]
    config_Step = np.zeros((numBeads, numRotors), dtype=np.float64)
    conversion_factor = ((2 * np.pi)/(21.0))
    for j in range(numBeads):
        for k in range(numRotors):
            config_Step[j, k] = np.cos(conversion_factor * diag_data[(mc_step + numSteps - numStepsToUse) * (numBeads * numRotors) + (j * numRotors) + k, 3])

    return config_Step, numBeads, numRotors, numSteps

plt.rcParams['mathtext.fontset']='stix'
fig = plt.figure(constrained_layout=False)

numStepsToUse = 100000
#numStepsToUse = 1000
gs = fig.add_gridspec(3,3,bottom=0.22)

#diagFolder = "/Users/shaeermoeed/Github/DVRPairMC/Results/PIGS/28_08_2024_17_14_21/g_0.1_T_0.1_P_60_N_10_l_10/Block_1"

diagFolder = "/Users/shaeermoeed/Github/DVRPairMC/Results/PIGS/29_08_2024_08_47_51/g_0.2_T_0.1_P_60_N_150_l_10/Block_1"
inFile = os.path.join(diagFolder, "MC Diagnostic Outputs.csv")
diag_data = np.loadtxt(inFile, skiprows=2, delimiter=",", dtype=int)

mcStep_1 = 1000
configs, numBeads, numRotors, numSteps = buildConfigMatrix(diag_data, mcStep_1, numStepsToUse)

cmap_selection = 'Spectral'

ax_1 = plt.subplot(gs[0])
sns.heatmap(configs, ax=ax_1, vmax=1.0, vmin=-1.0, cbar=False, annot=False, fmt=".2f", cmap=cmap_selection)
ax_1.set_xticks([])
ax_1.set_yticks([])

ax_1.text(0.5, 1.25, r'$t_{sim} \sim 10^3$',horizontalalignment='center', verticalalignment='top', transform=ax_1.transAxes, fontsize=15)

colorbar_data = ax_1.collections[0]

mcStep_2 = 10000
configs, numBeads, numRotors, numSteps = buildConfigMatrix(diag_data, mcStep_2, numStepsToUse)

ax_2 = plt.subplot(gs[1])
sns.heatmap(configs, ax=ax_2, vmax=1.0, vmin=-1.0, cbar=False, annot=False, fmt=".2f", cmap=cmap_selection)
ax_2.set_xticks([])
ax_2.set_yticks([])

ax_2.text(0.5, 1.25, r'$t_{sim} \sim 10^4$',horizontalalignment='center', verticalalignment='top', transform=ax_2.transAxes, fontsize=15)

mcStep_3 = numSteps-100
#mcStep_3 = 100
configs, numBeads, numRotors, numSteps = buildConfigMatrix(diag_data, mcStep_3, numStepsToUse)
print("g=0.2")
print(numSteps)

ax_3 = plt.subplot(gs[2])
sns.heatmap(configs, ax=ax_3, vmax=1.0, vmin=-1.0, cbar=False, annot=False, fmt=".2f", cmap=cmap_selection)
ax_3.set_xticks([])
ax_3.set_yticks([])

ax_3.text(1.4, 0.5, r'$g < g_c$',horizontalalignment='right', verticalalignment='center', transform=ax_3.transAxes, fontsize=15)
ax_3.text(0.5, 1.25, r'$t_{sim} \sim 10^5$',horizontalalignment='center', verticalalignment='top', transform=ax_3.transAxes, fontsize=15)

diagFolder = "/Users/shaeermoeed/Github/DVRPairMC/Results/PIGS/28_08_2024_17_59_09/g_0.5_T_0.1_P_60_N_150_l_10/Block_1"
inFile = os.path.join(diagFolder, "MC Diagnostic Outputs.csv")
diag_data = np.loadtxt(inFile, skiprows=2, delimiter=",", dtype=int)

configs, numBeads, numRotors, numSteps = buildConfigMatrix(diag_data, mcStep_1, numStepsToUse)
ax_1 = plt.subplot(gs[3])
sns.heatmap(configs, ax=ax_1, vmax=1.0, vmin=-1.0, cbar=False, annot=False, fmt=".2f", cmap=cmap_selection)
ax_1.set_xticks([])
ax_1.set_yticks([])

configs, numBeads, numRotors, numSteps = buildConfigMatrix(diag_data, mcStep_2, numStepsToUse)

ax_2 = plt.subplot(gs[4])
sns.heatmap(configs, ax=ax_2, vmax=1.0, vmin=-1.0, cbar=False, annot=False, fmt=".2f", cmap=cmap_selection)
ax_2.set_xticks([])
ax_2.set_yticks([])

configs, numBeads, numRotors, numSteps = buildConfigMatrix(diag_data, mcStep_3, numStepsToUse)
print("g=0.5")
print(numSteps)

ax_3 = plt.subplot(gs[5])
sns.heatmap(configs, ax=ax_3, vmax=1.0, vmin=-1.0, cbar=False, annot=False, fmt=".2f", cmap=cmap_selection)
ax_3.set_xticks([])
ax_3.set_yticks([])

ax_3.text(1.4, 0.5, r'$g \approx g_c$',horizontalalignment='right', verticalalignment='center', transform=ax_3.transAxes, fontsize=15)

diagFolder = "/Users/shaeermoeed/Github/DVRPairMC/Results/PIGS/29_08_2024_11_06_23/g_1.0_T_0.1_P_60_N_150_l_10/Block_1"
inFile = os.path.join(diagFolder, "MC Diagnostic Outputs.csv")
diag_data = np.loadtxt(inFile, skiprows=2, delimiter=",", dtype=int)

configs, numBeads, numRotors, numSteps = buildConfigMatrix(diag_data, mcStep_1, numStepsToUse)
ax_1 = plt.subplot(gs[6])
sns.heatmap(configs, ax=ax_1, vmax=1.0, vmin=-1.0, cbar=False, annot=False, fmt=".2f", cmap=cmap_selection)
ax_1.set_xticks([])
ax_1.set_yticks([])

ax_1.text(0.0,-0.18, r'$\longrightarrow$', horizontalalignment='left', verticalalignment='bottom', transform=ax_1.transAxes, fontsize=15)
ax_1.text(0.15,-0.34, r'$N$', horizontalalignment='left', verticalalignment='bottom', transform=ax_1.transAxes, fontsize=15)

ax_1.text(-0.11,0.0, r'$\longrightarrow$', horizontalalignment='left', verticalalignment='bottom', transform=ax_1.transAxes, fontsize=15, rotation=90)
ax_1.text(-0.19,0.26, r'$P$', horizontalalignment='left', verticalalignment='bottom', transform=ax_1.transAxes, fontsize=15, rotation=90)

configs, numBeads, numRotors, numSteps = buildConfigMatrix(diag_data, mcStep_2, numStepsToUse)

ax_2 = plt.subplot(gs[7])
sns.heatmap(configs, ax=ax_2, vmax=1.0, vmin=-1.0, cbar=False, annot=False, fmt=".2f", cmap=cmap_selection)
ax_2.set_xticks([])
ax_2.set_yticks([])

configs, numBeads, numRotors, numSteps = buildConfigMatrix(diag_data, mcStep_3, numStepsToUse)
print("g=1.0")
print(numSteps)

ax_3 = plt.subplot(gs[8])
sns.heatmap(configs, ax=ax_3, vmax=1.0, vmin=-1.0, cbar=False, annot=False, fmt=".2f", cmap=cmap_selection)
ax_3.set_xticks([])
ax_3.set_yticks([])

ax_3.text(1.4, 0.5, r'$g > g_c$',horizontalalignment='right', verticalalignment='center', transform=ax_3.transAxes, fontsize=15)

gs1 = fig.add_gridspec(nrows=1, ncols=1, top=0.15)
ax4 = plt.subplot(gs1[0])
plt.colorbar(colorbar_data, cax=ax4, orientation='horizontal')
ax4.set_xlabel(r'$cos(\phi)$', loc='right', fontsize=15)
plt.savefig(os.path.join(diagFolder, "heatmap_test_6.pdf"), dpi=1500)
#plt.show()