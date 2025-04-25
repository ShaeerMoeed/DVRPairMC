import numpy as np
import matplotlib.pyplot as plt

def parse_mc_data(file):

    mc_data = np.loadtxt(file, delimiter=",", skiprows=2)
    mc_g = mc_data[:,0]
    mc_corr_mean = mc_data[:,4]
    mc_corr_se = mc_data[:,5]
    mc_g, mc_corr_mean, mc_corr_se = zip(*sorted(zip(mc_g, mc_corr_mean, mc_corr_se)))
    mc_corr_mean = list(mc_corr_mean)
    mc_corr_se = list(mc_corr_se)

    mc_corr_mean = list(map(lambda x: x/149, mc_corr_mean))
    mc_corr_se = list(map(lambda x: x/149, mc_corr_se))

    '''
    mc_g_1 = []
    mc_corr_mean_1 = []
    mc_corr_se_1 = []
    for i in range(len(mc_g)):
        if mc_g[i] >= 0.5 and mc_g[i] <= 0.6:
            mc_g_1.append(mc_g[i])
            mc_corr_mean_1.append(mc_corr_mean[i])
            mc_corr_se_1.append(mc_corr_se[i])
    '''

    return mc_g, mc_corr_mean, mc_corr_se

dmrg_g, dmrg_corr = np.loadtxt("/Users/shaeermoeed/Github/DVRPairMC/Results/corr.txt", skiprows=4, unpack=True)
mc_results_nofit = "/Users/shaeermoeed/Github/DVRPairMC/Results/PIGS/19_07_2024_17_48_46/Parameter Sweep Cumulative No Fit.csv"
mc_g_nofit, mc_corr_mean_nofit, mc_corr_se_nofit = parse_mc_data(mc_results_nofit)
'''
plt.figure()
plt.rcParams['mathtext.fontset']='stix'
plt.plot(dmrg_g[10:], dmrg_corr[10:], label="DMRG", color="C3")
plt.scatter(mc_g_nofit, mc_corr_mean_nofit, label="MC", color="C0")
plt.errorbar(mc_g_nofit, mc_corr_mean_nofit, mc_corr_se_nofit, capsize=5, fmt="None", color='C0')
plt.ylabel(r"$C$", fontsize=20)
plt.xlabel(r"$g$", fontsize=20)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(prop={'size': 12})
plt.savefig("Correlation DMRG Comparison Unfixed.png", bbox_inches='tight', dpi=1200)
plt.close()
'''
mc_results_nofit_1 = "/Users/shaeermoeed/Github/DVRPairMC/Results/PIGS/01_08_2024_18_11_46/Parameter Sweep No Fit.csv"
mc_results_nofit_2 = "/Users/shaeermoeed/Github/DVRPairMC/Results/PIGS/20_07_2024_23_53_49/Parameter Sweep No Fit.csv"
mc_g_nofit_1, mc_corr_mean_nofit_1, mc_corr_se_nofit_1 = parse_mc_data(mc_results_nofit_1)
mc_g_nofit_2, mc_corr_mean_nofit_2, mc_corr_se_nofit_2 = parse_mc_data(mc_results_nofit_2)

plt.rcParams['mathtext.fontset']='stix'
fig, ax1 = plt.subplots()
left, bottom, width, height = [0.57, 0.17, 0.3, 0.3]
#ax2 = fig.add_axes([left, bottom, width, height])
ax1.plot(dmrg_g[10:], dmrg_corr[10:], label="DMRG", color="C3")
ax1.scatter(mc_g_nofit, mc_corr_mean_nofit, label="MC (P=60)", color="C0")
ax1.errorbar(mc_g_nofit, mc_corr_mean_nofit, mc_corr_se_nofit, capsize=5, fmt="None", color='C0')
ax1.scatter(mc_g_nofit_1, mc_corr_mean_nofit_1, label="MC (P=100)", color="black")
ax1.errorbar(mc_g_nofit_1, mc_corr_mean_nofit_1, mc_corr_se_nofit_1, capsize=5, fmt="None", color='black')
dmrg_critical_g = []
dmrg_critical_corr = []
for i in range(len(dmrg_g)):
    if dmrg_g[i] > 0.499 and dmrg_g[i] < 0.601:
        dmrg_critical_g.append(dmrg_g[i])
        dmrg_critical_corr.append(dmrg_corr[i])
mc_g_nofit_critical = []
mc_corr_nofit_critical = []
mc_corr_se_nofit_critical = []
for i in range(len(mc_g_nofit)):
    if mc_g_nofit[i] == 0.5 or mc_g_nofit[i] == 0.6:
        mc_g_nofit_critical.append(mc_g_nofit[i])
        mc_corr_nofit_critical.append(mc_corr_mean_nofit[i])
        mc_corr_se_nofit_critical.append(mc_corr_se_nofit[i])
#ax2.plot(dmrg_critical_g, dmrg_critical_corr, color="C3")
#ax2.scatter(mc_g_nofit_critical, mc_corr_nofit_critical, color="C0")
#ax2.errorbar(mc_g_nofit_critical, mc_corr_nofit_critical, mc_corr_se_nofit_critical, capsize=5, fmt="None", color='C0')
#ax2.scatter(mc_g_nofit_1, mc_corr_mean_nofit_1, label="MC (P=100)", color="black")
#ax2.errorbar(mc_g_nofit_1, mc_corr_mean_nofit_1, mc_corr_se_nofit_1, capsize=5, fmt="None", color='black')
#plt.scatter(mc_g_nofit_2, mc_corr_mean_nofit_2, label="MC (P=70)", color="C2")
#plt.errorbar(mc_g_nofit_2, mc_corr_mean_nofit_2, mc_corr_se_nofit_2, capsize=5, fmt="None", color='C2')
ax1.set_ylabel(r"$C$", fontsize=20)
ax1.set_xlabel(r"$g$", fontsize=20)
ax1.tick_params(axis='both', labelsize=12)
#handles, labels = ax1.get_legend_handles_labels()
#handles1, labels1 = ax2.get_legend_handles_labels()
#ax1.legend(handles+handles1, labels+labels1, prop={'size': 12}, loc="upper left")
ax1.legend(prop={'size': 12}, loc="upper left")
#ax2.set_xlim([0.49, 0.61])
#ax2.set_ylim([mc_corr_mean_nofit_1[0]-0.02, dmrg_critical_corr[-1]+0.02])
plt.savefig("Correlation DMRG Comparison Increased Beads.png", bbox_inches='tight', dpi=1200)
plt.close()

"""
mc_results_quadratic = "/Users/shaeermoeed/Github/DVRPairMC/Results/PIGS/19_07_2024_17_48_46/Parameter Sweep Cumulative quadratic Fit.csv"
mc_results_quadratic_1 = "/Users/shaeermoeed/Github/DVRPairMC/Results/PIGS/01_08_2024_18_11_46/Parameter Sweep quadratic Fit.csv"
mc_results_quadratic_2 = "/Users/shaeermoeed/Github/DVRPairMC/Results/PIGS/20_07_2024_23_53_49/Parameter Sweep quadratic Fit.csv"
mc_g_quadratic, mc_corr_mean_quadratic, mc_corr_se_quadratic = parse_mc_data(mc_results_quadratic)
mc_g_quadratic_1, mc_corr_mean_quadratic_1, mc_corr_se_quadratic_1 = parse_mc_data(mc_results_quadratic_1)
mc_g_quadratic_2, mc_corr_mean_quadratic_2, mc_corr_se_quadratic_2 = parse_mc_data(mc_results_quadratic_2)

plt.figure()
plt.rcParams['mathtext.fontset']='stix'
plt.plot(dmrg_g[10:], dmrg_corr[10:], label="DMRG", color="C3")
plt.scatter(mc_g_quadratic, mc_corr_mean_quadratic, label="MC (P=20,40,60)", color="C0")
plt.errorbar(mc_g_quadratic, mc_corr_mean_quadratic, mc_corr_se_quadratic, capsize=5, fmt="None", color='C0')
plt.scatter(mc_g_quadratic_1, mc_corr_mean_quadratic_1, label="MC (P=80,90,100)", color="C1")
plt.errorbar(mc_g_quadratic_1, mc_corr_mean_quadratic_1, mc_corr_se_quadratic_1, capsize=5, fmt="None", color='C1')
plt.scatter(mc_g_quadratic_2, mc_corr_mean_quadratic_2, label="MC (P=50,60,70)", color="C2")
plt.errorbar(mc_g_quadratic_2, mc_corr_mean_quadratic_2, mc_corr_se_quadratic_2, capsize=5, fmt="None", color='C2')
plt.ylabel(r"$C$", fontsize=20)
plt.xlabel(r"$g$", fontsize=20)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(prop={'size': 12})
plt.savefig("Correlation DMRG Comparison Quadratic Fit.png", bbox_inches='tight', dpi=1200)
plt.close()

mc_g_fixed = [0.5, 0.6]
mc_corr_mean_fixed = []
mc_corr_se_fixed = []
for i in range(len(mc_g_quadratic)):
    if mc_g_quadratic[i] == 0.5:
        for j in range(len(mc_g_quadratic_1)):
            if mc_g_quadratic_1[j] == 0.5:
                mc_corr_mean_fixed.append(mc_corr_mean_quadratic_1[j])
                mc_corr_se_fixed.append(mc_corr_se_quadratic_1[j])
    if mc_g_quadratic[i] == 0.6:
        for j in range(len(mc_g_quadratic_2)):
            if mc_g_quadratic_2[j] == 0.6:
                mc_corr_mean_fixed.append(mc_corr_mean_quadratic_2[j])
                mc_corr_se_fixed.append(mc_corr_se_quadratic_2[j])

plt.figure()
plt.rcParams['mathtext.fontset']='stix'
plt.plot(dmrg_g[10:], dmrg_corr[10:], label="DMRG", color="C3")
plt.scatter(mc_g_quadratic, mc_corr_mean_quadratic, label="MC Fit", color="C0")
plt.errorbar(mc_g_quadratic, mc_corr_mean_quadratic, mc_corr_se_quadratic, capsize=5, fmt="None", color='C0')
plt.scatter(mc_g_fixed, mc_corr_mean_fixed, color="black", label="MC Fit (Increased P)")
plt.errorbar(mc_g_fixed, mc_corr_mean_fixed, mc_corr_se_fixed, color="black", capsize=5, fmt="None")
plt.ylabel(r"$C$", fontsize=20)
plt.xlabel(r"$g$", fontsize=20)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(prop={'size': 12})
plt.savefig("Correlation DMRG Comparison Fixed New.png", bbox_inches='tight', dpi=1200)
plt.close()
"""