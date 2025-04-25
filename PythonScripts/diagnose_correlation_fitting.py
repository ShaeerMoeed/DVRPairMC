import numpy as np
import matplotlib.pyplot as plt

dmrg_g_2, dmrg_corr_2 = np.loadtxt("/Users/shaeermoeed/Github/DVRPairMC/Results/PIGS/22_03_2024_22_08_16/dmrg_corr.txt", skiprows=4, unpack=True)
dmrg_g = []
dmrg_corr = []
'''
for i in range(len(dmrg_g_2)):
    if dmrg_g_2[i] >= 0.5 and dmrg_g_2[i] <= 0.6:
        dmrg_g.append(dmrg_g_2[i])
        dmrg_corr.append(dmrg_corr_2[i])
'''
dmrg_g = dmrg_g_2
dmrg_corr = dmrg_corr_2

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

mc_resullts_2_linear = "/Users/shaeermoeed/Github/DVRPairMC/Results/PIGS/20_07_2024_23_53_49/Parameter Sweep linear.csv" # 200000 steps
mc_resullts_2_quadratic = "/Users/shaeermoeed/Github/DVRPairMC/Results/PIGS/20_07_2024_23_53_49/Parameter Sweep quadratic.csv" # 200000 steps
mc_resullts_2_nofit = "/Users/shaeermoeed/Github/DVRPairMC/Results/PIGS/20_07_2024_23_53_49/Parameter Sweep.csv" # 200000 steps

mc_g_2_linear, mc_corr_mean_2_linear, mc_corr_se_2_linear = parse_mc_data(mc_resullts_2_linear)
mc_g_2_quadratic, mc_corr_mean_2_quadratic, mc_corr_se_2_quadratic = parse_mc_data(mc_resullts_2_quadratic)
mc_g_2_nofit, mc_corr_mean_2_nofit, mc_corr_se_2_nofit = parse_mc_data(mc_resullts_2_nofit)

mc_resullts_1_linear = "/Users/shaeermoeed/Github/DVRPairMC/Results/PIGS/01_08_2024_18_11_46/Parameter Sweep linear.csv" # 300000
mc_resullts_1_quadratic = "/Users/shaeermoeed/Github/DVRPairMC/Results/PIGS/01_08_2024_18_11_46/Parameter Sweep quadratic.csv" # 300000
mc_resullts_1_nofit = "/Users/shaeermoeed/Github/DVRPairMC/Results/PIGS/01_08_2024_18_11_46/Parameter Sweep.csv" # 300000

mc_g_1_linear, mc_corr_mean_1_linear, mc_corr_se_1_linear = parse_mc_data(mc_resullts_1_linear)
mc_g_1_quadratic, mc_corr_mean_1_quadratic, mc_corr_se_1_quadratic = parse_mc_data(mc_resullts_1_quadratic)
mc_g_1_nofit, mc_corr_mean_1_nofit, mc_corr_se_1_nofit = parse_mc_data(mc_resullts_1_nofit)

plt.figure()
plt.rcParams['mathtext.fontset']='stix'
plt.plot(dmrg_g, dmrg_corr, label="DMRG", color="C10")

plt.scatter(mc_g_2_linear[1], mc_corr_mean_2_linear[1], label="MC (0.6) Linear", color='C1')
plt.errorbar(mc_g_2_linear[1], mc_corr_mean_2_linear[1], mc_corr_se_2_linear[1], capsize=10, fmt="None", color='C1')

plt.scatter(mc_g_2_quadratic[1], mc_corr_mean_2_quadratic[1], label="MC (0.6) Quadratic", color='C2')
plt.errorbar(mc_g_2_quadratic[1], mc_corr_mean_2_quadratic[1], mc_corr_se_2_quadratic[1], capsize=5, fmt="None", color='C2')

#plt.scatter(mc_g_2_nofit[1], mc_corr_mean_2_nofit[1], label="MC (0.6) No Fit", color='C6')
#plt.errorbar(mc_g_2_nofit[1], mc_corr_mean_2_nofit[1], mc_corr_se_2_nofit[1], capsize=5, fmt="None", color='C6')

plt.scatter(mc_g_1_linear[0], mc_corr_mean_1_linear[0], label="MC (0.5) Linear", color='C3')
plt.errorbar(mc_g_1_linear[0], mc_corr_mean_1_linear[0], mc_corr_se_1_linear[0], capsize=10, fmt="None", color='C3')

plt.scatter(mc_g_1_quadratic[0], mc_corr_mean_1_quadratic[0], label="MC (0.5) Quadratic", color='C4')
plt.errorbar(mc_g_1_quadratic[0], mc_corr_mean_1_quadratic[0], mc_corr_se_1_quadratic[0], capsize=5, fmt="None", color='C4')

#plt.scatter(mc_g_1_nofit[0], mc_corr_mean_1_nofit[0], label="MC (0.5) No Fit", color='C8')
#plt.errorbar(mc_g_1_nofit[0], mc_corr_mean_1_nofit[0], mc_corr_se_1_nofit[0], capsize=5, fmt="None", color='C8')

plt.ylabel(r"$C$", fontsize=20)
plt.xlabel(r"$g$", fontsize=20)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(prop={'size': 12})
plt.savefig("Correlation Fits Comparison 1.png", bbox_inches='tight', dpi=1200)

'''
plt.figure()
plt.rcParams['mathtext.fontset']='stix'
plt.plot(dmrg_g, dmrg_corr, label="DMRG", color="C10")

plt.scatter(mc_g_2_nofit, mc_corr_mean_2_nofit, label="MC (0.6) No Fit", color='C6')
plt.errorbar(mc_g_2_nofit, mc_corr_mean_2_nofit, mc_corr_se_2_nofit, capsize=5, fmt="None", color='C6')

plt.scatter(mc_g_1_nofit, mc_corr_mean_1_nofit, label="MC (0.5) No Fit", color='C8')
plt.errorbar(mc_g_1_nofit, mc_corr_mean_1_nofit, mc_corr_se_1_nofit, capsize=5, fmt="None", color='C8')

plt.ylabel(r"$C$", fontsize=20)
plt.xlabel(r"$g$", fontsize=20)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(prop={'size': 12})
plt.savefig("Correlation Data Comparison.png", bbox_inches='tight', dpi=1200)
'''

mc_results_3_nofit = "/Users/shaeermoeed/Github/DVRPairMC/Results/PIGS/diagnosing_correlation/Parameter Sweep.csv"
mc_results_3_linear = "/Users/shaeermoeed/Github/DVRPairMC/Results/PIGS/diagnosing_correlation/Parameter Sweep linear.csv"
mc_results_3_quadratic = "/Users/shaeermoeed/Github/DVRPairMC/Results/PIGS/diagnosing_correlation/Parameter Sweep quadratic.csv"
mc_results_3_lin_plus_quad = "/Users/shaeermoeed/Github/DVRPairMC/Results/PIGS/diagnosing_correlation/Parameter Sweep lin_plus_quad.csv"
mc_g_3_linear, mc_corr_mean_3_linear, mc_corr_se_3_linear = parse_mc_data(mc_results_3_linear)
mc_g_3_quadratic, mc_corr_mean_3_quadratic, mc_corr_se_3_quadratic = parse_mc_data(mc_results_3_quadratic)
mc_g_3_nofit, mc_corr_mean_3_nofit, mc_corr_se_3_nofit = parse_mc_data(mc_results_3_nofit)
mc_g_3_lin_plus_quad, mc_corr_mean_3_lin_plus_quad, mc_corr_se_3_lin_plus_quad = parse_mc_data(mc_results_3_lin_plus_quad)

plt.figure()
plt.rcParams['mathtext.fontset']='stix'
plt.plot(dmrg_g, dmrg_corr, label="DMRG", color="C10")

plt.scatter(mc_g_3_nofit, mc_corr_mean_3_nofit, label="MC No Fit", color='C0')
plt.errorbar(mc_g_3_nofit, mc_corr_mean_3_nofit, mc_corr_se_3_nofit, capsize=5, fmt="None", color='C0')

plt.scatter(mc_g_3_linear, mc_corr_mean_3_linear, label="MC linear", color='C1')
plt.errorbar(mc_g_3_linear, mc_corr_mean_3_linear, mc_corr_se_3_linear, capsize=5, fmt="None", color='C1')

plt.scatter(mc_g_3_quadratic, mc_corr_mean_3_quadratic, label="MC quadratic", color='C2')
plt.errorbar(mc_g_3_quadratic, mc_corr_mean_3_quadratic, mc_corr_se_3_quadratic, capsize=5, fmt="None", color='C2')

#plt.scatter(mc_g_3_lin_plus_quad, mc_corr_mean_3_lin_plus_quad, label="MC lin_plus_quad", color='C3')
#plt.errorbar(mc_g_3_lin_plus_quad, mc_corr_mean_3_lin_plus_quad, mc_corr_se_3_lin_plus_quad, capsize=5, fmt="None", color='C3')

plt.ylabel(r"$C$", fontsize=20)
plt.xlabel(r"$g$", fontsize=20)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(prop={'size': 12})
plt.savefig("Correlation Fits Comparison 2.png", bbox_inches='tight', dpi=1200)