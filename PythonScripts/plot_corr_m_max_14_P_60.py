import numpy as np
import matplotlib.pyplot as plt

def parse_mc_data(file, num_sites=149):

    mc_data = np.loadtxt(file, delimiter=",", skiprows=2)
    mc_g = mc_data[:,0]
    mc_corr_mean = mc_data[:,4]
    mc_corr_se = mc_data[:,5]
    mc_g, mc_corr_mean, mc_corr_se = zip(*sorted(zip(mc_g, mc_corr_mean, mc_corr_se)))
    mc_corr_mean = list(mc_corr_mean)
    mc_corr_se = list(mc_corr_se)

    mc_corr_mean = list(map(lambda x: x/num_sites, mc_corr_mean))
    mc_corr_se = list(map(lambda x: x/num_sites, mc_corr_se))

    return mc_g, mc_corr_mean, mc_corr_se

dmrg_g, dmrg_corr = np.loadtxt("/Users/shaeermoeed/Github/DVRPairMC/Results/DMRG/corr.txt", skiprows=4, unpack=True)
dmrg_g_m_14, dmrg_corr_m_14 = np.loadtxt("/Users/shaeermoeed/Github/DVRPairMC/Results/DMRG/corr_M_14.txt", skiprows=4, unpack=True)
dmrg_g_m_14_dvr, dmrg_corr_m_14_dvr = np.loadtxt("/Users/shaeermoeed/Github/DVRPairMC/Results/DMRG/corr_M_14_N_150_DVR_basis_DMRG.txt", skiprows=4, unpack=True)
dmrg_g_m_5_dvr, dmrg_corr_m_5_dvr = np.loadtxt("/Users/shaeermoeed/Github/DVRPairMC/Results/DMRG/corr_N_150_M_5_DVR.txt", skiprows=4, unpack=True)

mc_results_nofit = "/Users/shaeermoeed/Github/DVRPairMC/Results/PIGS/19_07_2024_17_48_46/Parameter Sweep Cumulative No Fit.csv"
mc_g_nofit, mc_corr_mean_nofit, mc_corr_se_nofit = parse_mc_data(mc_results_nofit)

mc_results_nofit_1 = "/Users/shaeermoeed/Github/DVRPairMC/Results/PIGS/01_08_2024_18_11_46/Parameter Sweep No Fit.csv"
mc_g_nofit_1, mc_corr_mean_nofit_1, mc_corr_se_nofit_1 = parse_mc_data(mc_results_nofit_1)

mc_results_nofit_2 = "/Users/shaeermoeed/Github/DVRPairMC/Results/PIGS/06_02_2025_21_17_59/Parameter Sweep No Fit.csv"
mc_g_nofit_2, mc_corr_mean_nofit_2, mc_corr_se_nofit_2 = parse_mc_data(mc_results_nofit_2)

plt.rcParams['mathtext.fontset']='stix'
fig, ax1 = plt.subplots()
left, bottom, width, height = [0.57, 0.17, 0.3, 0.3]
ax1.plot(dmrg_g[10:], dmrg_corr[10:], label="DMRG (M=5)", color="C3")
ax1.scatter(mc_g_nofit, mc_corr_mean_nofit, label="MC (P=60, M=10)", color="C0")
ax1.errorbar(mc_g_nofit, mc_corr_mean_nofit, mc_corr_se_nofit, capsize=5, fmt="None", color='C0')
#ax1.scatter(mc_g_nofit_1, mc_corr_mean_nofit_1, label="MC (P=100, M=14)", color="black")
#ax1.errorbar(mc_g_nofit_1, mc_corr_mean_nofit_1, mc_corr_se_nofit_1, capsize=5, fmt="None", color='black')
ax1.scatter(mc_g_nofit_2, mc_corr_mean_nofit_2, label="MC (P=60, M=14)", color="C2")
ax1.errorbar(mc_g_nofit_2, mc_corr_mean_nofit_2, mc_corr_se_nofit_2, capsize=5, fmt="None", color='C2')

ax1.set_ylabel(r"$C$", fontsize=20)
ax1.set_xlabel(r"$g$", fontsize=20)
ax1.tick_params(axis='both', labelsize=12)

ax1.legend(prop={'size': 12}, loc="upper left")

plt.savefig("Correlation DMRG Comparison M 14.png", bbox_inches='tight', dpi=1200)
plt.close()

plt.figure()
plt.plot(dmrg_g[10:], dmrg_corr[10:], '.', label="DMRG (M=5)", color="C0")
plt.plot(dmrg_g_m_14[10:], dmrg_corr_m_14[10:], '--', label="DMRG (M=14)", color="C3")
plt.ylabel(r"$C$", fontsize=20)
plt.xlabel(r"$g$", fontsize=20)
plt.tick_params(axis='both', labelsize=12)
plt.legend(prop={'size': 12}, loc="upper left")
plt.savefig("Correlation DMRG Comparison (M=14 vs M=5).png", bbox_inches='tight', dpi=1200)

plt.figure()
plt.plot(dmrg_g_m_14_dvr[10:], dmrg_corr_m_14_dvr[10:], '.', label="DMRG (DVR)", color="C0")
plt.plot(dmrg_g_m_14[10:], dmrg_corr_m_14[10:], '--', label="DMRG (M basis)", color="C3")
plt.ylabel(r"$C$", fontsize=20)
plt.xlabel(r"$g$", fontsize=20)
plt.tick_params(axis='both', labelsize=12)
plt.legend(prop={'size': 12}, loc="upper left")
plt.savefig("Correlation DMRG Comparison M basis vs DVR (M=14).png", bbox_inches='tight', dpi=1200)

plt.figure()
plt.plot(dmrg_g[10:], dmrg_corr[10:], linestyle='dotted', label="DMRG (M basis, M=5)", color="C0")
plt.plot(dmrg_g_m_5_dvr[10:], dmrg_corr_m_5_dvr[10:], linestyle='dashed', label="DMRG (DVR, M=5)", color="C2")
plt.plot(dmrg_g_m_14_dvr[10:], dmrg_corr_m_14_dvr[10:], linestyle='dashdot', label="DMRG (DVR, M=14)", color="C3")
plt.plot(dmrg_g_m_14[10:], dmrg_corr_m_14[10:], linestyle='solid', label="DMRG (M basis, M=14)", color="C4")
plt.ylabel(r"$C$", fontsize=20)
plt.xlabel(r"$g$", fontsize=20)
plt.tick_params(axis='both', labelsize=12)
plt.legend(prop={'size': 12}, loc="upper left")
plt.savefig("Correlation DMRG Comparison M basis vs DVR (M=5, 14).png", bbox_inches='tight', dpi=1200)

dmrg_g_m_5_dvr_N_5, dmrg_corr_m_5_dvr_N_5 = np.loadtxt("/Users/shaeermoeed/Github/DVRPairMC/Results/DMRG/5_Rotors_M_5_DVR/corr.txt", skiprows=4, unpack=True)
dmrg_g_m_14_dvr_N_5, dmrg_corr_m_14_dvr_N_5 = np.loadtxt("/Users/shaeermoeed/Github/DVRPairMC/Results/DMRG/5_Rotors_M_14_DVR/corr.txt", skiprows=4, unpack=True)

plt.figure()
plt.plot(dmrg_g_m_5_dvr_N_5[10:], dmrg_corr_m_5_dvr_N_5[10:], linestyle='dotted', label="DMRG (DVR, M=5)", color="C0")
plt.plot(dmrg_g_m_14_dvr_N_5[10:], dmrg_corr_m_14_dvr_N_5[10:], linestyle='dashed', label="DMRG (DVR, M=14)", color="C3")
plt.title("DMRG N=5")
plt.ylabel(r"$C$", fontsize=20)
plt.xlabel(r"$g$", fontsize=20)
plt.tick_params(axis='both', labelsize=12)
plt.legend(prop={'size': 12}, loc="upper left")
plt.savefig("Correlation DMRG Comparison DVR (M=5, 14) N=5.png", bbox_inches='tight', dpi=1200)

mc_results_nofit_3 = "/Users/shaeermoeed/Github/DVRPairMC/Results/PIGS/07_02_2025_19_54_08/Parameter Sweep No Fit.csv"
mc_g_nofit_3, mc_corr_mean_nofit_3, mc_corr_se_nofit_3 = parse_mc_data(mc_results_nofit_3, 4)

mc_results_nofit_4 = "/Users/shaeermoeed/Github/DVRPairMC/Results/PIGS/07_02_2025_19_58_16/Parameter Sweep No Fit.csv"
mc_g_nofit_4, mc_corr_mean_nofit_4, mc_corr_se_nofit_4 = parse_mc_data(mc_results_nofit_4, 4)

plt.rcParams['mathtext.fontset']='stix'
fig, ax1 = plt.subplots()
left, bottom, width, height = [0.57, 0.17, 0.3, 0.3]
ax1.plot(dmrg_g_m_14_dvr_N_5[10:], dmrg_corr_m_14_dvr_N_5[10:], label="DMRG (DVR, M=14)", color="C3")
ax1.scatter(mc_g_nofit_3, mc_corr_mean_nofit_3, label="MC (P=60, M=14)", color="C0")
ax1.errorbar(mc_g_nofit_3, mc_corr_mean_nofit_3, mc_corr_se_nofit_3, capsize=5, fmt="None", color='C0')
ax1.scatter(mc_g_nofit_4, mc_corr_mean_nofit_4, label="MC (P=60, M=10)", color="black")
ax1.errorbar(mc_g_nofit_4, mc_corr_mean_nofit_4, mc_corr_se_nofit_4, capsize=5, fmt="None", color='black')
plt.plot(dmrg_g_m_5_dvr_N_5[10:], dmrg_corr_m_5_dvr_N_5[10:], linestyle='dotted', label="DMRG (DVR, M=5)", color="C0")
ax1.set_ylabel(r"$C$", fontsize=20)
ax1.set_xlabel(r"$g$", fontsize=20)
ax1.tick_params(axis='both', labelsize=12)
ax1.legend(prop={'size': 12}, loc="upper left")
plt.savefig("Correlation DMRG DVR PIGS DVR Comparison N 5.png", bbox_inches='tight', dpi=1200)
plt.close()
