import numpy as np
import matplotlib.pyplot as plt
import os
from scipy.optimize import curve_fit 

def func_quad(x, a, b):

    return a*x*x + b

def func_cubic(x, a, b):

    return a*x*x*x + b

def fit_error(x_pts, y_pts, fit_func, num_pts_out=100):

    coefficients, cov = curve_fit(fit_func, x_pts, y_pts)
    err = np.sqrt(np.diag(cov)[1])

    fit_x = []
    fit_y = []
    for i in range(num_pts_out):
        x = max(x_pts) * (1.0 - (float(i)/float(num_pts_out-1.0)))
        y = fit_func(x, coefficients[0], coefficients[1])
        fit_x.append(x)
        fit_y.append(y)

    return fit_x, fit_y, err

primitive_results = "/Users/shaeermoeed/Github/DVRPairMC/primitive_gibbs_run.csv"
pair_results = "/Users/shaeermoeed/Github/DVRPairMC/pair_gibbs_run.csv"

g,p,tau_prim,e0_prim,e0_err_prim = np.loadtxt(primitive_results, skiprows=4, unpack=True, delimiter=",")
pair_data = np.loadtxt(pair_results, skiprows=2, unpack=False, delimiter=",")
tau_pair = pair_data[:, 3]
e0_pair = pair_data[:, -2]
e0_pair_err = pair_data[:, -1]

print(tau_pair)
print(e0_pair)

tau_pair, e0_pair, e0_pair_err = zip(*sorted(zip(tau_pair, e0_pair, e0_pair_err)))

print(tau_pair)
print(e0_pair)

tau_pair = tau_pair[:-2]
e0_pair = e0_pair[:-2]
e0_pair_err = e0_pair_err[:-2]
tau_prim = tau_prim[2:]
e0_prim = e0_prim[2:]
e0_err_prim = e0_err_prim[2:]

tau_fit_pair, e0_fit_pair, e0_fit_err_pair = fit_error(tau_pair, e0_pair, func_quad)
tau_fit_prim, e0_fit_prim, e0_fit_err_prim = fit_error(tau_prim, e0_prim, func_quad)

plt.figure()
plt.rcParams['mathtext.fontset']='stix'
plt.scatter(tau_pair, e0_pair, label="Pair", color="C0")
plt.errorbar(tau_pair, e0_pair, e0_pair_err, fmt="None", capsize=5, color="C0")
plt.plot(tau_fit_pair, e0_fit_pair, color="C0")
#plt.scatter([0.0], e0_fit_pair[-1], color="black")
plt.errorbar([0.0], e0_fit_pair[-1], e0_fit_err_pair, capsize=5, fmt="None", color="black")
plt.scatter(tau_prim, e0_prim, label="Primitive", color="C3")
plt.errorbar(tau_prim, e0_prim, e0_err_prim, fmt="None", capsize=5, color="C3")
plt.plot(tau_fit_prim, e0_fit_prim, color="C3")
plt.errorbar([0.0], e0_fit_prim[-1], e0_fit_err_prim, capsize=5, color="black", fmt="None")
#plt.scatter([0.0], e0_fit_prim[-1], color="black")
plt.plot(tau_fit_pair, [-73.39831131744273]*100, label="DMRG", color="black")
plt.xlabel(r"$\tau$", fontsize=20)
plt.ylabel(r"$E_0$", fontsize=20)
plt.xticks(fontsize=12)
plt.yticks(fontsize=12)
plt.legend(prop={'size': 12})
plt.savefig("MC_Pair_Primitive_Comparison.png", bbox_inches='tight', dpi=1200)
plt.show()