import numpy as np 
import matplotlib.pyplot as plt 
from scipy.optimize import curve_fit


def R_parall(resistances):
	reciprocal_sum = sum(1 / R for R in resistances)
	R_total = 1 / reciprocal_sum
	return R_total

#   03      07         06       05      08
#[1248673, 1303595, 1298975, 1753171, 1294474]
# 07  05  03  06  08  06+07  07+05  07+05+06+08+03  05+06+07  05+03

L_sigma = np.array([1303595, 1753171, 1248673, 1298975, 1294474, np.mean((1298975, 1303595)),np.mean((1303595, 1753171)), np.mean((1753171,1298975,1303595)), np.mean((1753171, 1248673))])
#L_sigma = np.array([1303595, 1753171, 1248673, 1298975, 1294474, R_parall((1298975, 1303595)),R_parall((1303595, 1753171)), R_parall((1753171,1298975,1303595)), R_parall((1753171, 1248673))])
w =  2*np.pi*300*10**(3)


T = np.log(np.array([2.4, 6.2, 4.8, 3.2, 2, 0.8, 1.4, 0.5, 2.8])/50.6)/(np.sqrt(L_sigma*w))
T2 = np.log(np.array([2.4, 6.2, 4.8, 3.2, 2, 0.8, 1.4, 0.5, 2.8])/50.6)/(np.sqrt(1303595*w))

inc_rapport = (np.array([2.4, 6.2, 4.8, 3.2, 2, 0.8, 1.4, 0.5, 2.8])/50.6)*np.sqrt( ((0.1)/(np.array([2.4, 6.2, 4.8, 3.2, 2, 0.8, 1.4, 0.5, 2.8])))**2 + ((0.1)/(50.6))**2 )
T_inc = ((50.6/np.array([2.4, 6.2, 4.8, 3.2, 2, 0.8, 1.4, 0.5, 2.8]))*inc_rapport)/(np.sqrt(L_sigma*w))
T_inc_2 = ((50.6/np.array([2.4, 6.2, 4.8, 3.2, 2, 0.8, 1.4, 0.5, 2.8]))*inc_rapport)/(np.sqrt(1303595*w))

L_e = np.array([1.22, 0.49, 0.69, 1.03, 1.48, 2.25, 1.71, 2.74, 1.18])*10**(-3)
L_e_inc = np.ones(len(L_e))*0.05*10**(-3)

fit_e = np.linspace(min(L_e), max(L_e), 100)


def modele(x, a, b):
	return a*x+b 



opt, cov = curve_fit(modele, L_e, T)
a_opt, b_opt = opt


opt2, cov2 = curve_fit(modele, L_e, T2)
a_opt2, b_opt2 = opt2


mu_theorique = 1.256*10**(-6)
mu_exp = ((a_opt*2)**2)
mu_exp2 = ((a_opt2*2)**2)


print("Ecart relatif variable = ", np.abs((mu_exp-mu_theorique)/(mu_theorique))*100)
print("mu0 = ", mu_exp, "+-", mu_exp*np.sqrt(2*((np.sqrt(np.diag(cov))[0])/(a_opt))**2 ))

print("Ecart relatif uniforme = ", np.abs((mu_exp2-mu_theorique)/(mu_theorique))*100)
print("mu0 = ", mu_exp2, "+-", mu_exp2*np.sqrt(2*((np.sqrt(np.diag(cov2))[0])/(a_opt2))**2 ))


plt.errorbar(L_e*10**3, T, xerr=L_e_inc, yerr=T_inc, elinewidth=1, ecolor="black", fmt='None')
plt.scatter(L_e*10**3, T, color="b", label="Mesure (σ variable)", s=20)
plt.plot(fit_e*10**3, modele(fit_e, opt[0], opt[1]), color="blue", linestyle="--", label="FIT y={}z+{}".format(round(a_opt, 5), round(b_opt, 5)))

plt.errorbar(L_e*10**3, T2, xerr=L_e_inc, yerr=T_inc_2, elinewidth=1, ecolor="red", fmt='None')
plt.scatter(L_e*10**3, T2, color="r", label="Mesure (σ uniforme)", s=10, alpha=0.75)
plt.plot(fit_e*10**3, modele(fit_e, opt2[0], opt2[1]), color="r", linestyle="--", label="FIT y={}z+{}".format(round(a_opt2, 5), round(b_opt2, 5)))
plt.plot(fit_e*10**(3), -fit_e*np.sqrt(1.256*10**(-6))/2 + opt[1], color="black", label="Theorie")
plt.xlabel("Epaisseur (mm)")
plt.ylabel("log(t)/sqrt(σw)  (SI)")
plt.legend()
plt.grid()
#plt.yscale('function')
plt.show()
