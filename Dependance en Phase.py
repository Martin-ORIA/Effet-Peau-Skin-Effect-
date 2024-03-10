import numpy as np 
import matplotlib.pyplot as plt 
from scipy.optimize import curve_fit
import cmath

#e = 0.51 mm
L_f = np.array([30, 40, 50, 60, 70, 80, 90, 100, 150, 200, 300, 400, 500, 600, 700, 800, 900, 1000])*10**3

L_V1 = np.array([200, 198, 198, 196, 196, 194, 194, 194, 190, 188, 182, 174, 168, 160, 152, 146, 138, 130])*10**(-3)

L_V2 = np.array([10.3, 12, 13.2, 14.1, 14.8, 15.3, 15.6, 16, 16.7, 16.8, 16.5, 15.7, 14.7, 13.7, 12.7, 11.9, 11.2, 10.3])*10**(-3)
L_Phase = np.array([5.68, 3.27, 2.46, 1.78, 1.26, 0.92, 0.7, 0.57, 0.15, 0, 0.1, 0.14, 0.16, 0.17,  0.18, 0.19, 0.18, 0.18])*10**(-6)




#e = 1.22 mm
L_f_07 = np.array([20, 30, 40, 50, 60, 70, 80, 90, 100, 150, 200, 300, 400, 500, 600,  700, 800, 900, 1000])*10**3
L_V1_07 = np.array([202, 200, 198, 198, 198, 196, 194, 194, 194, 190, 188, 182, 174, 168, 160, 152, 146, 138, 130])*10**(-3)
L_V2_07 = np.array([5.7, 6.56, 7.12, 7.3, 7.55, 7.65, 7.68, 7.7, 7.76, 7.52, 7.2, 6.48, 5.76, 5, 4.4, 3.84, 3.3, 3, 2.5])*10**(-3)
L_Phase_07 = np.array([6, 2.88, 1.56, 0.84, 0.56, 0.2, 0.02, 0.01, 0.2, 0.44, 0.48, 0.52, 0.51, 0.5, 0.5, 0.48, 0.47, 0.45, 0.45])*10**(-6)


#à vide
L_V1_00 = np.array([200, 200, 198, 198, 196, 196, 194, 194, 190, 188, 182, 174, 166, 158, 150, 142, 134, 128])*10**(-3)
L_V2_00 = np.array([13.2, 17.8, 21.8, 25.8, 29.6, 33.6, 37.4, 40.8, 60, 76.8, 108, 136, 156, 172, 184, 194, 198, 204])*10**(-3)
L_Phase_00 = np.array([8, 6.64, 4.8, 3.9, 3.4, 2.8, 2.5, 2.27, 1.46, 1.43, 0.66, 0.79, 0.66, 0.56, 0.5, 0.46, 0.42, 0.39])*10**(-9)




delta = np.sqrt(2/(1.256*10**(-6)*1303594*2*np.pi*L_f))
L = 2*np.pi*L_f*np.exp(-0.51*10**(-3)/delta)*np.exp(L_Phase*L_f*2*np.pi)

delta_07 = np.sqrt(2/(1.256*10**(-6)*1303594*2*np.pi*L_f_07))
L_07 = 2*np.pi*L_f_07*np.exp(-1.22*10**(-3)/delta_07)*np.exp(L_Phase_07*L_f_07*2*np.pi)

delta_00 = np.sqrt(2/(1.256*10**(-6)*1303594*2*np.pi*L_f))
L_00 = 2*np.pi*L_f*np.exp(-0/delta_00)*np.exp(L_Phase_00*L_f*2*np.pi)



L_Phase_07_inc = np.array([0.5, 0.25, 0.2, 0.15, 0.1, 0.05, 0.05, 0.005, 0.05, 0.08, 0.08, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1])*10**(-6)
L_Phase_inc = np.array([0.5, 0.25, 0.2, 0.15, 0.1, 0.1, 0.1, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05])*10**(-6)



L_rapport_07_inc = (L_V2_07/L_V1_07)*np.sqrt( ((0.3*10**(-3))/(L_V2_07))**2 + ((1*10**(-3))/(L_V1_07))**2 )
L_rapport_inc = (L_V2/L_V1)*np.sqrt( ((0.3*10**(-3))/(L_V2))**2 + ((1*10**(-3))/(L_V1))**2 )


L_f_07_inc = np.ones(len(L_f_07))*0.01*10**3
L_f_inc = np.ones(len(L_f))*0.01*10**3


L = ((L_V2/L_V1)[0]/L[0])*L
L_07 = ((L_V2_07/L_V1_07)[0]/L_07[0])*L_07
L_00 = ((L_V2_00/L_V1_00)[0]/L_00[0])*L_00







def fiting_method(X, a, delta, e):
	x, y = X

	return (2j*np.pi*a*x*np.exp(-1.22*10**(-3) * delta*np.sqrt(2*np.pi*x) ) * np.exp(1j*(2*np.pi*x*y -(1.22*10**(-3))*np.sqrt(2*np.pi*x)*delta )) +e).real



aa = -1
opt, cov = curve_fit(fiting_method, (L_f_07[:aa], L_Phase_07[:aa]), (L_V2_07/L_V1_07)[:aa], maxfev=10000, sigma=L_rapport_07_inc[:aa], absolute_sigma=True)
a, delta, e = opt
mu_theorique = 1.256*10**(-6)
mu_exp = (delta**2)*2/1303594
delta_inc = np.sqrt(np.diag(cov))[1]

print("#############################")
print("fitting pour e=1.22mm mu0 = ", mu_exp,"+- ", mu_exp*np.sqrt(2*((delta_inc)/(delta))**2))
print("Ecart relatif = ", np.abs((mu_exp-mu_theorique)/(mu_theorique))*100)




f_indexes = np.linspace(0, len(L_f_07) - 1, 1000)
L_f_fit = np.interp(f_indexes, np.arange(len(L_f_07)), L_f_07)

Phase_indexes = np.linspace(0, len(L_Phase_07) - 1, 1000)
L_Phase_fit = np.interp(Phase_indexes, np.arange(len(L_Phase_07)), L_Phase_07)



fig, ax = plt.subplots()
ax.scatter(L_f_07, (L_V2_07/L_V1_07), s=10, color="black", label="e=1.22mm")
ax.plot(L_f_fit, fiting_method((L_f_fit, L_Phase_fit), a, delta, e), linewidth=1, color="r", label="FIT b={}".format(round(delta, 5)))
ax.errorbar(L_f_07, (L_V2_07/L_V1_07), xerr=L_f_07_inc, yerr=L_rapport_07_inc, elinewidth=1, ecolor="black", fmt='None', capsize=2)

plt.xlabel("Frequence (Hz)")
plt.ylabel("Coefficient de transmission")
plt.legend()
plt.grid()


#POUR AFFICHER LA PHASE
# fig, ax = plt.subplots()
# ax.plot(L_Phase_fit, fiting_method((L_f_fit, L_Phase_fit), a, delta, e), linewidth=2, color="g", label="FIT")
# ax.errorbar(L_Phase_07, (L_V2_07/L_V1_07), xerr=L_Phase_07_inc, yerr=L_rapport_07_inc, elinewidth=1, ecolor="black", fmt='None', capsize=2)
# ax.scatter(L_Phase_07, (L_V2_07/L_V1_07), s=10, color="black")
# plt.xlabel("Frequence (Hz)")
# plt.ylabel("Coefficient de transmission")
# plt.grid()

a, delta, e = 0, 0, 0
#####################
print("#############################")
def fiting_method2(X, a, delta, e):
	x, y = X

	return (2j*np.pi*a*x*np.exp(-0.51*10**(-3) * delta*np.sqrt(2*np.pi*x) ) * np.exp(1j*(2*np.pi*x*y -(0.51*10**(-3))*np.sqrt(2*np.pi*x)*delta )) +e).real


p_0 = 1
aa = -1
opt, cov = curve_fit(fiting_method2, (L_f[:aa], L_Phase[:aa]), (L_V2/L_V1)[:aa], maxfev=10000, sigma=L_rapport_inc[:aa], absolute_sigma=True)
a, delta, e = opt

mu_theorique = 1.256*10**(-6)
mu_exp = (delta**2)*2/1753171
delta_inc = np.sqrt(np.diag(cov))[1]

print("fitting pour e=0.51mm mu0 = ",mu_exp,  "+- ", mu_exp*np.sqrt(2*((delta_inc)/(delta))**2))
print("Ecart relatif = ", np.abs((mu_exp-mu_theorique)/(mu_theorique))*100)
print("#############################")

f_indexes = np.linspace(0, len(L_f) - 1, 1000)
L_f_fit = np.interp(f_indexes, np.arange(len(L_f)), L_f)

Phase_indexes = np.linspace(0, len(L_Phase) - 1, 1000)
L_Phase_fit = np.interp(Phase_indexes, np.arange(len(L_Phase)), L_Phase)

fig, ax = plt.subplots()
ax.scatter(L_f, (L_V2/L_V1), s=10, color="black", label="e=0.51mm")
ax.plot(L_f_fit, fiting_method2((L_f_fit, L_Phase_fit), a, delta, e), linewidth=1, color="blue", label="FIT b={}".format(round(delta,5)))
ax.errorbar(L_f, (L_V2/L_V1), xerr=L_f_inc, yerr=L_rapport_inc, elinewidth=1, ecolor="black", fmt='None', capsize=2)

plt.xlabel("Frequence (Hz)")
plt.ylabel("Coefficient de transmission")
plt.legend()
plt.grid()


#POUR AFFICHER LA PHASE
# fig, ax = plt.subplots()
# ax.plot(L_Phase_fit, fiting_method((L_f_fit, L_Phase_fit), a, delta, e), linewidth=2, color="g", label="FIT")
# ax.errorbar(L_Phase, (L_V2/L_V1), xerr=L_Phase_inc, yerr=L_rapport_inc, elinewidth=1, ecolor="black", fmt='None', capsize=2)
# ax.scatter(L_Phase, (L_V2/L_V1), s=10, color="black")
# plt.xlabel("Frequence (Hz)")
# plt.ylabel("Coefficient de transmission")
# plt.grid()
# plt.show()


fig2, ax2 = plt.subplots()
ax2.errorbar(L_f_07, L_Phase_07*L_f_07*2*np.pi, xerr=L_f_07_inc, yerr=L_Phase_07_inc*L_f_07*2*np.pi, elinewidth=0.5, ecolor="black", fmt='None', capsize=2)
ax2.scatter(L_f_07, L_Phase_07*L_f_07*2*np.pi, label="e=1.2mm", c="b", s=15)
ax2.errorbar(L_f, L_Phase*L_f*2*np.pi, xerr=L_f_inc, yerr=L_Phase_inc*L_f*2*np.pi, elinewidth=0.5, ecolor="black", fmt='None', capsize=2)
ax2.scatter(L_f, L_Phase*L_f*2*np.pi, label="e=0.51mm", c="r", s=15)
ax2.errorbar(L_f, L_Phase_00*L_f*2*np.pi, xerr=L_f_inc, yerr=L_Phase_inc*L_f*2*np.pi, elinewidth=0.5, ecolor="black", fmt='None', capsize=2)
ax2.scatter(L_f, L_Phase_00*L_f*2*np.pi, label="à vide", c='g', s=15)
#ax2.axhline(-np.pi/2, color="r", linestyle="--")
ax2.set_xlabel("frequence (Hz)")
ax2.set_ylabel("Dephasage (rad)")
plt.legend()
plt.grid()
#ax.scatter(L_f, L_V2)
plt.show()