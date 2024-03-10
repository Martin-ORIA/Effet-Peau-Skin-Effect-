import numpy as np 
import matplotlib.pyplot as plt 
from scipy.optimize import curve_fit

L_f = np.array([30, 40, 50, 60, 70, 80, 90, 100, 150, 200, 300, 400, 500, 600, 700, 800, 900, 1000])*10**(-3)
L_V1 = np.array([200, 198, 198, 196, 196, 194, 194, 194, 190, 188, 182, 174, 168, 160, 152, 146, 138, 130])*10**(-3)
L_V2 = np.array([10.3, 12, 13.2, 14.1, 14.8, 15.3, 15.6, 16, 16.7, 16.8, 16.5, 15.7, 14.7, 13.7, 12.7, 11.9, 11.2, 10.3])*10**(-3)

L_V1_00 = np.array([200, 200, 198, 198, 196, 196, 194, 194, 190, 188, 182, 174, 166, 158, 150, 142, 134, 128])*10**(-3)
L_V2_00 = np.array([13.2, 17.8, 21.8, 25.8, 29.6, 33.6, 37.4, 40.8, 60, 76.8, 108, 136, 156, 172, 184, 194, 198, 204])*10**(-3)
L_Phase_00 = np.array([8, 6.64, 4.8, 3.9, 3.4, 2.8, 2.5, 2.27, 1.46, 1.43, 0.66, 0.79, 0.66, 0.56, 0.5, 0.46, 0.42, 0.39])*10**(-6)


L_f_07 = np.array([30, 40, 50, 60, 70, 80, 90, 100, 150, 200, 300, 400, 500, 600,  700, 800, 900, 1000])*10**(-3)
L_V1_07 = np.array([200, 198, 198, 198, 196, 194, 194, 194, 190, 188, 182, 174, 168, 160, 152, 146, 138, 130])*10**(-3)
L_V2_07 = np.array([6.56, 7.12, 7.3, 7.55, 7.65, 7.68, 7.7, 7.76, 7.52, 7.2, 6.48, 5.76, 5, 4.4, 3.84, 3.3, 3, 2.5])*10**(-3)

L_f_inc = np.ones(len(L_f))*10**(3)*0.01
L_V_inc = np.ones(len(L_f))*1*10**(-3)
L_rapp_inc = (L_V2/L_V2_00)*np.sqrt( ( (L_V_inc)/(L_V2) )**2 +( (L_V_inc)/(L_V2_00) )**2 )
L_rapp2_inc = (L_V2_07/L_V2_00)*np.sqrt( ( (L_V_inc)/(L_V2_07) )**2 +( (L_V_inc)/(L_V2_00) )**2 )


x1 = 0
x2 = len(L_f)


def transmission(x, alpha, a, b):
	return a*np.exp(-alpha*np.sqrt(x*2*np.pi))+b




alpha_theo_five = 0.51*10**(-3)*np.sqrt(1.256*10**(-6)*1303594)/2
alpha_theo_seven = 1.22*10**(-3)*np.sqrt(1.256*10**(-6)*1303594)/2


opt, cov = curve_fit(transmission, L_f[x1:x2], (L_V2/L_V2_00)[x1:x2], sigma=L_rapp_inc, absolute_sigma=True)
opt2, cov2 = curve_fit(transmission, L_f_07[x1:x2], (L_V2_07/L_V2_00)[x1:x2], sigma=L_rapp2_inc, absolute_sigma=True)


alpha, a, b = opt
alpha2, a2, b2 = opt2

#print(alpha, np.sqrt(np.diag(cov))[0])
#print(alpha2, np.sqrt(np.diag(cov2))[0])

alpha_inc = np.sqrt(np.diag(cov))[0]
alpha2_inc = np.sqrt(np.diag(cov2))[0]
mu0 = (alpha*2/(0.51*np.sqrt(1753171) ))**2
mu02 = (alpha2*2/(1.22*np.sqrt(1303594) ))**2
mu0_the = 1.256*10**(-6)

print("mu0 = ", mu0, "  ", mu0*np.sqrt(2*( (alpha_inc)/(alpha) )**2 ))
print("mu0 2 = ", mu02, "  ", mu02*np.sqrt(2*( (alpha2_inc)/(alpha2) )**2 ))


fig, ax = plt.subplots()

ax.errorbar(L_f[x1:x2]*10**(3), (L_V2/L_V2_00)[x1:x2], xerr=L_f_inc[x1:x2], yerr=L_rapp_inc[x1:x2], elinewidth=1, ecolor="black", fmt='None')
ax.scatter(L_f[x1:x2]*10**(3), (L_V2/L_V2_00)[x1:x2], label="e=0.51mm", color="red", s=10)
ax.plot(np.linspace(min(L_f[x1:x2]), max(L_f[x1:x2]), 100)*10**(3), transmission(np.linspace(min(L_f[x1:x2]), max(L_f[x1:x2]), 100), alpha,a, b), color="red", linestyle="--", label="FIT : alpha={} SI".format(round(alpha*10**(-3), 4)))


ax.errorbar(L_f[x1:x2]*10**(3), (L_V2_07/L_V2_00)[x1:x2], xerr=L_f_inc[x1:x2], yerr=L_rapp2_inc[x1:x2], elinewidth=1, ecolor="black", fmt='None')
ax.scatter(L_f[x1:x2]*10**(3), (L_V2_07/L_V2_00)[x1:x2], label="e=1.22mm", color="mediumblue", s=10)
ax.plot(np.linspace(min(L_f_07[x1:x2]), max(L_f_07[x1:x2]), 100)*10**(3), transmission(np.linspace(min(L_f_07[x1:x2]), max(L_f_07[x1:x2]), 100), alpha2,a2, b2), color="navy", linestyle="--", label="FIT : alpha={} SI".format(round(alpha2*10**(-3), 4)))

ax.set_title("Coefficient de transmission")
ax.set_xlabel("frequence kHz")
ax.set_ylabel("Rapport de tension V2/V0")

plt.legend()
plt.grid()
plt.show()