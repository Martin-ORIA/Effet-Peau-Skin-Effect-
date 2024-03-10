import numpy as np 
import matplotlib.pyplot as plt 
from scipy.optimize import curve_fit
from colour import Color
red = Color("lightcoral")
blue = Color("lightblue")
colors1 = list(red.range_to(Color("darkred"),3))
colors2 = list(blue.range_to(Color("midnightblue"),3))

colors = colors1 + colors2

I_00 = np.array([0.2, 0.4, 0.59, 0.79, 1.00, 1.20, 1.39, 1.61, 1.80, 2.00, 2.20, 2.40, 2.60, 2.81, 3.04])
V_00 = np.array([22, 50, 78, 103, 129, 156, 181, 209, 235, 262, 287, 315, 340, 366, 396])*10**(-6)

I_00 = np.array([0.2, 0.59, 1.00, 2.20, 2.60, 3.04])
V_00 = np.array([22, 78, 129, 287, 340, 396])*10**(-6)

#0.69 mm
I_03 = np.array([0.25, 0.51, 1.01, 2.51, 2.75, 3.04])
V_03 = np.array([0.032, 0.068, 0.132, 0.322, 0.354, 0.391])*10**(-3)

#1.22 mm
I_07 = np.array([0.25, 0.51, 1.00, 2.50, 2.75, 3.04])
V_07 = np.array([0.0186, 0.038, 0.072, 0.175, 0.193, 0.213])*10**(-3)

#1.03 mm
I_06 = np.array([0.26, 0.50, 1.01, 2.51, 2.75, 3.04])
V_06 = np.array([0.023, 0.044, 0.086, 0.209, 0.229, 0.253])*10**(-3)

#0.49 mm
I_05 = np.array([0.26, 0.50, 1.01, 2.49, 2.75, 3.04])
V_05 = np.array([0.033, 0.065, 0.130, 0.320, 0.354, 0.390])*10**(-3)

#1.48 mm
I_08 = np.array([0.24, 0.50, 1.01, 2.51, 2.76, 3.04])
V_08 = np.array([0.017, 0.032, 0.061, 0.148, 0.162, 0.178])*10**(-3)

inc_I = np.ones(len(I_08))*0.02
inc_V = np.ones(len(V_08))*5*10**(-6)


L_e = [0.49, 0.69, 0.7, 1.03, 1.22, 1.48]
I = [I_05, I_03, I_00, I_06, I_07, I_08]
V = [V_05, V_03, V_00, V_06, V_07, V_08]
L_color=["y", "b", "g", "r", "c", "m"]
L_label=["e=0.49 mm", "e=0.69 mm", "e=0.70 mm", "e=1.03 mm", "e=1.22 mm", "e=1.48 mm"]


def modele(x, a, b):
    return a*x + b

L_R = []
L_R_inc = []

def fitting_plotting(I, V, e, lab, col):

    opt, p_cov = curve_fit(modele, I, V, sigma=inc_V, absolute_sigma=True)

    a_opt, b_opt = opt
    fit_I = np.linspace(min(I), max(I), 100)
    a_opt_inc = np.sqrt(np.diag(p_cov))[0]
    L_R.append(1/((1/a_opt)*((np.log(2))/(2*e*10**(-3)*np.pi))))
    L_R_inc.append(1/((1/a_opt_inc)*((np.log(2))/(2*e*10**(-3)*np.pi))))

    print("e = ", e)
    print("sigma =", (1/np.array((1/((1/a_opt)*((np.log(2))/(2*e*10**(-3)*np.pi))))))*1, "+-", (np.array((1/((1/a_opt_inc)*((np.log(2))/(2*e*10**(-3)*np.pi)))))/(np.array((1/((1/a_opt)*((np.log(2))/(2*e*10**(-3)*np.pi)))))**2)), "SI")
    print("#################")

    plt.scatter(I, V, color=col, s=10, label=lab)
    plt.errorbar(I, V, xerr=inc_I, yerr=inc_V, elinewidth=1, ecolor=col, fmt='None')
    plt.plot(fit_I, a_opt*fit_I+b_opt, color=col, linestyle="--")

    return 0

for i in range(len(I)):
    fitting_plotting(I[i], V[i], L_e[i], L_label[i], str(colors[i]))
print(type(colors[1]))
plt.legend()
plt.grid()
plt.xlabel("I (A)")
plt.ylabel("V (MicroV)")
plt.title("Détermination de la conductivité")
plt.show()