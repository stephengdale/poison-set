import numpy as np
import yaml
import psi4 
import matplotlib.pyplot as plt
import scipy.optimize as opt
from scipy.optimize import fsolve


test_set = "P30-5"

reactions = yaml.load(open(f"{test_set}/Reactions.yaml"),
                      Loader=yaml.UnsafeLoader)

psi4.set_options({'basis': 'def2-svp',
                  'scf_type': 'df',
                  'reference' : 'uhf',
                  'save_jk' : True,
                  'soscf' : False,
                  'df_ints_io' : 'save',
                  'e_convergence': 1e-8,
                  'd_convergence': 1e-8}) 

def single_calc(mol, alpha = 0.0):
    custom = {"name" : "global_hybrid",
              "x_functionals" : {"GGA_X_PBE" : {"alpha" : 1.0-alpha}},
              "x_hf" : {"alpha" : 0.0+alpha},
              "c_functionals" : {"GGA_C_PBE" : {}}
              }

    psi4.set_output_file(f"{mol[:-4]}.out")
    xyz = open("P30-5/" + mol, 'r')
    xyz = "".join(xyz.readlines()[1:])
    xyz = psi4.geometry(xyz)
    return psi4.energy('scf', dft_functional=custom, molecule=xyz)

def full_calc(alpha = 0.0):
    data = {}
    for i in reactions[r][1:]:
        data[i[1]] = single_calc(i[1], alpha) 
    target = 0
    for i in reactions[r][1:]:
        target += i[0]*data[i[1]]
    return target*627.5

for r in reactions:
    check = open(f'{test_set}.txt', 'r')
    ref = reactions[r][0]
    skip = 0
    for line in check:
        if line.split()[0] == r:
            skip = 1
    check.close()
    if skip == 0:
        keep = open(f'{test_set}.txt', 'a')
        keep.write(f'{r:10s} {full_calc(0.0)-ref:10.3f} {full_calc(0.05)-ref:10.3f} {full_calc(0.1)-ref:10.3f} {full_calc(0.15)-ref:10.3f} {full_calc(0.2)-ref:10.3f} {full_calc(0.4)-ref:10.3f} {full_calc(0.6)-ref:10.3f} {full_calc(0.8)-ref:10.3f} {full_calc(1.0)-ref:10.3f}\n')
        keep.close()

# x = np.arange(-0.5, 1.1, 0.1)
# starting_guess = 0.
# 
# def Quadratic(x, a, b, c):
#     return a*x**2 + b*x + c
# 
# read = open('alphaGraph.txt', 'r')
# keep = open('alphaExtrapolated.txt', 'w')
# for r in read.readlines():
#     print(r)
#     read2 = open('alphaGraph2.txt', 'r')
#     for r2 in read2.readlines():
#         if r2.split()[0] == r.split()[0]:
#             dat = [float(i) for i in r2.split()[1:]]
# 
#             Params, Pcov = opt.curve_fit(Quadratic, [0., 0.05, 0.10, 0.15, 0.20], dat, maxfev=5000)
#             a, b, c = tuple(Params)
#             QFit = Quadratic(x, a, b, c)
#             def Q(x):
#                 return a*x**2 + b*x + c - reactions[r.split()[0]][0]
#             print(a, b, c, reactions[r.split()[0]][0], fsolve(Q, starting_guess))
#             alpha2 = fsolve(Q, starting_guess)[0]
# 
#     read3 = open('alphaGraph3.txt', 'r')
#     for r3 in read3.readlines():
#         if r3.split()[0] == r.split()[0]:
#             dat = [float(i) for i in r3.split()[1:]]
# 
#             Params, Pcov = opt.curve_fit(Quadratic, [0., 0.05, 0.10, 0.15, 0.20, 0.40, 0.60, 0.80, 1.0], dat, maxfev=5000)
#             a, b, c = tuple(Params)
#             QFit = Quadratic(x, a, b, c)
#             def Q(x):
#                 return a*x**2 + b*x + c - reactions[r.split()[0]][0]
#             print(a, b, c, reactions[r.split()[0]][0], fsolve(Q, starting_guess))
#             alpha3 = fsolve(Q, starting_guess)[0]
# 
#     dat = [float(i) for i in r.split()[1:]]
# 
#     Params, Pcov = opt.curve_fit(Quadratic, [0., 0.4, 0.6, 0.8, 1.0], dat, maxfev=5000)
#     a, b, c = tuple(Params)
#     QFit = Quadratic(x, a, b, c)
#     def Q(x):
#         return a*x**2 + b*x + c - reactions[r.split()[0]][0]
#     print(a, b, c, reactions[r.split()[0]][0], fsolve(Q, starting_guess))
#     alpha1 = fsolve(Q, starting_guess)[0]
#     keep.write(f'{r.split()[0]:10s} {alpha1*100:10.3f} {alpha2*100:10.3f} {alpha3*100:10.3f}\n')
# 
#     fig, ax = plt.subplots()
# 
#     ax.scatter([0., 0.4, 0.6, 0.8, 1.0], dat, s=30)
#     ax.axhline(y=reactions[r.split()[0]][0], color='r')
#     ax.plot(x, QFit, 'b')
# 
#     ax.set_title(r.split()[0])
#     ax.set_xlabel('alpha')
#     ax.set_ylabel('delta E (kcal/mol)')
# 
#     plt.tight_layout()
#     plt.savefig(r.split()[0] + '.png')
#     plt.close()
#     print(r.split()[0] + '.png')



