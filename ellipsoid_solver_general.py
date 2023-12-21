import numpy as np
import scipy.optimize as optimize
import time

mu = 1/10
n = 6
r = [np.sqrt(mu)/(1+np.sqrt(mu)), -np.sqrt(mu)/(1+np.sqrt(mu)), np.sqrt(mu)/(1-np.sqrt(mu)), -np.sqrt(mu)/(1-np.sqrt(mu))]

#rmatrix = np.array([[0, r[0], r[0], r[0], r[0]],
#                 [r[0], 0, r[0], r[0], r[1]],
##                 [r[0], r[0], 0, r[1], r[0]],
 #                [r[0], r[0], r[1], 0, r[0]],
 #                [r[0], r[1], r[0], r[0], 0]])
# Armin's numerical solution
rmatrix = np.array([[0,1/2,1/4,-1/2,1/4], [1/2,0,1/2,1/2,1/2], [1/4,1/2,0,1/4,-1/2], [-1/2,1/2,1/4,0,1/4], [1/4,1/2,-1/2,1/4,0]])

def find_b(rmat):
    b = np.empty((5,5))
    for i in range(len(b)):
        for j in range(len(b)):
            if (i != j) & (i > 0) & (j > 0):
                b[i][j] = (rmat[i][j]-rmat[i][0]*rmat[j][0])/(np.sqrt(1-rmat[i][0]*rmat[i][0])*np.sqrt(1-rmat[j][0]*rmat[j][0]))
            else:
                b[i][j] = 0

    return b

def find_arr(b):
    arr = [1/(1+b[1][2]), 1/(1-b[1][2]), 1/(1+b[1][3]), 1/(1-b[1][3]), 1/(1+b[1][4]), 1/(1-b[1][4]), 1/(1+b[2][3]), 1/(1-b[2][3]), 1/(1+b[2][4]), 1/(1-b[2][4]), 1/(1+b[3][4]), 1/(1-b[3][4])]

    if np.abs(np.round(b[1][2], 4)) == 1:
        arr[0], arr[1] = 0, 0
    if np.abs(np.round(b[1][3], 4)) == 1:
        arr[2], arr[3] = 0, 0
    if np.abs(np.round(b[1][4], 4)) == 1:
        arr[4], arr[5] = 0, 0
    if np.abs(np.round(b[2][3],4)) == 1:
        arr[6], arr[7] = 0, 0
    if np.abs(np.round(b[2][4],4)) == 1:
        arr[8], arr[9] = 0, 0
    if np.abs(np.round(b[3][4],4)) == 1:
        arr[10], arr[11] = 0, 0

    return arr
# (non transformed) ellipse condition + 2 leftover conditions + conditions imposed when beta = \pm 1
def ellipse(z, arr, neq):
    # z = [z23, z24, z25, z34, z35, z45]
    eq1 = (1/neq/2)*(arr[0]*z[0]*z[0] + arr[1]*np.power(z[3]-z[1], 2))
    eq2 = (1/neq/2)*(arr[2]*z[1]*z[1] + arr[3]*np.power(z[3]-z[0], 2))
    eq3 = (1/neq/2)*(arr[4]*z[2]*z[2] + arr[5]*np.power(z[4]-z[0], 2))
    eq4 = (1/neq/2)*(arr[6]*z[3]*z[3] + arr[7]*np.power(z[1]-z[0], 2))
    eq5 = (1/neq/2)*(arr[8]*z[4]*z[4] + arr[9]*np.power(z[2]-z[0], 2))
    eq6 = (1/neq/2)*(arr[10]*z[5]*z[5] + arr[11]*np.power(z[2]-z[1], 2))

    return eq1 + eq2 + eq3 + eq4 + eq5 + eq6 - 1

def catch(b):
    cons = ()
    if b[1][4] == 1:
        cons += {'type':'eq','fun':lambda z: z[4]-z[0]}
        cons += {'type':'eq','fun':lambda z: z[5]-z[1]}
    return cons
# error to set to zero
def obj_func(z, arr, num_ones):
    z23, z24, z25, z34, z35, z45 = z
    t1 = 1 - arr[0]*z23*z23/2 - arr[1]*np.power(z34-z24,2)/2
    t2 = 1 - arr[2]*z24*z24/2 - arr[3]*np.power(z34-z23,2)/2
    t3 = 1 - arr[4]*z25*z25/2 - arr[5]*np.power(z35-z23,2)/2
    t4 = 1 - arr[6]*z34*z34/2 - arr[7]*np.power(z24-z23,2)/2
    t5 = 1 - arr[8]*z35*z35/2 - arr[9]*np.power(z25-z23,2)/2
    t6 = 1 - arr[10]*z45*z45/2 - arr[11]*np.power(z25-z24,2)/2

    return t1*t1 + t2*t2 + t3*t3 + t4*t4 + t5*t5 + t6*t6 - num_ones

def multistart_minimize(arr, n_eq, n_ones, n_trials=200):

    best = None
    cons = ({'type': 'eq', 'fun': lambda z: ellipse(z, arr, n_eq)}, {'type': 'eq', 'fun': lambda z: z[0] - z[2] - z[3] + z[5]},
            {'type': 'eq', 'fun': lambda z: z[1] - z[2] - z[3] + z[4]}) + catch(b)
    bnds = ((-2, 2), (-2, 2), (-2, 2), (-2, 2), (-2, 2), (-2, 2))
    for i in range(n_trials):
        init = np.random.uniform(-2., 2., (6))  # fix this
        now = optimize.minimize(obj_func, init, args=(arr, n_ones), bounds=bnds, constraints=cons)
        if now.success and (best is None or best.fun > now.fun):
            best = now
    return best

# i want to run multistart_minimize()
'''b = find_b(rmatrix)
arr = find_arr(b)
num_ones = len(np.where(np.abs(np.round(b, 4)) == 1)[0])/2
neq = n - num_ones # number of non-trivial ellipse equations
t1=time.time()
cons = ({'type': 'eq', 'fun': lambda z: ellipse(z, arr, neq)}, {'type': 'eq', 'fun': lambda z: z[0] - z[2] - z[3] + z[5]},
        {'type': 'eq', 'fun': lambda z: z[1] - z[2] - z[3] + z[4]}) + catch(b)
bnds = ((-2, 2), (-2, 2), (-2, 2), (-2, 2), (-2, 2), (-2, 2))
init = np.random.uniform(-2., 2., (6))  
res = optimize.minimize(obj_func, init, args=(arr, num_ones), bounds=bnds, constraints=cons)
print("Result", res)'''


'''for p1 in range(len(r)):
    l12 = r[p1]
    for p2 in range(len(r)):
        l13 = r[p2]
        for p3 in range(len(r)):
            l14 = r[p3]
            for p4 in range(len(r)):
                l15 = r[p4]
                for p5 in range(len(r)):
                    l23 = r[p5]
                    for p6 in range(len(r)):
                        l24 = r[p6]
                        for p7 in range(len(r)):
                            l25 = r[p7]
                            for p8 in range(len(r)):
                                l34 = r[p8]
                                for p9 in range(len(r)):
                                    l35 = r[p9]
                                    for p10 in range(len(r)):
                                        l45 = r[p10]
                                        rmatrix = np.array([[0,l12,l13,l14,l15],
                                                            [l12,0,l23,l24,l25],
                                                            [l13,l23,0,l34,l35],
                                                            [l14,l24,l34,0,l45],
                                                            [l15,l25,l35,l45,0]])
                                        b = find_b(rmatrix)
                                        arr = find_arr(b)
                                        num_ones = len(np.where(np.abs(np.round(b, 4)) == 1)[0]) / 2
                                        neq = n - num_ones  # number of non-trivial ellipse equations
                                        cons = ({'type': 'eq', 'fun': lambda z: ellipse(z, arr, neq)},
                                                {'type': 'eq', 'fun': lambda z: z[0] - z[2] - z[3] + z[5]},
                                                {'type': 'eq', 'fun': lambda z: z[1] - z[2] - z[3] + z[4]}) + catch(b)
                                        bnds = ((-2, 2), (-2, 2), (-2, 2), (-2, 2), (-2, 2), (-2, 2))
                                        init = np.random.uniform(-2., 2., (6))
                                        res = optimize.minimize(obj_func, init, args=(arr, num_ones), bounds=bnds,
                                                               constraints=cons)
                                        if res.fun < 1.e-11:
                                            print(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10)
                                            print("results ", res)

'''
b = find_b(rmatrix)
arr = find_arr(b)
num_ones = len(np.where(np.abs(np.round(b, 4)) == 1)[0]) / 2
neq = n - num_ones  # number of non-trivial ellipse equations
res = multistart_minimize(arr, neq, num_ones)
print("Result", res)