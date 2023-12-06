import numpy as np
import scipy.optimize as optimize

mu = 1/10
n = 6
r = [np.sqrt(mu)/(1+np.sqrt(mu)),-np.sqrt(mu)/(1+np.sqrt(mu)), np.sqrt(mu)/(1-np.sqrt(mu)), -np.sqrt(mu)/(1-np.sqrt(mu))]

#rmat = np.array([[0, r[0], r[1], r[1], r[3]],
#                 [r[0], 0, -r[1], -r[1], r[1]],
#                 [r[1], r[3], 0, r[2], r[2]],
#                 [r[1], -r[1], r[2], 0, r[1]],
#                 [r[3], r[1], r[2], r[1], 0]])
# Armin's numerical solution
rmat = np.array([[0,1/2,1/4,-1/2,1/4], [1/2,0,1/2,1/2,1/2], [1/4,1/2,0,1/4,-1/2], [-1/2,1/2,1/4,0,1/4], [1/4,1/2,-1/2,1/4,0]])
b = np.empty((5,5))
for i in range(len(b)):
    for j in range(len(b)):
        if (i != j) & (i > 0) & (j > 0):
            b[i][j] = (rmat[i][j]-rmat[i][0]*rmat[j][0])/(np.sqrt(1-rmat[i][0]*rmat[i][0])*np.sqrt(1-rmat[j][0]*rmat[j][0]))
        else:
            b[i][j] = 0

arr = [1/(1+b[1][2]), 1/(1-b[1][2]), 1/(1+b[1][3]), 1/(1-b[1][3]), 1/(1+b[1][4]), 1/(1-b[1][4]), 1/(1+b[2][3]), 1/(1-b[2][3]), 1/(1+b[2][4]), 1/(1-b[2][4]), 1/(1+b[3][4]), 1/(1-b[3][4])]
num_ones = len(np.where(np.abs(np.round(b, 4)) == 1)[0])/2
neq = n - num_ones # number of non-trivial ellipse equations

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

# (non transformed) ellipse condition + 2 leftover conditions + conditions imposed when beta = \pm 1
def ellipse(z):
    # z = [z23, z24, z25, z34, z35, z45]
    eq1 = (1/neq/2)*(arr[0]*z[0]*z[0] + arr[1]*np.power(z[3]-z[1], 2))
    eq2 = (1/neq/2)*(arr[2]*z[1]*z[1] + arr[3]*np.power(z[3]-z[0], 2))
    eq3 = (1/neq/2)*(arr[4]*z[2]*z[2] + arr[5]*np.power(z[4]-z[0], 2))
    eq4 = (1/neq/2)*(arr[6]*z[3]*z[3] + arr[7]*np.power(z[1]-z[0], 2))
    eq5 = (1/neq/2)*(arr[8]*z[4]*z[4] + arr[9]*np.power(z[2]-z[0], 2))
    eq6 = (1/neq/2)*(arr[10]*z[5]*z[5] + arr[11]*np.power(z[2]-z[1], 2))

    return eq1 + eq2 + eq3 + eq4 + eq5 + eq6 - 1

def catch():
    cons = ()
    if b[1][4] == 1:
        cons += {'type':'eq','fun':lambda z: z[4]-z[0]}
        cons += {'type':'eq','fun':lambda z: z[5]-z[1]}
    return cons
# error to set to zero
def obj_func(z):
    z23, z24, z25, z34, z35, z45 = z
    t1 = 1 - arr[0]*z23*z23/2 - arr[1]*np.power(z34-z24,2)/2
    t2 = 1 - arr[2]*z24*z24/2 - arr[3]*np.power(z34-z23,2)/2
    t3 = 1 - arr[4]*z25*z25/2 - arr[5]*np.power(z35-z23,2)/2
    t4 = 1 - arr[6]*z34*z34/2 - arr[7]*np.power(z24-z23,2)/2
    t5 = 1 - arr[8]*z35*z35/2 - arr[9]*np.power(z25-z23,2)/2
    t6 = 1 - arr[10]*z45*z45/2 - arr[11]*np.power(z25-z24,2)/2

    #print(t1)
    #print(t2)
    #print(t3)
    ##print(t4)
    #print(t5)
    #print(t6)
    return t1*t1 + t2*t2 + t3*t3 + t4*t4 + t5*t5 + t6*t6 - num_ones

def multistart_minimize(n_trials=500):
    best = None
    cons = ({'type': 'eq', 'fun': lambda z: ellipse(z)}, {'type': 'eq', 'fun': lambda z: z[0] - z[2] - z[3] + z[5]},
            {'type': 'eq', 'fun': lambda z: z[1] - z[2] - z[3] + z[4]}) + catch()
    bnds = ((-2, 2), (-2, 2), (-2, 2), (-2, 2), (-2, 2), (-2, 2))
    for i in range(n_trials):
        print('trial ' + str(i))
        init = np.random.uniform(-2., 2., (6))  # fix this
        now = optimize.minimize(obj_func, init, bounds=bnds, constraints=cons)
        if now.success and (best is None or best.fun > now.fun):
            best = now
    return best

res = multistart_minimize()
print("Result", res)