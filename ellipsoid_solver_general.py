import numpy as np
import scipy.optimize as optimize

mu = 1/9
r = [np.sqrt(mu)/(1+np.sqrt(mu)),-np.sqrt(mu)/(1+np.sqrt(mu)), np.sqrt(mu)/(1-np.sqrt(mu)), -np.sqrt(mu)/(1-np.sqrt(mu))]
rmat = np.array([[0,1/2,1/4,-1/2,1/4],
                 [1/2,0,1/2,1/2,1/2],
                 [1/4,1/2,0,1/4,-1/2],
                 [-1/2,1/2,1/4,0,1/4],
                 [1/4,1/2,-1/2,1/4,0]])
# rmat = np.array([[0, 1/4, 1/2, 1/2, -1/2],
#                 [1/4, 0, -1/4, 1/4, 1/2],
#                 [1/2, -1/4, 0, 1/4, 1/4],
#                 [1/2, 1/4, 1/4, 0, -1/2],
#                 [-1/2, 1/2, 1/4, -1/2, 0]])
b = np.empty((5,5))
for i in range(len(b)):
    for j in range(len(b)):
        if (i != j) & (i > 0) & (j > 0):
            b[i][j] = (rmat[i][j]-rmat[i][0]*rmat[j][0])/(np.sqrt(1-rmat[i][0]*rmat[i][0])*np.sqrt(1-rmat[j][0]*rmat[j][0]))
        else:
            b[i][j] = 0

arr = [1/(1+b[1][2]), 1/(1-b[1][2]), 1/(1+b[1][3]), 1/(1-b[1][3]), 1/(1+b[1][4]), 1/(1-b[1][4]), 1/(1+b[2][3]), 1/(1-b[2][3]), 1/(1+b[2][4]), 1/(1-b[2][4]), 1/(1+b[3][4]), 1/(1-b[3][4])]
num_ones = np.size(np.where(np.abs(np.round(arr,4))==1))
print(num_ones)
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

a1 = (1/10)*(arr[0]+arr[1]+arr[2]+arr[3]+arr[4]+arr[5])
a2 = (1/10)*(arr[0]+arr[2]+arr[5]+arr[6]+arr[9]+arr[11])
a3 = (1/10)*(arr[1]+arr[2]+arr[5]+arr[7]+arr[8]+arr[11])
a4 = (1/10)*(arr[0]+arr[3]+arr[5]+arr[7]+arr[9]+arr[10])
a5 = (1/5)*(arr[0]+arr[2]+arr[5])
a6 = (-1/5)*(arr[1]+arr[2]+arr[5])
a7 = (-1/5)*(arr[0]+arr[3]+arr[5])
a8 = (-1/5)*(arr[2]+arr[5]+arr[11])
a9 = (-1/5)*(arr[0]+arr[5]+arr[9])
a10 = (1/5)*(arr[5]-arr[7])

qmat = np.array([[a1, a5/2, a6/2, a7/2, 0],[a5/2, a2, a8/2, a9/2, 0],[a6/2, a8/2, a3, a10/2, 0],[a7/2, a9/2, a10/2, a4, 0],[0,0,0,0, -1]])
eval, evec = np.linalg.eig(qmat)
qd = np.diag(eval)
print(qd)
# print(np.round(np.matmul(np.linalg.inv(evec), np.matmul(qmat, evec)),4))

# final ellipse condition
def constraint(z):
    # transform z = [z25, z34, z35, z45] to diagonalized coord.
    zd = np.matmul(np.linalg.inv(evec), np.append(z, [1]))
    zd25, zd34, zd35, zd45, one = zd
    return zd25*zd25/qd[0][0] + zd34*zd34/qd[1][1] + zd35*zd35/qd[2][2] + zd45*zd45/qd[3][3] - 1

# error to set to zero
def obj_func(z):
    z23, z24, z25, z34, z35, z45 = z
    #t1,t2,t3,t4 = 0,0,0,0

    if np.abs(np.round(b[1][2], 4)) == 1:
        t1 = 0
        z25 = z35
    else:
        t1 = np.power(1 - np.power(z25+z34-z45, 2)/(2*(1+b[1][2])) - np.power(z35-z25, 2)/(2*(1-b[1][2])), 2)

    if np.abs(np.round(b[1][3], 4)) == 1:
        t2 = 0
        z25 = z45
    else:
        t2 = np.power(1 - np.power(z25+z34-z35, 2)/(2*(1+b[1][3])) - np.power(z45-z25, 2)/(2*(1-b[1][3])), 2)

    if np.round(b[1][4], 4) == 1:
        t3 = np.power(1 - np.power(z25, 2)/4, 2)
    elif np.round(b[1][4], 4) == -1:
        t3 = np.power(1 - np.power(z35-z25-z34+z45, 2)/4, 2)
    else:
        t3 = np.power(1 - z25*z25/(2*(1+b[1][4])) - np.power(z35-z25-z34+z45, 2)/(2*(1-b[1][4])), 2)

    if np.round(b[2][3], 4) == 1:
        t4 = np.power(1 - np.power(z34, 2)/4, 2)
    elif np.round(b[2][3], 4) == -1:
        t4 = np.power(1 - np.power(z45-z35, 2)/4, 2)
    else:
        t4 = np.power(1 - np.power(z34, 2)/(2*(1+b[2][3])) - np.power(z45-z35, 2)/(2*(1-b[2][3])), 2)

    if np.round(b[2][4], 4) == 1:
        t5 = np.power(1 - np.power(z35, 2)/4, 2)
    elif np.round(b[2][4], 4) == -1:
        t5 = np.power(1 - np.power(z45-z34,2)/4, 2)
    else:
        t5 = np.power(1 - z35*z35/(2*(1+b[2][4])) - np.power(z45-z34, 2)/(2*(1-b[2][4])), 2)

    if np.round(b[3][4], 4) == 1:
        t6 = np.power(1 - np.power(z45,2)/4, 2)
    elif np.round(b[3][4], 1) == -1:
        t6 = np.power(1 - np.power(z35-z34,2)/4, 2)
    else:
        t6 = np.power(1 - z45*z45/(2*(1+b[3][4])) - np.power(z35-z34,2)/(2*(1-b[3][4])), 2)
    return t1 + t2 + t3 + t4 + t5 + t6

init = [0, 0, 0, np.sqrt(2)]
cons = ({'type':'eq', 'fun': lambda z: constraint(z)})
bnds = ((-2, 2),(-2, 2),(-2, 2),(-2, 2))
res = optimize.minimize(obj_func, init, bounds=bnds, constraints = cons)
print("Result", res)