import numpy as np

def findxy(z):
    # z = [z23,z24,z25,z34,z35,z45]
    zbar = np.array([z[3]-z[1], z[3]-z[0], z[4]-z[0], z[1]-z[0], z[2]-z[1]])
    x2 = (z[0]-zbar[0])/2
    x3 = (z[0]+zbar[0])/2
    x4 = (z[1]+zbar[1])/2
    x5 = (z[2]+zbar[2])/2
    print(x2)
    print(x3)
    print(x4)
    print(x5)
    #consistency check
    if (z[3] != x3 + x4) or (zbar[3] != -1*x3 + x4):
        print("inconsistent: z[3] = " + str(z[3]) + ", zbar[3] = " + str(zbar[3]) + ", x3 = " + str(x3) + ", x4 = " + str(x4))
        return False
    if (z[4] != x3 + x5) or (zbar[4] != -1*x3 + x5):
        print("inconsistent: z[4] = " + str(z[4]) + ", zbar[4] = " + str(zbar[4]) + ", x3 = " + str(x3) + ", x5 = " + str(x5))
        return False
    if (z[5] != x4 + x5) or (zbar[5] != -1*x4 + x5):
        print("inconsistent: z[5] = " + str(z[5]) + ", zbar[5] = " + str(zbar[5]) + ", x4 = " + str(x4) + ", x5 = " + str(x5))
        return False

    return [x2, x3, x4, x5]


def find_n(x, y_par, cos, sin):
    y = [y_par[i]*np.sqrt(1-x[i]*x[i]) for i in range(len(x))]
    n2 = np.array([sin[0]*x[0], sin[0]*y[0], cos[0]])
    n3 = np.array([sin[1]*x[1], sin[1]*y[1], cos[1]])
    n4 = np.array([sin[2]*x[2], sin[2]*y[2], cos[2]])
    n5 = np.array([sin[3]*x[3], sin[3]*y[3], cos[3]])
    return [n2, n3, n4, n5]

z = [0.144039498, 1.1847, 1.0438, -0.5868, -0.727759, 0.3129]
x = findxy(z)

y_par = np.array([1,1,-1,1])
cos = np.array([1/2,1/4,-1/2,1/4])
sin = np.array([np.sqrt(1-cos[i]*cos[i]) for i in range(len(cos))])
n = find_n(x, y_par, cos, sin)

print(np.matmul(n[0], n[1]))
print(np.matmul(n[0], n[2]))
print(np.matmul(n[0], n[3]))
print(np.matmul(n[1], n[2]))
print(np.matmul(n[1], n[3]))
print(np.matmul(n[2], n[3]))


