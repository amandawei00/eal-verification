import numpy as np

# these quantities are determined by mu

alp2 = 3/10
rho2 = 4/10

alp = np.sqrt(alp2)
rho = np.sqrt(rho2)

ref_mat = np.array([[rho2, alp2, alp2],[alp2, rho2, alp2],[alp2, alp2, rho2]])

def gen_m(r):
    # r: [(0)b, (1)c, (2)qab, (3)rab, (4)sab, (5)tab, (6)uab, (7)vab, (8)wab]
    # b and c \in {0,1} are sign choices for beta and gamma
    # the rest are phase choices. For now, in the conditions we take positive signs for both off-diagonal unitary conditions
    # the variables should be assigned so cos(p-r-v+x)=cos(t-u-w+x)=cos(p-q-s+t)=-5/8
    bab, cab, qab, rab, sab, tab, uab, vab, wab = r[0], r[1], r[2], r[3], r[4], r[5], r[6], r[7], r[8]
    xab = tab + rab + vab - qab - sab
    pab = tab - uab - wab + rab + vab

    return np.array([[rho*np.exp(1.j*pab), alp*np.exp(1.j*qab), np.power(-1, bab)*alp*np.exp(1.j*rab)],
                     [alp*np.exp(1.j*sab), rho*np.exp(1.j*tab), np.power(-1, cab)*alp*np.exp(1.j*uab)],
                     [np.power(-1, bab)*alp*np.exp(1.j*vab), np.power(-1,cab)*alp*np.exp(1.j*wab), rho*np.exp(1.j*xab)]])


def is_consistent(r1, r2):
    # r1: M1 -> 2
    # r2: M1 -> 3
    # this method checks to see if the product of these two transformation M2 -> 3 also has the correct magnitudes

    m = np.asarray(np.matmul(np.asmatrix(gen_m(r1)).getH(), gen_m(r2)))
    diff = np.zeros((3, 3))
    for i in range(len(m)):
        for j in range(len(m[0])):
            if i == j:
                diff[i][j] = np.round(np.power(np.abs(m[i][j]), 2),6)
            else:
                diff[i][j] = np.round(np.power(np.abs(m[i][j]), 2),6)
    return diff

# search for solutions first assuming parity is zero
b = 0
c = 0

# generate random r vectors
for i in range(1000):
    r1 = np.concatenate(([b,c], np.random.rand(7)))
    r2 = np.concatenate(([b,c], np.random.rand(7)))
    #r3 = np.concatenate(([b, c], np.random.rand(7)))
    #r4 = np.concatenate(([b, c], np.random.rand(7)))
    #r5 = np.concatenate(([b, c], np.random.rand(7)))
    #r6 = np.concatenate(([b, c], np.random.rand(7)))
    #r7 = np.concatenate(([b, c], np.random.rand(7)))
    #r8 = np.concatenate(([b, c], np.random.rand(7)))
    #r9 = np.concatenate(([b, c], np.random.rand(7))

    #print(is_consistent(r1,r2))
    if (is_consistent(r1, r2) == ref_mat).all():
        print("solution " + str(i) + "------------------------------------")
        print(r1)
        print(r2)

