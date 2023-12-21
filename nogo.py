import numpy as np


mu = 1/10
a1 = np.sqrt(mu)/(1+np.sqrt(mu))
a2 = np.sqrt(mu)/(1-np.sqrt(mu))
nij = [a1, -1*a1, a2, -1*a2]

def get_r(iv):
   # each position in vector v takes on one of four possible values in nij
   # v: r12, r13, r14, r15, r23, r24, r25, r34, r35, r45
   v = [nij[iv[0]], nij[iv[1]], nij[iv[2]], nij[iv[3]], nij[iv[4]], nij[iv[5]], nij[iv[6]], nij[iv[7]], nij[iv[8]], nij[iv[9]]]
   return np.array([[0, v[0], v[1], v[2], v[4]],
                    [v[0], 0, v[4], v[5], v[6]],
                    [v[1], v[4], 0, v[7], v[8]],
                    [v[2], v[5], v[7], 0, v[9]],
                    [v[4], v[6], v[8], v[9], 0]])

def gen_cp(r):
    s1 = np.sqrt(1 - r[0][1]*r[0][1])
    s2 = np.sqrt(1 - r[0][2]*r[0][2])
    s3 = np.sqrt(1 - r[0][3]*r[0][3])
    s4 = np.sqrt(1 - r[0][4]*r[0][4])

    cp3 = (r[1][2] - r[0][1]*r[0][2])/s1/s2
    cp4 = (r[1][3] - r[0][1]*r[0][3])/s1/s3
    cp5 = (r[1][4] - r[0][1]*r[0][4])/s1/s4

    return [[cp3, cp4, cp5], [cp3, cp4, -1*cp5], [cp3, -1*cp4, cp5], [cp3, -1*cp4, -1*cp5],
            [-1*cp3, cp4, cp5], [-1*cp3, cp4, -1*cp5], [-1*cp3, -1*cp4, cp5], [-1*cp3, -1*cp4, -1*cp5]]

def gen_sp(cp):
    sp3 = np.sqrt(1-cp[0]*cp[0])
    sp4 = np.sqrt(1-cp[1]*cp[1])
    sp5 = np.sqrt(1-cp[2]*cp[2])

    return [[sp3, sp4, sp5], [sp3, sp4, -1*sp5], [sp3, -1*sp4, sp5], [sp3, -1*sp4, -1*sp5],
            [-1*sp3, sp4, sp5], [-1*sp3, sp4, -1*sp5], [-1*sp3, -1*sp4, sp5], [-1*sp3, -1*sp4, -1*sp5]]

for r12 in range(len(nij)):
    for r13 in range(len(nij)):
        for r14 in range(len(nij)):
            for r15 in range(len(nij)):
                for r23 in range(len(nij)):
                    for r24 in range(len(nij)):
                        for r25 in range(len(nij)):
                            for r34 in range(len(nij)):
                                for r35 in range(len(nij)):
                                    for r45 in range(len(nij)):
                                        r = get_r([r12, r13, r14, r15, r23, r24, r25, r34, r35, r45])
                                        cp = gen_cp(r)

                                        st2 = np.sqrt(1 - r[0][1]*r[0][1])
                                        st3 = np.sqrt(1 - r[0][2]*r[0][2])
                                        st4 = np.sqrt(1 - r[0][3]*r[0][3])
                                        st5 = np.sqrt(1 - r[0][4]*r[0][4])
                                        for i in cp:
                                            sp = gen_sp(i)
                                            for j in sp:
                                                if np.round(st3*st4*(i[0]*i[1]+j[0]*j[1])+r[0][2]*r[0][3], 6) == np.round(r[2][3], 6):
                                                    if np.round(st3*st5*(i[0]*i[2]+j[0]*j[2])+r[0][2]*r[0][4], 6) == np.round(r[2][4], 6):
                                                        if np.round(st4*st5*(i[1]*i[2]+j[1]*j[2])+r[0][3]*r[0][4], 6) == np.round(r[3][4], 6):
                                                            print("consistent solution found")
print("no solution found")