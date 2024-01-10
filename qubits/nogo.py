import numpy as np


mu = 1/9

# allowed magnitude for ni \cdot nj (equation 1)
a1 = np.sqrt(mu)/(1+np.sqrt(mu))  # \sqrt{mu}/(1+\sqrt{mu})
a2 = np.sqrt(mu)/(1-np.sqrt(mu))  # \sqrt{mu}/(1-\sqrt{mu})

# including sign
nij = [a1, -1*a1, a2, -1*a2]

mu19solution = np.array([[-np.sqrt(15)/4,0,-1/4],
                         [-1*np.sqrt(3/5)/2,-1*np.sqrt(3/5), -1/2],
                         [-1*np.sqrt(3/5)/2, -1*np.sqrt(3/5), 1/2],
                         [-3*np.sqrt(3/5)/4, np.sqrt(3/5), 1/4]])
def get_r(iv):
   # generates a matrix of ni \cdot nj values, like the one shown in equation 9
   v = [nij[iv[0]], nij[iv[1]], nij[iv[2]], nij[iv[3]], nij[iv[4]], nij[iv[5]], nij[iv[6]], nij[iv[7]], nij[iv[8]], nij[iv[9]]]
   return np.array([[0, v[0], v[1], v[2], v[4]],
                    [v[0], 0, v[4], v[5], v[6]],
                    [v[1], v[4], 0, v[7], v[8]],
                    [v[2], v[5], v[7], 0, v[9]],
                    [v[4], v[6], v[8], v[9], 0]])

def gen_cp(r):

    # the cosines in equation 4 define the corresponding sines by sqrt(1-cos^2(theta)).
    s1 = np.sqrt(1 - r[0][1]*r[0][1])
    s2 = np.sqrt(1 - r[0][2]*r[0][2])
    s3 = np.sqrt(1 - r[0][3]*r[0][3])
    s4 = np.sqrt(1 - r[0][4]*r[0][4])

    # these are the cos(phi_i) in equation 6 for i = 3, 4, 5
    cp3 = (r[1][2] - r[0][1]*r[0][2])/s1/s2
    cp4 = (r[1][3] - r[0][1]*r[0][3])/s1/s3
    cp5 = (r[1][4] - r[0][1]*r[0][4])/s1/s4

    # since the s1, s2, s3, s4 may be positive or negative, there are 2^3 possible sign combinations for the cos(phi_i)
    # for a given table from the get_r() method, and these are returned from this method
    return [[cp3, cp4, cp5], [cp3, cp4, -1*cp5], [cp3, -1*cp4, cp5], [cp3, -1*cp4, -1*cp5],
            [-1*cp3, cp4, cp5], [-1*cp3, cp4, -1*cp5], [-1*cp3, -1*cp4, cp5], [-1*cp3, -1*cp4, -1*cp5]]

def gen_sp(cp):
    # for one of the combinations of cos(phi_i) (i=3,4,5) from above, the associated sin(phi_i) are defined here
    sp3 = np.sqrt(1-np.round(cp[0]*cp[0],8))
    sp4 = np.sqrt(1-np.round(cp[1]*cp[1], 8))
    sp5 = np.sqrt(1-np.round(cp[2]*cp[2],8))

    # again, they are only defined up to sign so there are 2^3 possible sign combinations
    return [[sp3, sp4, sp5], [sp3, sp4, -1*sp5], [sp3, -1*sp4, sp5], [sp3, -1*sp4, -1*sp5],
            [-1*sp3, sp4, sp5], [-1*sp3, sp4, -1*sp5], [-1*sp3, -1*sp4, sp5], [-1*sp3, -1*sp4, -1*sp5]]

def mat_equal(a, b):
    return np.allclose(np.round(a, 8), np.round(b, 8))

# for loops generate all possible 4^10 possible combinations of n_{i} \cdot n_{j}
# len(nij) = 4
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
                                        # produces the table of ni \cdot nj values
                                        r = get_r([r12, r13, r14, r15, r23, r24, r25, r34, r35, r45])
                                        cp = gen_cp(r)

                                        # these are the same sines defined in lines 25-28... i probably only needed to define these once but ok
                                        st2 = np.sqrt(1 - r[0][1]*r[0][1])
                                        st3 = np.sqrt(1 - r[0][2]*r[0][2])
                                        st4 = np.sqrt(1 - r[0][3]*r[0][3])
                                        st5 = np.sqrt(1 - r[0][4]*r[0][4])
                                        # tests for a consistent set of phi_i's for each possible sign combination of the cos(phi_i) and sin(phi_i)'s
                                        for i in cp:
                                            sp = gen_sp(i)
                                            for j in sp:
                                                # quantities are rounded to 6 decimal places
                                                # vector i is now one possible [cos(phi3), cos(phi4), cos(phi5)] from the gen_cp() method
                                                # vector j is now one possible [sin(phi3), sin(phi4), sin(phi5)] from the gen_sp() method
                                                # these are testing equality in equation 8. the first tests equality for index i=3, j=4 (unfortunate choice of variabls, im sorry)
                                                if np.round(st3*st4*(i[0]*i[1]+j[0]*j[1])+r[0][2]*r[0][3], 8) == np.round(r[2][3], 8):
                                                    # tests equality for index i=3, j=5
                                                    if np.round(st3*st5*(i[0]*i[2]+j[0]*j[2])+r[0][2]*r[0][4], 8) == np.round(r[2][4], 8):
                                                        # tests equality for index i=4, j=5
                                                        if np.round(st4*st5*(i[1]*i[2]+j[1]*j[2])+r[0][3]*r[0][4], 8) == np.round(r[3][4], 8):
                                                            # if all are consistent then a solution has been found
                                                            n2 = [st2, 0, r[0][1]]
                                                            n3 = [st3*i[0], st3*j[0], r[0][2]]
                                                            n4 = [st4*i[1], st4*j[1], r[0][3]]
                                                            n5 = [st5*i[2], st5*j[2], r[0][4]]
                                                            sol = np.array([n2, n3, n4, n5])
                                                            if mat_equal(mu19solution, sol):
                                                                print(np.round(sol, 5))
# otherwise, no solution found
print("no solution found")
'''
# for loops generate all possible 4^10 possible combinations of n_{i} \cdot n_{j}
# len(nij) = 4

# produces the table of ni \cdot nj values
r = np.array([[1, 1/2, 1/4, -1/2, 1/4],
              [1/2, 1, 1/2, 1/2, 1/2],
              [1/4, 1/2, 1, 1/4, -1/2],
              [-1/2, 1/2, 1/4, 1, 1/4],
              [1/4, 1/2, -1/2, 1/4, 1]])
cp = gen_cp(r)

# these are the same sines defined in lines 25-28... i probably only needed to define these once but ok
st2 = np.sqrt(1 - r[0][1]*r[0][1])
st3 = np.sqrt(1 - r[0][2]*r[0][2])
st4 = np.sqrt(1 - r[0][3]*r[0][3])
st5 = np.sqrt(1 - r[0][4]*r[0][4])
# tests for a consistent set of phi_i's for each possible sign combination of the cos(phi_i) and sin(phi_i)'s
for i in cp:
    sp = gen_sp(i)
    for j in sp:
        # quantities are rounded to 6 decimal places
        # vector i is now one possible [cos(phi3), cos(phi4), cos(phi5)] from the gen_cp() method
        # vector j is now one possible [sin(phi3), sin(phi4), sin(phi5)] from the gen_sp() method
        # these are testing equality in equation 8. the first tests equality for index i=3, j=4 (unfortunate choice of variabls, im sorry)
        if np.round(st3*st4*(i[0]*i[1]+j[0]*j[1])+r[0][2]*r[0][3], 6) == np.round(r[2][3], 6):
            # tests equality for index i=3, j=5
            if np.round(st3*st5*(i[0]*i[2]+j[0]*j[2])+r[0][2]*r[0][4], 6) == np.round(r[2][4], 6):
                # tests equality for index i=4, j=5
                if np.round(st4*st5*(i[1]*i[2]+j[1]*j[2])+r[0][3]*r[0][4], 6) == np.round(r[3][4], 6):
                    # if all are consistent then a solution has been found
                    print("n2: cos(theta) = " + str(r[0][1]) + ", sin(theta) = " + str(st2))
                    print("n3: cos(theta) = " + str(r[0][2]) + ", sin(theta) = " + str(st3) + ", cos(phi) = " + str(i[0]) + ", sin(phi) = " + str(j[0]))
                    print("n4: cos(theta) = " + str(r[0][3]) + ", sin(theta) = " + str(st4) + ", cos(phi) = " + str(i[1]) + ", sin(phi) = " + str(j[1]))
                    print("n5: cos(theta) = " + str(r[0][4]) + ", sin(theta) = " + str(st5) + ", cos(phi) = " + str(i[2]) + ", sin(phi) = " + str(j[2]))
                    print("------------------------------------------------------------------------------------------------------------------------------")
'''