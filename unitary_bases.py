import numpy as np

# linear combination of matrices in mat with coefficients in vec
# vec and mat should be same length vectors
def lin_comb(vec,mat):
    sum = np.zeros_like(mat[0])
    for i in range(len(vec)):
        sum += vec[i]*mat[i]
    return sum

def is_unitary(mat):
    print(np.matmul(np.asmatrix(mat).getH(), mat))

# pauli matrices (2x2)--------------------------------------------------------------------
id_p = np.array([[1.j,0],[0,1.j]])
sigx = np.array([[0,1],[1,0]])
sigy = np.array([[0,-1.j],[1.j,0]])
sigz = np.array([[1,0],[0,-1]])

#pauli = np.array([id_p, sigx, sigy, sigz])
#test_pauli = lin_comb(np.array([1/2, 1/2, 1/2, 1/2]), pauli)
#print(is_unitary(test_pauli))

# gell mann (3x3)--------------------------------------------------------------------
l1 = np.array([[0,1,0],[1,0,0],[0,0,0]])
l2 = np.array([[0,-1.j,0],[1.j,0,0],[0,0,0]])
l3 = np.array([[1,0,0],[0,-1,0],[0,0,0]])
l4 = np.array([[0,0,1],[0,0,0],[1,0,0]])
l5 = np.array([[0,0,-1.j],[0,0,0],[1.j,0,0]])
l6 = np.array([[0,0,0],[0,0,1],[0,1,0]])
l7 = np.array([[0,0,0],[0,0,-1.j],[0,1.j,0]])
l8 = np.array([[1,0,0],[0,1,0],[0,0,-2]])/np.sqrt(3)

l9=l3-l8
#h1 = 0.5*(l1+l2+l4+l5+l6+l7)
#h2 = 0.5*(l1-l2+l4-l5+l6-l7)
print(np.linalg.eig(-l1+l9))

# gellmann = np.array([l0, l1, l2, l3, l4, l5, l6, l7, l8])
# test_gellmann = lin_comb(np.ones(9)/3,gellmann)


# schwinger (3x3)--------------------------------------------------------------------
