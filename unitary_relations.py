import numpy as np
import matplotlib.pyplot as plt

# 5th roots of unity values
a = (-1+np.sqrt(5))/4
b = np.sqrt((5+np.sqrt(5))/8)
c = -(1+np.sqrt(5))/4
d = np.sqrt((5-np.sqrt(5))/8)

# matrices for the SU(2)xSU(2) ~ SO(4) mapping
m = np.array([[1,1.j,0,0],[0,0,1.j,1],[0,0,1.j,-1],[1,-1.j,0,0]])/np.sqrt(2)
mT = np.asmatrix(m).getH()

# pauli
id   = np.array([[1, 0], [0, 1]])
sigx = np.array([[0, 1], [1, 0]])
sigy = np.array([[0, -1.j], [1.j, 0]])
sigz = np.array([[1, 0], [0, -1]])
pauli = [sigx, sigy, sigz]
# checks to see if two matrices are the same
def isequal(mat1,mat2):
    mat1 = np.round(mat1, 8)
    mat2 = np.round(mat2, 8)
    return np.allclose(mat1,mat2)

# given a set, this method will check if the matrix is in the set
def isin(mat,set):
    for i in set:
        if isequal(mat,i):
            return True
    return False

def closed_under_mult(set):
    for i in range(len(set)):
        for j in range(i+1,len(set)):
            if not isin(np.matmul(set[i],set[j]), set):
                print("("+str(i)+", "+str(j)+")")
                return False
    return True

def dist(vec):
    return np.sqrt(np.matmul(vec,vec))
'''
# 6 lines in R4 --------------------------------------------------------------------------------------------------------
# unitaries:
p = 0.5*(1+np.sqrt(5))
norm = np.sqrt(1+p*p)

o1 = (1/norm)*np.array([[0, 1-p*1.j], [1+1.j*p, 0]])
o2 = (1/norm)*np.array([[0, 1+1.j*p], [1-1.j*p, 0]])
o3 = (1/norm)*np.array([[-1.j, p], [p, -1.j]])
o4 = (1/norm)*np.array([[-1.j, -p], [-p, -1.j]])
o5 = (1/norm)*np.array([[1.j*p, -1.j], [1.j, 1.j*p]])
o6 = (1/norm)*np.array([[-1.j*p, -1.j], [1.j, -1.j*p]])
oset = [o1, o2, o3, o4, o5, o6]

obloch1 = np.array([1,p,0])/norm
obloch2 = np.array([1,-p,0])/norm
obloch3 = np.array([1,0,0])
obloch4 = np.array([-1,0,0])
obloch5 = np.array([0,1,0])
obloch6 = np.array([0,1,0])
obloch_set = [obloch1, obloch2, obloch3, obloch4, obloch5, obloch6]
'''
# define set of 5-cyclic matrices which generate the four equiangular lines in R4 --------------------------------------
generator = np.array([[a, -b, 0, 0], [b, a, 0, 0], [0, 0, c, -d], [0, 0, d, c]])
cyclic5    = [np.linalg.matrix_power(generator, i) for i in range(5)]
su2cyclic5 = [np.matmul(m, np.matmul(cyclic5[i], mT)) for i in range(len(cyclic5))] # this is indeed cyclic

# 5 unitaries in SU(2) which rotate Bell state phi^(+) to equiangular lines --------------------------------------------
u1 = np.array([[1.j, -1.j], [1.j, 1.j]])/np.sqrt(2)
u2 = np.array([[d+a*1.j, b-1.j*c], [b+1.j*c, -1*d+1.j*a]])/np.sqrt(2)
u3 = np.array([[-1*d+a*1.j, -b-1.j*c], [-b+1.j*c, d+a*1.j]])/np.sqrt(2)
u4 = np. array([[b+c*1.j, -d-a*1.j], [-d+a*1.j, -b+c*1.j]])/np.sqrt(2)
u5 = np.array([[-b+c*1.j, d-a*1.j], [d+a*1.j, b+c*1.j]])/np.sqrt(2)
u = [u1, u2, u3, u4, u5]

bloch1 = np.array([0, 1, 0])
bloch2 = np.array([b,c,d])/np.sqrt(b*b+c*c+d*d)
bloch3 = np.array([-b,c,-d])/np.sqrt(b*b+c*c+d*d)
bloch4 = np.array([-d,a,b])/np.sqrt(b*b+a*a+d*d)
bloch5 = np.array([d,a,-b])/np.sqrt(b*b+a*a+d*d)

bloch_set = [bloch1, bloch2, bloch3, bloch4, bloch5]

u1_inv = np.array([[-1.j, -1.j], [1.j, -1.j]])/np.sqrt(2)
'''
u2_inv = np.array([[d-a*1.j,b-c*1.j],[b+c*1.j,-a*1.j-d]])/np.sqrt(2)
u3_inv = np.array([[-d-a*1.j,-b-c*1.j],[-b+c*1.j,d-a*1.j]])/np.sqrt(2)
u4_inv = np.array([[b-c*1.j,-d-a*1.j],[-d+a*1.j,-b-c*1.j]])/np.sqrt(2)
'''
# shift the u_i so that the first one is the identity -----------------------------------------------------------------
v = [np.matmul(u[i], u1_inv) for i in range(len(u))]
bigv = [np.kron(v[i], np.array([[1, 0], [0, 1]])) for i in range(len(v))]
realv = [np.matmul(mT, np.matmul(bigv[i], m)) for i in range(len(bigv))]
evecs = np.transpose(np.linalg.eig(v[2]).eigenvectors)

evec1 = np.linalg.eig(v[0]).eigenvectors
evec2 = np.linalg.eig(v[1]).eigenvectors
evec3 = np.linalg.eig(v[2]).eigenvectors
evec4 = np.linalg.eig(v[3]).eigenvectors
evec5 = np.linalg.eig(v[4]).eigenvectors

''' # check that the states do not form a group
n = 1000000

gen = realv[2]
upow = gen
for i in range(1, n+1):
    if (np.round(upow, 8) == np.identity(4)).all():
        print("identity found at "+str(i)+"th power")
        #print(np.round(upow,4))
    #gen = np.matmul(gen, realv[1])
    #upow = np.matmul(upow, gen)
    upow = np.linalg.matrix_power(gen,i)'''

# bloch vectors --------------------------------------------------------------------------------------------------------

shifted_bloch1 = np.array([0, 0, 0])
shifted_bloch2 = np.array([-b-d, a-c, b-d])/np.sqrt((b+d)*(b+d)+(a-c)*(a-c)+(b-d)*(b-d))
shifted_bloch3 = np.array([b+d, a-c, -b+d])/np.sqrt((b+d)*(b+d)+(a-c)*(a-c)+(-b+d)*(-b+d))
shifted_bloch4 = np.array([-b+d, -a+c, -b-d])/np.sqrt((-b+d)*(-b+d)+(-a+c)*(-a+c)+(-b-d)*(-b-d))
shifted_bloch5 = np.array([b-d, -a+c, b+d])/np.sqrt((b-d)*(b-d)+(-a+c)*(-a+c)+(b+d)*(b+d))


shifted_bloch_set = [shifted_bloch1, shifted_bloch2, shifted_bloch3, shifted_bloch4, shifted_bloch5]

transpose_bloch1 = np.array([0,0,0])
transpose_bloch2 = np.array([b+d,- (a-c), -(b-d)])*(2/np.sqrt(15))
transpose_bloch3 = np.array([-b-d, -a+c, b-d])*(2/np.sqrt(15))
transpose_bloch4 = np.array([b-d, a-c, b+d])*(2/np.sqrt(15))
transpose_bloch5 = np.array([d-b, a-c, -b-d])*(2/np.sqrt(15))

transpose_bloch_set = [transpose_bloch1, transpose_bloch2, transpose_bloch3, transpose_bloch4, transpose_bloch5]

#-----------------------------------------------------------------------------------------------------------------------
'''b = qutip.Bloch()
b.make_sphere()
b.add_vectors(transpose_bloch_set)
b.add_vectors(shifted_bloch_set)
b.render()
plt.show()'''