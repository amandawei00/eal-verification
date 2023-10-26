import numpy as np
from qutip import *
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# 5th roots of unity values
a = (-1+np.sqrt(5))/4
b = np.sqrt((5+np.sqrt(5))/8)
c = -(1+np.sqrt(5))/4
d = np.sqrt((5-np.sqrt(5))/8)

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

o1_inv = np.linalg.inv(o1)
o2_inv = np.linalg.inv(o2)
o3_inv = np.linalg.inv(o3)
o4_inv = np.linalg.inv(o4)
o5_inv = np.linalg.inv(o5)
o6_inv = np.linalg.inv(o6)
o_invset = [o1_inv, o2_inv, o3_inv, o4_inv, o5_inv, o6_inv]

go12 = np.matmul(o2, o1_inv)
go13 = np.matmul(o3, o1_inv)
go14 = np.matmul(o4, o1_inv)
go15 = np.matmul(o5, o1_inv)
go16 = np.matmul(o6, o1_inv)

go21 = np.matmul(o1, o2_inv)
go23 = np.matmul(o3, o2_inv)
go24 = np.matmul(o4, o2_inv)
go25 = np.matmul(o5, o2_inv)
go26 = np.matmul(o6, o2_inv)

go31 = np.matmul(o1, o3_inv)
go32 = np.matmul(o2, o3_inv)
go34 = np.matmul(o4, o3_inv)
go35 = np.matmul(o5, o3_inv)
go36 = np.matmul(o6, o3_inv)

go41 = np.matmul(o1, o4_inv)
go42 = np.matmul(o2, o4_inv)
go43 = np.matmul(o3, o4_inv)
go45 = np.matmul(o5, o4_inv)
go46 = np.matmul(o6, o4_inv)

go51 = np.matmul(o1, o5_inv)
go52 = np.matmul(o2, o5_inv)
go53 = np.matmul(o3, o5_inv)
go54 = np.matmul(o4, o5_inv)
go56 = np.matmul(o6, o5_inv)

go61 = np.matmul(o1, o6_inv)
go62 = np.matmul(o2, o6_inv)
go63 = np.matmul(o3, o6_inv)
go64 = np.matmul(o4, o6_inv)
go65 = np.matmul(o5, o6_inv)

goset = [go12, go13, go14, go15, go16,
         go21, go23, go24, go25, go26,
         go31, go32, go34, go35, go36,
         go41, go42, go43, go45, go46,
         go51, go52, go53, go54, go56,
         go61, go62, go63, go64, go65]
'''
print("g12 = ")
print(go12)
print("g13 = ")
print(go13)
print("g14 = ")
print(go14)
print("g15 = ")
print(go15)
print("g16 = ")
print(go16)'''
# 5 lines in R4---------------------------------------------------------------------------------------------------------
u1 = np.array([[1.j,-1.j],[1.j,1.j]])/np.sqrt(2)
u2 = np.array([[d+a*1.j, b-1.j*c],[b+1.j*c,-1*d+1.j*a]])/np.sqrt(2)
u3 = np.array([[-1*d+a*1.j,-b-1.j*c],[-b+1.j*c,d+a*1.j]])/np.sqrt(2)
u4 = np. array([[b+c*1.j,-d-a*1.j],[-d+a*1.j,-b+c*1.j]])/np.sqrt(2)
u5 = np.array([[-b+c*1.j,d-a*1.j],[d+a*1.j, b+c*1.j]])/np.sqrt(2)

bloch1 = np.array([0, 1, 0])
bloch2 = np.array([b,c,d])/np.sqrt(b*b+c*c+d*d)
bloch3 = np.array([-b,c,-d])/np.sqrt(b*b+c*c+d*d)
bloch4 = np.array([-d,a,b])/np.sqrt(b*b+a*a+d*d)
bloch5 = np.array([d,a,-b])/np.sqrt(b*b+a*a+d*d)

bloch_set = [bloch1, bloch2, bloch3, bloch4, bloch5]

u1_inv = np.array([[-1.j, -1.j], [1.j, -1.j]])/np.sqrt(2)
u2_inv = np.array([[d-a*1.j,b-c*1.j],[b+c*1.j,-a*1.j-d]])/np.sqrt(2)
u3_inv = np.array([[-d-a*1.j,-b-c*1.j],[-b+c*1.j,d-a*1.j]])/np.sqrt(2)
u4_inv = np.array([[b-c*1.j,-d-a*1.j],[-d+a*1.j,-b-c*1.j]])/np.sqrt(2)

g12 = 0.5*np.array([[a+c+(b-d)*1.j,a-c-1.j*(b+d)],[-(a-c)-1.j*(b+d),a+c-1.j*(b-d)]])
g13 = np.matmul(u3,u1_inv)
g14 = 0.5*np.array([[a+c-1.j*(b+d),-(a-c)+1.j*(d-b)],[a-c+1.j*(d-b),a+c+1.j*(b+d)]])
g15 = 0.5*np.array([[a+c+1.j*(b+d),-(a-c)+1.j*(b-d)],[a-c+1.j*(b-d),a+c-1.j*(b+d)]])

g23 = np.matmul(u3,u2_inv)
g24 = 0.5*np.array([[2*a*c-2*a*b*1.j,-a*a+b*b+c*c+d*d+2*a*d*1.j],[a*a-b*b-c*c-d*d+2*a*d*1.j, 2*a*c+2*a*b*1.j]])
g25 = 0.5*np.array([[2*a*c+2*c*d*1.j, c*c-a*a-b*b-d*d+2*b*c*1.j],[a*a+b*b-c*c+d*d+2*b*c*1.j, 2*a*c-2*c*d*1.j]])

g34 = np.matmul(u4, u3_inv)
g35 = np.matmul(u5, u3_inv)

g45 = 0.5*np.array([[a*a-b*b+c*c-d*d+2*1.j*(a*d+c*b),2*1.j*(a*b-c*d)],[2*1.j*(a*b-c*d), a*a-b*b+c*c-d*d-2*1.j*(a*d+c*b)]])
'''
print("G12 = ", g12)
print("G13 = ", g13)
print("G14 = ", g14)
print("G15 = ", g15)

print("G23 = ", g23)
print("G24 = ", g24)
print("G25 = ", g25)

print("G34 = ", g34)
print("G35 = ", g35)

print("G45 = ", g45)
'''
# shifted bases---------------------------------------------------------------------------------------------------------
v1 = np.matmul(u1, u1_inv)
v2 = np.matmul(u2, u1_inv)
v3 = np.matmul(u3, u1_inv)
v4 = np.matmul(u4, u1_inv)
v5 = np.matmul(u5, u1_inv)

shifted_bloch1 = np.array([0,0,0])
shifted_bloch2 = np.array([-b-d, a-c, b-d])/np.sqrt((-b-d)*(-b-d)+(a-c)*(a-c)+(b-d)*(b-d))
shifted_bloch3 = np.array([b+d, a-c, a+c])/np.sqrt((b+d)*(b+d)+(a-c)*(a-c)+(a+c)*(a+c))
shifted_bloch4 = np.array([-b+d, -a+c, -b-d])/np.sqrt((-b+d)*(-b+d)+(-a+c)*(-a+c)+(-b-d)*(-b-d))
shifted_bloch5 = np.array([b-d, -a+c, b+d])/np.sqrt((b-d)*(b-d)+(-a+c)*(-a+c)+(b+d)*(b+d))

shifted_bloch_set = [shifted_bloch1, shifted_bloch2, shifted_bloch3, shifted_bloch4, shifted_bloch5]


b = qutip.Bloch()
b.make_sphere()
b.add_vectors(obloch_set)
b.render()
plt.show()