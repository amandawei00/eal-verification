import numpy as np

#a and b must both be n-dim vectors
def dot(a,b):
    s = 0
    for i in range(len(a)):
        s += a[i]*b[i]
    return s

# input: set of any number of nx1 vectors
# returns yes and prints the scalar product if the vectors are equiangular, returns no if not
def equiangular(a):
    ref = round(np.abs(dot(a[1],a[2])),5)
    print("angle: " + str(ref))
    for i in range(len(a)-1):
        for j in range(i+1, len(a)):
            inverse = [-1 * b for b in a[j]]
            if round(np.abs(dot(a[i], a[j])), 5) != ref and round(np.abs(dot(a[i], inverse)), 5) != ref:
                print(dot(a[i], a[j]))
                print("not equiangular")
                return False
            print("vectors " + str(i) + ", " + str(j) + ": " + str(round(np.abs(dot(a[i],a[j])),5)))

    print("yes equiangular")

# test of equiangular lines
# a = ((1,0,0),(0,1,0),(0,0,1)) # orthoggonal basis in d=3
# a = ((1,0),(0.5, np.sqrt(3)/2),(-0.5, np.sqrt(3)/2)) # set of equiangular lines in d=2

# corners of the icosahedron, maximal equiangular lines in d=3
# p = 0.5*(1 + np.sqrt(5))

# a = ((0, 1, p), (0, 1, -1 * p), (-1, p, 0), (-1, -p, 0), (p, 0, 1), (-1 * p, 0, 1))
# equiangular(a)

a = 1/2
b = np.sqrt(3)/2

u1 = np.array([1,0,-1, 0])/np.sqrt(3)
u2 = np.array([a, b, -1, 0])/np.sqrt(3)
u3 = np.array([-1*a, b, -1, 0])/np.sqrt(3)
u4 = np.array([-1,0,-1*a,-1*b])/np.sqrt(3)
u5 = np.array([-1*a,-1*b,-1*a,-1*b])/np.sqrt(3)
u6 = np.array([-1*a,b,-1*a,-1*b])/np.sqrt(3)
u7 = np.array([1,0,a,-1*b])/np.sqrt(3)
u8 = np.array([a,b,-1*a,-1*b])/np.sqrt(3)
u9 = np.array([-1*a,b,a,-1*b])/np.sqrt(3)

set = [u1,u2,u3,u4,u5,u6,u7,u8,u9]
# check that the u_i are normal vectors

print("magnitude of vectors u1-u9:")
for i in set:
    print("{:.2f}".format(dot(i,i)))


# check that they satisfy equiangular condition
print("overlap of all pairs of lines")
for i in range(len(set)):
    for j in range(i+1,len(set)):
        print("{:.2f}".format(dot(set[i], set[j])))

'''        

# check that the lines resolve the identity
print("synthesis operator: (should be a multiple of the identity)")
s = np.transpose(np.array([u1, u2, u3, u4, u5]))
sc = np.transpose(s)
print(s)
print(np.matmul(s,sc))


# unitary matrices producing the maximally entangled states (matrices 3 and 5 need to be switched. check the states)
mat1 = np.asmatrix(np.array([[1j, -1j], [1j, 1j]])/np.sqrt(2))
mat2 = np.asmatrix(np.array([[d+a*1j, b-c*1j],[b+c*1j, -1*d+a*1j]])/np.sqrt(2))
mat3 = np.asmatrix(np.array([[-1*b+c*1j, d-a*1j], [d+a*1j, b+c*1j]])/np.sqrt(2))
mat4 = np.asmatrix(np.array([[b+c*1j, -1*d-a*1j], [-1*d+a*1j, -1*b+c*1j]])/np.sqrt(2))
mat5 = np.asmatrix(np.array([[-1*d+a*1j, -1*b-c*1j], [-1*b+c*1j, d+a*1j]])/np.sqrt(2))
# checks that the matrices really are unitary matrices
print("check that the matrices are unitary matrices")
print(np.matmul(mat1.getH(), mat1))
print(np.matmul(mat2.getH(), mat2))
print(np.matmul(mat3.getH(), mat3))
print(np.matmul(mat4.getH(), mat4))
print(np.matmul(mat5.getH(), mat5))'''
