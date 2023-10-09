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

# u1 = np.array([1,0,0,0])
# u2 = np.array([1/4, np.sqrt(15)/4, 0,0])
# u3 = np.array([1/4, np.sqrt(3)/(4*np.sqrt(5)), 3/np.sqrt(10),0])
# u4 = np.array([1/4, np.sqrt(15)/20, 0.5*np.sqrt(1/10), np.sqrt(7/8)])
# u5 = np.array([1/4, 3/(4*np.sqrt(15)), np.sqrt(10)/20, np.sqrt(2/7)/4])

a = (np.sqrt(5)-1)/4
b = np.sqrt(1-a*a)
c = -1*(1+np.sqrt(5))/4
d = np.sqrt(1-c*c)

u1 = np.array([1,0,1,0])/np.sqrt(2)
u2 = np.array([a, b, c, d])/np.sqrt(2)
u3 = np.array([a, -1*b, c, -1*d])/np.sqrt(2)
u4 = np.array([c, -1*d, a, b])/np.sqrt(2)
u5 = np.array([c, d, a, -1*b])/np.sqrt(2)
# check that the u_i are normal vectors
print("magnitude of vectors u1-u5:")
print("{:.2f}".format(dot(u1,u1)))
print("{:.2f}".format(dot(u2,u2)))
print("{:.2f}".format(dot(u3,u3)))
print("{:.2f}".format(dot(u4,u4)))
print("{:.2f}".format(dot(u5,u5)))

# check that they satisfy equiangular condition
print("overlap of all pairs of lines")
print("{:.2f}".format(dot(u1, u2)))
print("{:.2f}".format(dot(u1, u3)))
print("{:.2f}".format(dot(u1, u4)))
print("{:.2f}".format(dot(u1, u5)))

print("{:.2f}".format(dot(u2, u3)))
print("{:.2f}".format(dot(u2, u4)))
print("{:.2f}".format(dot(u2, u5)))

print("{:.2f}".format(dot(u3, u4)))
print("{:.2f}".format(dot(u3, u5)))

print("{:.2f}".format(dot(u4, u5)))

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
print(np.matmul(mat5.getH(), mat5))
