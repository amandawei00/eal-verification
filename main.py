import numpy as np

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
            print("vectors " + str(i) + ", " + str(j) + ": " + str(round(np.abs(dot(a[i], a[j])), 5)))

    print("yes equiangular")

# input: set of vectors
def identity_resolution(a):
    mat = np.asmatrix(a)
    return np.matmul(mat.getH(), mat).round(4)

# input: generator of the group (r), and a fiducial vector (fid)
def set_generator(r, fid):
    v1 = fid
    set = [v1]
    # condition that the application of g doesn't complete the cycle
    while not (np.matmul(r, fid).round(8) == v1.round(8)).all():
        new_vec = np.matmul(r, fid)
        set.append(new_vec)
        fid = new_vec

    return set

# returns set of tuples, coordinates of roots of unity coord. on unit circle
def nth_roots(n):
    return 0 #(not done)
# test of equiangular lines
'''a = ((1,0,0),(0,1,0),(0,0,1)) # orthoggonal basis in d=3
# a = ((1,0),(0.5, np.sqrt(3)/2),(-0.5, np.sqrt(3)/2)) # set of equiangular lines in d=2

# corners of the icosahedron, maximal equiangular lines in d=3
# p = 0.5*(1 + np.sqrt(5))

# a = ((0, 1, p), (0, 1, -1 * p), (-1, p, 0), (-1, -p, 0), (p, 0, 1), (-1 * p, 0, 1))
# equiangular(a)'''

# 5 lines in R4 -------------------------------------------------
a = (-1+np.sqrt(5))/4
b = np.sqrt((5+np.sqrt(5))/8)
c = -(1+np.sqrt(5))/4
d = np.sqrt((5-np.sqrt(5))/8)
r5_generator = np.array([[a, -1*b, 0, 0],
               [b, a, 0, 0],
               [0, 0, c, -1*d],
               [0, 0, d, c]])
r5_fid = np.array([1, 0, 1, 0])/np.sqrt(2)

# 7 lines in R6 -------------------------------------------------

p = np.cos(2*np.pi/7)
q = np.sin(2*np.pi/7)
r = np.cos(4*np.pi/7)
s = np.sin(4*np.pi/7)
t = np.cos(6*np.pi/7)
u = np.sin(6*np.pi/7)

r7_generator = np.array([[p, -1*q, 0, 0, 0, 0],
                         [q, p, 0, 0, 0, 0],
                         [0, 0, r, -1*s, 0, 0],
                         [0, 0, s, r, 0, 0],
                         [0, 0, 0, 0, t, -1*u],
                         [0, 0, 0, 0, u, t]])
r7_fid = np.array([1, 0, 1, 0, 1, 0])/np.sqrt(3)

# 17 lines in R16 -------------------------------------------------
ax = np.cos(2*np.pi/17)
ay = np.sin(2*np.pi/17)
bx = np.cos(4*np.pi/17)
by = np.sin(4*np.pi/17)
cx = np.cos(6*np.pi/17)
cy = np.sin(6*np.pi/17)
dx = np.cos(8*np.pi/17)
dy = np.sin(8*np.pi/17)
ex = np.cos(10*np.pi/17)
ey = np.sin(10*np.pi/17)
fx = np.cos(12*np.pi/17)
fy = np.sin(12*np.pi/17)
gx = np.cos(14*np.pi/17)
gy = np.sin(14*np.pi/17)
hx = np.cos(16*np.pi/17)
hy = np.sin(16*np.pi/17)
r17_generator = np.array([[ax,-1*ay,0,0, 0,0,0,0, 0,0,0,0, 0,0,0,0],
                          [ay, ax, 0, 0, 0,0,0,0, 0,0,0,0, 0,0,0,0],
                          [0,0,bx,-1*by, 0,0,0,0, 0,0,0,0, 0,0,0,0],
                          [0, 0, by, bx, 0,0,0,0, 0,0,0,0, 0,0,0,0],
                          [0,0,0,0, cx,-1*cy,0,0, 0,0,0,0, 0,0,0,0],
                          [0,0,0,0, cy, cx, 0, 0, 0,0,0,0, 0,0,0,0],
                          [0,0,0,0, 0,0,dx,-1*dy, 0,0,0,0, 0,0,0,0],
                          [0,0,0,0, 0, 0, dy, dx, 0,0,0,0, 0,0,0,0],
                          [0,0,0,0, 0,0,0,0, ex,-1*ey,0,0, 0,0,0,0],
                          [0,0,0,0, 0,0,0,0, ey, ex, 0, 0, 0,0,0,0],
                          [0,0,0,0, 0,0,0,0, 0,0,fx,-1*fy, 0,0,0,0],
                          [0,0,0,0, 0,0,0,0, 0, 0, fy, fx, 0,0,0,0],
                          [0,0,0,0, 0,0,0,0, 0,0,0,0, gx,-1*gy,0,0],
                          [0,0,0,0, 0,0,0,0, 0,0,0,0, gy, gx, 0, 0],
                          [0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,hx,-1*hy],
                          [0,0,0,0, 0,0,0,0, 0,0,0,0, 0,0,hy,hx]])
r17_fid = np.array([1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0])/np.sqrt(8)

set = set_generator(r7_generator, r7_fid)

# check that the u_i are normal vectors

print("magnitude of vectors:")
for i in set:
    print("{:.4f}".format(np.matmul(i, i)))


# check that they satisfy equiangular condition
print("overlap of all pairs of lines")
for i in range(len(set)):
    for j in range(i+1, len(set)):
        print("("+str(i)+","+str(j)+"): "+"{:.4f}".format(np.matmul(set[i], set[j])))

# check that they resolve the identity
print(identity_resolution(set))