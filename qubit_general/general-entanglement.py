import numpy as np
from scipy.optimize import minimize

a = np.cos(2*np.pi/5)
b = np.sin(2*np.pi/5)
c = np.cos(4*np.pi/5)
d = np.sin(4*np.pi/5)

theta = 0
b0 = np.array([1, 0, 0, 0])
b1 = np.array([0, 0, 0, 1])

paulix = np.array([[0,1],[1,0]])
pauliy = np.array([[0,-1.j],[1.j,0]])
pauliz = np.array([[1,0],[0,-1]])
id = np.array([[1,0],[0,1]])

# bloch vectors
v1 = np.array([-b-d,a-c,b-d])
v2 = np.array([b+d,a-c,-b+d])
v3 = np.array([-b+d,-a+c,-b-d])
v4 = np.array([b-d,-a+c,b+d])

v1 = v1/np.sqrt(np.matmul(v1,v1))
v2 = v2/np.sqrt(np.matmul(v2,v2))
v3 = v3/np.sqrt(np.matmul(v3,v3))
v4 = v4/np.sqrt(np.matmul(v4,v4))

# given angles and bloch vecotr, produces ui tensor vi matrices
def uivi(alpha,beta,v):
    u1 = np.cos(alpha)*id + 1.j*np.sin(alpha)*(v[0]*paulix+v[1]*pauliy+v[2]*pauliz)
    u2 = np.cos(beta)*id + 1.j*np.sin(beta)*(v[0]*paulix+v[1]*pauliy+v[2]*pauliz)

    return np.tensordot(u1, u2,0).transpose((0, 2, 1, 3)).reshape((4, 4))

# returns |\langle \phi i | \phi j \rangle|^2 given ui\tensor vi and uj\tensor vj
def prod(ui,uj):
    uiuj = np.asarray(np.matmul(np.asmatrix(ui).getH(), uj))
    return np.power(np.abs(0.5*(np.power(np.cos(theta),2) * (np.matmul(np.transpose(b0), np.matmul(uiuj, b0)))
                + np.power(np.sin(theta), 2) * (np.matmul(np.transpose(b1), np.matmul(uiuj, b1)))
                + np.cos(theta)*np.sin(theta) * (np.matmul(np.transpose(b0),np.matmul(uiuj,b1)))
                                                           + np.matmul(np.transpose(b1),np.matmul(uiuj,b0)))),2)

# x0: a1(0), a2(1), a3(2), a4(3), b1(4), b2(5), b3(6), b4(7)
#x0 = [th1,2*th1,3*th1,4*th1,th1,2*th1,3*th1,4*th1]

def mu(u):
    return np.power(np.abs(prod(uivi(u[0], u[4], v1), uivi(u[1], u[5], v2))), 2)


n=100
min=None
for i in range(n):
    x0 = x0 = np.random.rand(8)
    res = minimize(mu, x0=x0, bounds=((0,2*np.pi),(0,2*np.pi),(0,2*np.pi),(0,2*np.pi),(0,2*np.pi),(0,2*np.pi),(0,2*np.pi),(0,2*np.pi)),
                                constraints=({'type':'eq','fun': lambda x: np.power(np.abs(prod(uivi(x[0],x[4],v1), np.identity(4))), 2) - np.power(np.abs(prod(uivi(x[0], x[4],v1), uivi(x[1],x[5], v2))), 2)},
                               {'type':'eq','fun': lambda x: np.power(np.abs(prod(uivi(x[1],x[5],v2), np.identity(4))), 2) - np.power(np.abs(prod(uivi(x[0], x[4],v1), uivi(x[1],x[5], v2))), 2)},
                               {'type':'eq','fun': lambda x: np.power(np.abs(prod(uivi(x[2],x[6],v3), np.identity(4))), 2) - np.power(np.abs(prod(uivi(x[0], x[4],v1), uivi(x[1],x[5], v2))), 2)},
                               {'type':'eq','fun': lambda x: np.power(np.abs(prod(uivi(x[3],x[7],v4), np.identity(4))), 2) - np.power(np.abs(prod(uivi(x[0], x[4],v1), uivi(x[1],x[5], v2))), 2)},
                               {'type':'eq','fun': lambda x: np.power(np.abs(prod(uivi(x[0],x[4],v1), uivi(x[1],x[5],v2))),2) - np.power(np.abs(prod(uivi(x[0], x[4],v1), uivi(x[1],x[5], v2))), 2)},
                               {'type':'eq','fun': lambda x: np.power(np.abs(prod(uivi(x[0],x[4],v1), uivi(x[2],x[6],v3))),2) - np.power(np.abs(prod(uivi(x[0], x[4],v1), uivi(x[1],x[5], v2))), 2)},
                               {'type':'eq','fun': lambda x: np.power(np.abs(prod(uivi(x[0],x[4],v1), uivi(x[3],x[7],v4))),2) - np.power(np.abs(prod(uivi(x[0], x[4],v1), uivi(x[1],x[5], v2))), 2)},
                               {'type':'eq','fun': lambda x: np.power(np.abs(prod(uivi(x[1],x[5],v2), uivi(x[2],x[6],v3))), 2) - np.power(np.abs(prod(uivi(x[0], x[4],v1), uivi(x[1],x[5], v2))), 2)},
                               {'type':'eq','fun': lambda x: np.power(np.abs(prod(uivi(x[1],x[5],v2), uivi(x[3],x[7],v4))),2) - np.power(np.abs(prod(uivi(x[0], x[4],v1), uivi(x[1],x[5], v2))), 2)},
                               {'type':'eq','fun': lambda x: np.power(np.abs(prod(uivi(x[2],x[6],v3), uivi(x[3],x[7],v4))), 2) - np.power(np.abs(prod(uivi(x[0], x[4],v1), uivi(x[1],x[5], v2))), 2)}))
    if min == None:
        min = res
    elif res.fun < min.fun:
        min = res

print(min)





