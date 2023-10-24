import numpy as np

a = (-1+np.sqrt(5))/4
b = np.sqrt((5+np.sqrt(5))/8)
c = -(1+np.sqrt(5))/4
d = np.sqrt((5-np.sqrt(5))/8)

u1 = np.array([[1.j,-1.j],[1.j,1.j]])
u2 = np.array([[d+a*1.j, b-1.j*c],[b+1.j*c,-1*d+1.j*a]])
u3 = np.array([[-1*d+a*1.j,-b-1.j*c],[-b+1.j*c,d+a*1.j]])
u4 = np. array([[b+c*1.j,-d-a*1.j],[-d+a*1.j,-b+c*1.j]])
u5 = np.array([[-b+c*1.j,d-a*1.j],[d+a*1.j, b+c*1.j]])

u1_inv = 0.5*np.array([[-1.j,-1.j],[1.j,-1.j]])
u2_inv = 0.5*np.array([[d-a*1.j,b-c*1.j],[b+c*1.j,-a*1.j-d]])
u3_inv = 0.5*np.array([[-d-a*1.j,-b-c*1.j],[-b+c*1.j,d-a*1.j]])
u4_inv = 0.5*np.array([[b-c*1.j,-d-a*1.j],[-d+a*1.j,-b-c*1.j]])

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