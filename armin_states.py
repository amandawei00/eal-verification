import numpy as np

# numabers in Armin's states
av = 0.176777
bv = 0.654329
cv = 0.145375
dv = 0.139536
ev = 0.0289026
fv = 0.451417
gv = 0.513943
hv = 0.327545
kv = 0.573293
mv = 0.181104
nv = 0.297881
pv = 0.26725
qv = 0.55512

# ARmin's states
v1 = np.array([-1, 0, 0, 0])
v2 = np.array([av, dv, cv, bv])*np.sqrt(2)
v3 = np.array([av, -gv, -fv, -ev])*np.sqrt(2)
v4 = np.array([-av, mv, -kv, hv])*np.sqrt(2)
v5 = np.array([av, qv, -pv, -nv])*np.sqrt(2)

# numbers in my states -------------------
aw = (-1+np.sqrt(5))/4
bw = np.sqrt((5+np.sqrt(5))/8)
cw = -(1+np.sqrt(5))/4
dw = np.sqrt((5-np.sqrt(5))/8)

# my states----------------------
w1 = np.array([1,0,1,0])/np.sqrt(2)
w2 = np.array([aw,bw,cw,dw])/np.sqrt(2)
w3 = np.array([aw,-bw,cw,-dw])/np.sqrt(2)
w4 = np.array([cw,-dw,aw,bw])/np.sqrt(2)
w5 = np.array([cw,dw,aw,-bw])/np.sqrt(2)

p = np.array([[np.matmul(v1,w1), np.matmul(v1,w2), np.matmul(v1,w3), np.matmul(v1,w4), np.matmul(v1,w5)],
              [np.matmul(v2,w1), np.matmul(v2,w2), np.matmul(v2,w3), np.matmul(v2,w4), np.matmul(v2,w5)],
              [np.matmul(v3,w1), np.matmul(v3,w2), np.matmul(v3,w3), np.matmul(v3,w4), np.matmul(v3,w5)],
              [np.matmul(v4,w1), np.matmul(v4,w2), np.matmul(v4,w3), np.matmul(v4,w4), np.matmul(v4,w5)],
              [np.matmul(v5,w1), np.matmul(v5,w2), np.matmul(v5,w3), np.matmul(v5,w4), np.matmul(v5,w5)]])
pinv = np.linalg.inv(p)
