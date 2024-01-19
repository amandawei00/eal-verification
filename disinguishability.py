import numpy as np
from scipy.optimize import

a = np.cos(2*np.pi/5)
b = np.sin(2*np.pi/5)
c = np.cos(4*np.pi/5)
d = np.cos(4*np.pi/5)

comp = np.array([np.array([1,0]), np.array([0,1])])
e0 = comp[0]
e1 = comp[1]

def shift_basis(theta):
    return np.array([np.cos(theta)*e0+np.sin(theta)*e1, -np.cos(theta)*e0+np.sin(theta)*e1])

meredes = np.array([np.array([1,0]), np.array([np.cos(np.pi/3), np.sin(np.pi/3)]),np.array([np.cos(2*np.pi/3), np.sin(2*np.pi/3)])])

# probability of obtaining component e of state state. it is the magnitude squared of the associated coefficient
def p(state, e, bas):
       return np.round(np.conjugate(np.matmul(state, bas[e])) * np.matmul(state, bas[e]),5)

def norm(state):
    return state/np.sqrt(np.matmul(np.conjugate(state),state))

# total probability of observing outcome out given possible states states in basis bas
def sum_prob(states, out, bas):
    sum = 0
    for i in range(len(states)):
        sum += p(states[i], out, bas)
    return np.round((1 / len(states)) * sum,5)

# given outcome out in basis bas, returns array of probabilities that it was the associated entry
# assuming states are handed out in equal distributions
def back_prob(states, out, bas):
    return np.real(np.array([p(states[i], out, bas)*(1/len(states))/sum_prob(states,out,bas) for i in range(len(states))]))


alice = np.array([(e0+e1)/np.sqrt(2), norm(((a-1.j*d)*e0+(c-1.j*b)*e1)/np.sqrt(2)),
                  norm(((a+d*1.j)*e0+(c+1.j*b)*e1)/np.sqrt(2)), norm(((c-1.j*b)*e0+(a+1.j*d)*e1)/np.sqrt(2)),
                  norm(((c+1.j*b)*e0+(a-1.j*d)*e1)/np.sqrt(2))])
bob = np.array([(e0+e1)/np.sqrt(2), norm(((c-1.j*d)*e0+(a-1.j*b)*e1)/np.sqrt(2)),
                  norm(((c+d*1.j)*e0+(a+1.j*b)*e1)/np.sqrt(2)), norm(((a-1.j*b)*e0+(c+1.j*d)*e1)/np.sqrt(2)),
                  norm(((a+1.j*b)*e0+(c-1.j*d)*e1)/np.sqrt(2))])

# probability for alice to measure 0, 1
palice0 = sum_prob(alice, 0, comp)
palice1 = sum_prob(alice,1,comp)

# probability for bob to measure 0, 1
pbob0 = sum_prob(bob,0,comp)
pbob1 = sum_prob(bob,1,comp)


# prob of each state (given in array) given that alice measures outcome out
# prob of each state (given in array) given that bob measures outcome out
prob_arr_alice = back_prob(alice,1, shift_basis(np.pi/4))
prob_arr_bob = back_prob(bob,1, shift_basis(3*np.pi/2))

# dist1 and dist2 are alice and bob's beliefs given their measurements
# i think describing the marginal probability distribution is the correct way to go.
# so Bob's distribution really should be p(state i |alice's measure and bob's measure)


# alice measures hers in basis defined by theta, produces a distribution from this
# bob measures his in a basis defined by phi, produces a distribution
# given the two distributions, alice and bob choose states 1,2,3,4,5 with probability p1,p2,p3,p4,p5 to be optimized
# the number of times they are correct is the distinguishability
def f(theta, phi1, phi2, p0, p1, p2, p3, p4):
    # random state, one of the 5
    state = np.random.choice(5)

    # alice and bob's measurements will land on one of the two basis states, determined randomly by meas
    # fix this. in some bases it may not be 50/50
    meas  = np.array([np.random.choice(2), np.random.choice(2)])

    alice_dist = back_prob(alice,meas[0], shift_basis(theta))
    if meas[0]==0:
        bob_dist = back_prob(bob, meas[1], shift_basis(phi1))
    elif meas[0]==1:
        bob_dist = back_prob(meas[1], shift_basis(phi2))

    # associated to each state is a probability: p(state i| meas[0], meas[1]) = (4/5)*p(meas[0],meas[1]|state i)
    dist = 4 * np.matmul(bob_dist, alice_dist)
    prob_choice = np.array([p0,p1,p2,p3,p4])

    choice = np.random.choice(5, np.matmul(dist,prob_choice))

    if choice == state:
        return 1
    else:
        return 0

def distinguishability(theta,phi1,phi2,p0,p1,p2,p3,p4):
    # n trials
    n = 1000
    count = 0
    for i in range(n):
        count += f(theta, phi1, phi2, p0, p1, p2, p3, p4)

    return count/n

# optimize distinguishability
# chck operation of distinguishability first
