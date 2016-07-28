#python script um den gradienten der Abstandsfunktion zu berechnen. Die Abstandsfunktion bezieht sich auf die Punkte (den nächsten) und nicht die Linie


import numpy as np
import matplotlib.pyplot as plt
import math

# Number of legs
N = 4

# Number of measurements
M = 10
# Measurement points
rc = [(np.random.rand(), np.random.rand()) for i in range(M)]

# Plot measurement points
plt.plot([e[0] for e in rc], [e[1] for e in rc], "+")

# Leg length
ell = 0.3

# Regularisation weights
W = [0.02, 0.2, 0.2, 2]

def get_q(x, j):
    q_x = x[0]
    q_y = x[1]

    angle_total = 0.0    
    
    for k in range(j):
        angle_total += x[k + 2]
        q_x += ell*math.cos(angle_total)
        q_y += ell*math.sin(angle_total)
        
    return q_x, q_y
    
    
def get_dist(a, b):
    return (a[0] - b[0])**2 + (a[1] - b[1])**2


def get_nearest_node_idx(x, r):
    return np.argmin([get_dist(get_q(x, j), r) for j in range(N + 1)])


def get_obj_fun(x, r, z):
    qz = [get_q(x, j) for j in z]
    omega = np.sum([0.5*weight*phi**2 for phi, weight in zip(x[2:], W)])
    return np.sum([0.5*get_dist(qz[i], r[i]) for i in range(M)]) + omega
   

def get_exact_grad(x, r, z):
    grad = [0.0 for j in range(len(x))]
    
    # Prepare cum angles, cum sin and cum cos    
    cum_angle = [x[2]]
    for phi in x[3:]:
        cum_angle.append(cum_angle[-1] + phi)
    cum_sin = [-ell*math.sin(a) for a in cum_angle]
    cum_cos = [ ell*math.cos(a) for a in cum_angle]
    
    qz = [get_q(x, j) for j in z]
    dqrx = [qz[i][0] - r[i][0] for i in range(M)]
    dqry = [qz[i][1] - r[i][1] for i in range(M)]
    
    for i in range(M):
        # d/dx0
        grad[0] += dqrx[i]
        # d/dy0
        grad[1] += dqry[i]
        for j in range(min(len(x) - 2, z[i])):
            sin_acc = sum([cum_sin[k] for k in range(j, z[i])])
            cos_acc = sum([cum_cos[k] for k in range(j, z[i])])
            grad[2 + j] += dqrx[i]*sin_acc + dqry[i]*cos_acc
    
    # Regularisation gradient
    for j in range(len(x) - 2):
        grad[2 + j] += x[2 + j]*W[j]

    return grad
   
   
def get_approximate_diff(fun, x):
    eps_rel = 1.0e-7
    eps_abs = 1.0e-9
    grad = []
    for j in range(len(x)):
        if abs(x[j]*eps_rel) < eps_abs:
            delta = eps_abs
        else:
            delta = eps_rel*x[j]
        x_new = [e for e in x]
        x_new[j] += delta
        grad.append((fun(x_new) - fun(x))/delta)
    return grad


xc = [0.0, 0.0, 0.3, -0.3, 0.3, 0.3]
zc = [get_nearest_node_idx(xc, e) for e in rc]
alpha = 2.0e-2
for i in range(50):
    grad_approx = get_approximate_diff(lambda a: get_obj_fun(a, rc, zc), xc)
    grad_exact = get_exact_grad(xc, rc, zc)
    print(sum([abs(approx - exact) for approx, exact in zip(grad_approx, grad_exact)]))
    
    xc = [e - alpha*g for e, g in zip(xc, grad_approx)]
    qc = [get_q(xc, j) for j in range(N + 1)]
    plt.plot([e[0] for e in qc], [e[1] for e in qc], "-o")

