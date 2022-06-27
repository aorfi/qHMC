import matplotlib as plt
import numpy as np
from qutip import *
import itertools

def qutip_tensors(N):
    # Define sigma matrices
    si = qeye(2)
    sx = sigmax()
    sy = sigmay()
    sz = sigmaz()

    # Create list of sigma matrices acting on one spin
    sx_list = []
    sy_list = []
    sz_list = []
    for n in range(N):
        op_list = []
        for m in range(N):
            op_list.append(si)
        op_list[n] = sx
        sx_list.append(tensor(op_list))
        op_list[n] = sy
        sy_list.append(tensor(op_list))
        op_list[n] = sz
        sz_list.append(tensor(op_list))
        op_list[n] = si
    I = tensor(op_list)
    return sx_list, sy_list, sz_list, I

def ham_classical(N):
    sz_list = qutip_tensors(N)[2]
    Hc =0
    for n in range(N-1):
        Hc += -sz_list[n]*sz_list[n+1] 
    #Periodic boundary conditions
    Hc += -sz_list[N-1]*sz_list[0]
    return Hc

def mixing_prob(N,beta):
    # Can't divide qobjs so need to turn them into np arrays
    sz_list = qutip_tensors(N)[2]
    I = qutip_tensors(N)[3]
    p = []
    for i in range(1,N-1):
        cosh_term = (Qobj.expm(beta*(sz_list[i-1]+sz_list[i+1]))+Qobj.expm(-beta*(sz_list[i-1]+sz_list[i+1]))).diag()
        exp_term = (Qobj.expm(-beta*sz_list[i]*(sz_list[i-1]+sz_list[i+1]))).diag()
        p.append(exp_term/(N*cosh_term))

    # # Periodic boundary conditions
    # # First element
    # cosh_term0 = (Qobj.expm(beta*(sz_list[-1]+sz_list[1]))+Qobj.expm(-beta*(sz_list[-1]+sz_list[1]))).diag()
    # exp_term0 = (Qobj.expm(-beta*sz_list[0]*(sz_list[-1]+sz_list[1]))).diag()
    # p0 = exp_term0/(N*cosh_term0)
    # p.insert(0,p0)
    # # Last element
    # cosh_term_last = (Qobj.expm(beta*(sz_list[-2]+sz_list[0]))+Qobj.expm(-beta*(sz_list[-2]+sz_list[0]))).diag()
    # exp_term_last = (Qobj.expm(-beta*sz_list[-1]*(sz_list[-2]+sz_list[0]))).diag()
    # p_last = exp_term_last / (N*cosh_term_last)
    # p.append(p_last)
    # Change back to qobj 

    # Open boundary conditions
    # First element
    cosh_term0 = (Qobj.expm(beta*(sz_list[1]))+Qobj.expm(-beta*(sz_list[1]))).diag()
    exp_term0 = (Qobj.expm(-beta*sz_list[0]*(sz_list[1]))).diag()
    p0 = exp_term0/(N*cosh_term0)
    p.insert(0,p0)
    # Last element
    cosh_term_last = (Qobj.expm(beta*(sz_list[-2]))+Qobj.expm(-beta*(sz_list[-2]))).diag()
    exp_term_last = (Qobj.expm(-beta*sz_list[-1]*(sz_list[-2]))).diag()
    p_last = exp_term_last / (N*cosh_term_last)
    p.append(p_last)
    # Change back to qobj 
    pQobj = []
    for i in range(len(p)):
        pQ= Qobj(np.eye(2**N)*p[i], dims = 2*[N*[2]])
        pQobj.append(pQ)
    return pQobj

def glauber_mixing(N,p):
    sx_list = qutip_tensors(N)[0]
    # define mixing matrix
    M = 0
    I = qutip_tensors(N)[3]
    for n in range(N):
        M += p[n]*sx_list[n]
    M= M
    # diagonal terms 
    p_sum = sum(p)
    M += (1-p_sum)*I
    return M

def parent_ham(Hc,M,beta,N):
    I = qutip_tensors(N)[3]
    Hp = N*(I-Qobj.expm(-beta*Hc/2)*M*Qobj.expm(beta*Hc/2))
    return Hp


beta = 1
N = 3
sz_list = qutip_tensors(N)[2]
H = 0 
for n in range(N-1):
    H += -sz_list[n]*sz_list[n+1]
print(H)
# sx_list = qutip_tensors(N)[0]
# p = mixing_prob(N,beta)
# M = 0
# for n in range(N):
#     M += p[n]*sx_list[n]

# # M = glauber_mixing(N,p)
# print(M)
# Hp = parent_ham(H,M,beta,N)
# print(Hp)


