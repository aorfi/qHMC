using LinearAlgebra
using SparseArrays
using Arpack
using JLD2

function ising_energy(N, couplings,h , config)
    # config is in [0,2^N] and spin in [1,N]
    eng = 0
    for SpinIndex in (0:N-1)
        if SpinIndex == N-1
            Si = 2*((config>>SpinIndex)&1)-1
            eng += h*Si
            Si_next = 2*((config>>(0))&1)-1
            eng += -couplings[N]*Si*Si_next
            break
        end
        Si = 2*((config>>SpinIndex)&1)-1
        eng += h*Si
        Si_next = 2*((config>>(SpinIndex+1))&1)-1
        eng += -couplings[N-(SpinIndex+1)]*Si*Si_next
    end
    return eng
end

function ising_ham(N, couplings, h)
    dim = 2^N
    H = zeros(dim,dim)
    for ket in (0:dim-1)
        Diagonal = 0
        for SpinIndex in (0:N-1)
            # Boundary 
            if SpinIndex == N-1
                Si = 2*((ket>>SpinIndex)&1)-1
                Diagonal += h*Si 
                Si_next = 2*((ket>>(0))&1)-1
                Diagonal += -couplings[N]*Si*Si_next
                break
            end
            Si = 2*((ket>>SpinIndex)&1)-1
            Diagonal += h*Si 
            Si_next = 2*((ket>>(SpinIndex+1))&1)-1
            Diagonal += -couplings[N-(SpinIndex+1)]*Si*Si_next
        end
        H[ket+1,ket+1] += Diagonal
    end
    return H #|> sparse
end 


function mixing_ham(N)
    dim = (2)^N
    Hm = zeros(dim,dim)
    for ket in (0:dim-1)
        ket_binary = bitstring(ket)
        for SpinIndex in (0:N-1)
            bit = Int(2)^(SpinIndex)
            bra = ket ⊻ bit
            # Sigma_x term
            Hm[bra+1,ket+1] += 1 
        end
    end
    return Hm
end

function ham(alpha, eta, N)
    H = alpha*ising_ham(N) + eta*mixing_ham(N)
    return H
end

function mixing_matrix(N,beta,alpha, eta)
    H = ham(alpha, eta, N)
    U = exp(-1im*H)
    prob = (abs.(U)).^2
    dim = 2^N
    M = zeros(dim,dim)
    for col in (0:dim-1)
        E_sp = ising_energy(N, couplings,h , col)
        diag = 0 
        for row in (0:dim-1)
            if row != col
                E_s = ising_energy(N, couplings,h,row)
                # MH
                # m = prob[col+1,row+1]*min(1,exp(-beta*(E_sp-E_s)))
                # Glaub
                m = prob[col+1,row+1]*(1/(1+exp(-beta*(E_s-E_sp))))
                M[col+1,row+1] = m
                diag += m
            end
        M[col+1,col+1] = 1-diag
        end
    end
    return M 
end


N = 3
alpha = 0 
eta = 1

couplings = ones(N)

# couplings[end] = 0 # OBC
display(couplings)
H = ising_ham(N,couplings)
display(H)
display(ising_energy(N, couplings,6))




# beta = 6
# N_values = (5:12)
# alpha = 1
# eta = 0.3
# N_max = last(N_values)
# gap_all = zeros(length(N_values))
# for j in (1:length(N_values))
#     N = N_values[j]
#     println(" Working on N = ",N)
#     M = mixing_matrix(N,beta,alpha, eta)
#     e,v  = eigs(M, nev = 3, which=:LM)
#     # gap_all[j] = abs(1-e[3])
#     gap_all[j] = 1-abs(e[2])
# end
# save_object("Data/qHMC/alphaEtaParam/alpha1eta0.3beta6", gap_all)