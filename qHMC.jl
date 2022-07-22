using LinearAlgebra
using SparseArrays
using Arpack
using JLD2
# using PyPlot
using Random, Distributions
Random.seed!(123)
include("GlauberMixingLocal.jl")
include("GlauberMixingGeneral.jl")

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

function ham(alpha, eta, N,couplings, h)
    H = alpha*ising_ham(N, couplings, h) + eta*mixing_ham(N)
    return H
end

function mixing_matrix(N,couplings,h,beta,alpha, eta) 
    H = ham(alpha, eta, N,couplings, h)
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

# N=3
# alpha = 5
# couplings = ones(N)
# eta = 1
# H = ham(alpha, eta, N,couplings, 0)
# U = exp(-1im*H)
# prob = abs.(U).^2
# # display(prob)
# beta = 6
# M = mixing_matrix(N,couplings,0,beta,alpha, eta)
# display(M)
# M_glab = mixing_glab(N,couplings,0,beta)
# display(M_glab)

# beta = 6
# N_values = (2:10)
# alpha = 0
# num_values = 500
# eta_values = range(0,10, length=num_values)
# h=0

# gap_all = zeros(length(N_values), num_values)
# for i in (1:num_values)
#     eta = eta_values[i]
#     println(" Working on i = ",i)
#     for j in (1:length(N_values))
#         N = N_values[j]
#         println(" Working on N = ",N)
#         couplings = ones(N)
#         couplings[end] = 0 
#         M = mixing_matrix(N,couplings,0,beta,alpha, eta)
#         try
#             e,v  = eigs(M, nev = 2, which=:LM)
#             gap_all[j,i] = 1-abs(e[2])
#         catch
#             println("Arpack method out of iteration")
#             e,v  = eigen(M)
#             gap_all[j,i] = 1-abs(e[end-1])
#         end
#     end
# end
# save_object("Data/qHMC/alphaEtaParam/alpha"*string(alpha)*"scalingBeta6", gap_all)
# # save_object("Data/GridSearch/alphaEtaParam/AlphaZero/"*string(num_values)*"N"*string(N)*"beta"*string(beta)*"e", eigens)


# beta = 6
# N_values = (5:12)
# alpha = 0
# epsilon = 0.01
# eta = 2.9
# h=0

# e2 = zeros(length(N_values))
# e3 = zeros(length(N_values))

# for j in (1:length(N_values))
#     N = N_values[j]
#     println("  N = ",N)
#     couplings = ones(N)
#     couplings[end] = 0 
#     M = mixing_matrix(N,couplings,0,beta,alpha, eta)
#     try
#         e,v  = eigs(M, nev = 3, which=:LM)
#         e2[j] = 1-abs(e[2])
#         e3[j] = 1-abs(e[3])
#     catch
#         println("Arpack method out of iteration")
#         e,v  = eigen(M)
#         e2[j] = 1-abs(e[end-1])
#         e3[j] = 1-abs(e[end-2])
#     end
# end

# save_object("Data/qHMC/alphaEtaParam/alpha"*string(alpha)*"eta"*string(eta)*"beta6e2", e2)
# save_object("Data/qHMC/alphaEtaParam/alpha"*string(alpha)*"eta"*string(eta)*"beta6e3", e3)
