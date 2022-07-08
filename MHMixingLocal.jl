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

function mixing(N,couplings,h,beta)
    dim = 2^N
    M = zeros(dim,dim)
    for ket in (0:dim-1)
        E_sp = ising_energy(N,couplings,h,ket)
        diag = 0
        for SpinIndex in (0:N-1)
            bit = Int(2)^(SpinIndex)
            bra = ket ⊻ bit
            E_s = ising_energy(N,couplings,h,bra)
            # Sigma_x term
            m = (1/N)*min(1,exp(-beta*(E_s-E_sp)))
            M[bra+1,ket+1] = m
            diag += (1/N)*min(1,exp(beta*(E_s-ising_energy(N,couplings,h,ket))))
        end 
        M[ket+1,ket+1] = 1-diag
    end
    return M #|> sparse
end


N = 10
h=5
couplings = ones(N)
couplings[end] = 0
temp = 10 .^ (range(-2.5,stop=2.5,length=50))
beta_values = 1 ./ temp
gap_all = zeros(length(beta_values))
for j in (1:length(beta_values))
    beta = beta_values[j]
    println(" Working on beta = ",beta)
    M = mixing(N,couplings,h,beta)
    e,v  = eigs(M, nev = 3, which=:LM)
    gap_all[j] = 1-abs(e[2])
    # e,v  = eigen(M)
    # gap_all[j]= 1-abs(e[end-1])
end
save_object("Data/MHLoc/OBCN10h"*string(h), gap_all)

# beta = 5
# N_values = (5:13)
# h=0
# gap_all = zeros(length(N_values))
# for j in (1:length(N_values))
#     N = N_values[j]
#     couplings = ones(N)
#     # couplings[end] = 0
#     println(" Working on N = ",N)
#     M = mixing(N,couplings,h,beta)
#     e,v  = eigs(M, nev = 2, which=:LM)
#     gap_all[j] = 1-abs(e[2])
#     print("Gap: ",1-abs(e[2]))
# end
# save_object("Data/MHLoc/PBCBeta"*string(beta), gap_all)
# beta = 5
# N_values = (5:13)
# h=0
# gap_all = zeros(length(N_values))
# for j in (1:length(N_values))
#     N = N_values[j]
#     couplings = ones(N)
#     couplings[end] = 0
#     println(" Working on N = ",N)
#     M = mixing(N,couplings,h,beta)
#     e,v  = eigs(M, nev = 2, which=:LM)
#     gap_all[j] = 1-abs(e[2])
#     print("Gap: ",1-abs(e[2]))
# end
# save_object("Data/MHLoc/OBCBeta"*string(beta), gap_all)



