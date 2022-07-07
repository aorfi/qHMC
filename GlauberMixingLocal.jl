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
        p_sum = 0
        for SpinIndex in (0:N-1)
            bit = Int(2)^(SpinIndex)
            bra = ket ⊻ bit
            E_s = ising_energy(N,couplings,h,bra)
            # Sigma_x term
            M[bra+1,ket+1] = (1/N)*(1/(1+exp(beta*(E_s-E_sp))))
            p_sum += (1/N)*accept_prob(N, beta, ket,SpinIndex+1,h)
        end
        # Diagonal term
        M[ket+1,ket+1] += 1-p_sum
    end
    return M |> sparse
end



beta = 6
N = 3
h=0
couplings = ones(N)
couplings[end] = 0
M = mixing(N,couplings,h,beta)
display(M)
display(sum(M[1,:]))
# e,v  = eigs(M, nev = 3, which=:LM)
# display(e)

# beta = 6
# h=0
# N_values = (5:15)
# N_max = last(N_values)
# gap_all = zeros(length(N_values))
# for j in (1:length(N_values))
#     N = N_values[j]
#     println(" Working on N = ",N)
#     M = mixing(N,beta,h)
#     e,v  = eigs(M, nev = 3, which=:LM)
#     # gap_all[j] = abs(1-e[3])
#     gap_all[j] = abs(1-e[2])
# end
# save_object("Data/GlaubLoc/CmpareGlaubLocOBCgapBeta6", gap_all)

# N = 12
# h = 0
# temp = 10 .^ (range(-2.5,stop=2.5,length=50))
# beta_values = 1 ./ temp
# gap_all = zeros(length(beta_values))
# e1 = zeros(length(beta_values))
# e2 = zeros(length(beta_values))
# e3 = zeros(length(beta_values))
# for j in (1:length(beta_values))
#     beta = beta_values[j]
#     println(" Working on beta = ",beta)
#     M = mixing(N,beta,h)
#     e,v  = eigs(M, nev = 3, which=:LM)
#     e1[j] = abs(e[1])
#     e2[j] = abs(e[2])
#     e3[j] = abs(e[3])
#     # print("\n largest ", e[1])
#     # print(" second ", e[2])
#     # print(" third ",e[3])
#     gap_all[j] = abs(1-e[2])
# end
# save_object("Data/GlaubLoc/GlaubLocOBCgapN12", gap_all)
# save_object("Data/GlaubLoc/GlaubLocF0.1OBCgapN10e1", e1)
# save_object("Data/GlaubLoc/GlaubLocF0.1OBCgapN10e2", e2)
# save_object("Data/GlaubLoc/GlaubLocF0.1OBCgapN10e3", e3)