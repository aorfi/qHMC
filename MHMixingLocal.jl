using LinearAlgebra
using SparseArrays
using Arpack
using JLD2

function ising_energy(N, config,h)
    # config is in [0,2^N] and spin in [1,N]
    eng = 0
    for SpinIndex in (0:N-1)
        # PBC
        if SpinIndex == N-1
            Si = 2*((config>>SpinIndex)&1)-1
            eng += h*Si
            Si_next = 2*((config>>(1))&1)-1
            eng += -Si*Si_next
            break
        end
        # OBC 
        # if SpinIndex == N-1
        #     Si = 2*((config>>SpinIndex)&1)-1
        #     eng += h*Si
        #     break
        # end
        Si = 2*((config>>SpinIndex)&1)-1
        eng += h*Si
        Si_next = 2*((config>>(SpinIndex+1))&1)-1
        eng += -Si*Si_next
    end
    return eng
end

function mixing(N, beta,h)
    dim = (2)^N
    M = zeros(dim,dim)
    for ket in (0:dim-1)
        E_sp = ising_energy(N,ket,h)
        diag = 0
        for SpinIndex in (0:N-1)
            bit = Int(2)^(SpinIndex)
            bra = ket ⊻ bit
            E_s = ising_energy(N,bra,h)
            # Sigma_x term
            m = (1/N)*min(1,exp(-beta*(E_s-E_sp)))
            M[bra+1,ket+1] = m
            diag += (1/N)*min(1,exp(beta*(E_s-ising_energy(N,ket,h))))
        end 
        M[ket+1,ket+1] = 1-diag
    end
    return M |> sparse
end




# beta = 6
# N_values = (2:13)
# h=0
# gap_all = zeros(length(N_values))
# for j in (1:length(N_values))
#     N = N_values[j]
#     println(" Working on N = ",N)
#     M = mixing(N,beta,h)
#     e,v  = eigs(M, nev = 3, which=:LR)
#     gap_all[j] = abs(1-e[2])
# end
# save_object("Data/MHLoc/MHLocPBCgapBeta6", gap_all)

N = 12
h=0
temp = 10 .^ (range(-2.5,stop=2.5,length=50))
beta_values = 1 ./ temp
gap_all = zeros(length(beta_values))
for j in (1:length(beta_values))
    beta = beta_values[j]
    println(" Working on beta = ",beta)
    M = mixing(N,beta,h)
    e,v  = eigs(M, nev = 3, which=:LR)
    gap_all[j] = abs(1-e[2])
    # gap_all[j] = abs(1-e[3])
end
save_object("Data/MHLoc/MHLocGapOBCN12", gap_all)


# N = 10
# h = 0.1
# temp = 10 .^ (range(-2.5,stop=2.5,length=50))
# beta_values = 1 ./ temp
# gap_all = zeros(length(beta_values))
# e1 = zeros(length(beta_values))
# e2 = zeros(length(beta_values))
# e3 = zeros(length(beta_values))
# for j in (1:length(beta_values))
#     beta = beta_values[j]
#     println(" Working on beta = ",beta)
#     M = mixing_field(N,beta,h)
#     e,v  = eigs(M, nev = 3, which=:LR)
#     e1[j] = abs(e[1])
#     e2[j] = abs(e[2])
#     e3[j] = abs(e[3])
#     # print("\n largest ", e[1])
#     # print(" second ", e[2])
#     # print(" third ",e[3])
#     # gap_all[j] = abs(1-e[2])
# end
# save_object("Data/MHLoc/MHLocLocF0.1OBCgapN10e1", e1)
# save_object("Data/MHLoc/MHLocLocF0.1OBCgapN10e2", e2)
# save_object("Data/MHLoc/MHLocLocF0.1OBCgapN10e3", e3)