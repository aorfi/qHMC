using LinearAlgebra
using SparseArrays
using Arpack
using JLD2

function ising_energy(N, config,h)
    # config is in [0,2^N] and spin in [1,N]
    eng = 0
    for SpinIndex in (0:N-1)
        # PBC
        # if SpinIndex == N-1
        #     Si = 2*((config>>SpinIndex)&1)-1
        #     eng += h*Si
        #     Si_next = 2*((config>>(1))&1)-1
        #     eng += -Si*Si_next
        #     break
        # end
        # OBC 
        if SpinIndex == N-1
            Si = 2*((config>>SpinIndex)&1)-1
            eng += h*Si
            break
        end
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
    for col in (0:dim-1)
        E_sp = ising_energy(N,col,h)
        diag = 0 
        for row in (0:dim-1)
            if row != col
                E_s = ising_energy(N,row,h)
                m = (1/2^N)*min(1,exp(-beta*(E_sp-E_s)))
                M[col+1,row+1] = m
                diag += m
            end
        M[col+1,col+1] = 1-diag
        end
    end
    return M
end


# beta = 300
# N_values = (2:13)
# h=0
# gap_all = zeros(length(N_values))
# for j in (1:length(N_values))
#     N = N_values[j]
#     println(" Working on N = ",N)
#     M = mixing(N,beta,h)
#     e,v  = eigen(M)
#     gap_all[j] = abs(1-e[end-1])
#     print("Gap: ",abs(1-e[end-1]))
# end
# save_object("Data/MH/OBCMHgapBeta300", gap_all)


# N = 10
# temp = 10 .^ (range(-2.5,stop=2.5,length=50))
# beta_values = 1 ./ temp
# gap_all = zeros(length(beta_values))
# for j in (1:length(beta_values))
#     beta = beta_values[j]
#     println(" Working on beta = ",beta)
#     M = mixing(N,beta)
#     e,v  = eigs(M, nev = 3, which=:LR)
#     gap_all[j] = abs(e[1]-e[2])
#     # gap_all[j] = abs(1-e[3])
# end
# save_object("Data/MH/MHgapOBCN10", gap_all)