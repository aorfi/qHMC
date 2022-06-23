using LinearAlgebra
using SparseArrays
using Arpack
using JLD2

function sigma_z(N, config, spin)
    # config is in [0,2^N] and spin in [1,N]
    spin_index = spin -1
    return 2*((config>>spin_index)&1)-1
end

function accept_prob_matrix(N, beta, config,spin)
    spin_forward = spin + 1
    spin_back = spin - 1
    #PBC
    if spin == N
        spin_forward = 1
    end
    if spin == 1
        spin_back = N
    end
    p_top = exp(-beta*sigma_z(N,config,spin)*(sigma_z(N,config,spin_forward)+sigma_z(N,config,spin_back)))
    p_bot = exp(-beta*(sigma_z(N,config,spin_forward)+sigma_z(N,config,spin_back)))+exp(beta*(sigma_z(N,config,spin_forward)+sigma_z(N,config,spin_back)))
    return p_top/p_bot
end

function mixing_local(N, beta)
    dim = (2)^N
    M = zeros(dim,dim)
    for ket in (0:dim-1)
        ket_binary = bitstring(ket)
        Diagonal = Int64(0)
        p_sum = 0
        for SpinIndex in (0:N-1)
            bit = Int(2)^(SpinIndex)
            bra = ket ⊻ bit
            # Sigma_x term
            M[bra+1,ket+1] = (1/N)*accept_prob(N, beta, bra,SpinIndex+1)
            p_sum += (1/N)*accept_prob(N, beta, ket,SpinIndex+1)
        end
        # Diagonal term
        M[ket+1,ket+1] += 1-p_sum
    end
    return M |> sparse
end


# beta = 6
# N_values = (2:15)
# N_max = last(N_values)
# gap_all = zeros(length(N_values))
# for j in (1:length(N_values))
#     N = N_values[j]
#     println(" Working on N = ",N)
#     M = mixing(N,beta)
#     e,v  = eigs(M, nev = 2, which=:LR)
#     gap_all[j] = abs(e[1]-e[2])
# end
# save_object("Data/gapBeta6", gap_all)

# N = 10
# temp = 10 .^ (range(-2.5,stop=2.5,length=50))
# beta_values = 1 ./ temp
# gap_all = zeros(length(beta_values))
# for j in (1:length(beta_values))
#     beta = beta_values[j]
#     println(" Working on beta = ",beta)
#     M = mixing(N,beta)
#     e,v  = eigs(M, nev = 2, which=:LR)
#     gap_all[j] = abs(e[1]-e[2])
# end
# save_object("Data/gapN10", gap_all)