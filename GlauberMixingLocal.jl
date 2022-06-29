using LinearAlgebra
using SparseArrays
using Arpack
using JLD2

function sigma_z(N, config, spin)
    # config is in [0,2^N] and spin in [1,N]
    spin_index = spin -1
    return 2*((config>>spin_index)&1)-1
end

function accept_prob(N, beta, config,spin,h)
    spin_forward = spin + 1
    spin_back = spin - 1
    #OBC
    if spin == N
        p_top = exp(beta*sigma_z(N,config,spin)*(sigma_z(N,config,spin_back))-beta*h*sigma_z(N,config,spin))
        p_bot = exp(-beta*(sigma_z(N,config,spin_back)-h))+exp(beta*(sigma_z(N,config,spin_back)-h))
    elseif spin == 1
        p_top = exp(beta*sigma_z(N,config,spin)*(sigma_z(N,config,spin_forward))-beta*h*sigma_z(N,config,spin))
        p_bot = exp(-beta*(sigma_z(N,config,spin_forward)-h))+exp(beta*(sigma_z(N,config,spin_forward)-h))
    else
        p_top = exp(beta*sigma_z(N,config,spin)*(sigma_z(N,config,spin_forward)+sigma_z(N,config,spin_back))-beta*h*sigma_z(N,config,spin))
        p_bot = exp(-beta*(sigma_z(N,config,spin_forward)+sigma_z(N,config,spin_back)-h))+exp(beta*(sigma_z(N,config,spin_forward)+sigma_z(N,config,spin_back)-h))
    end
    #PBC
    # if spin == N
    #     spin_forward =1
    #     end
    # if spin == 1
    #     spin_back = N
    # end
    # p_top = exp(-beta*sigma_z(N,config,spin)*(sigma_z(N,config,spin_forward)+sigma_z(N,config,spin_back))-beta*h*sigma_z(N,config,spin))
    # p_bot = exp(-beta*(sigma_z(N,config,spin_forward)+sigma_z(N,config,spin_back)-h))+exp(beta*(sigma_z(N,config,spin_forward)+sigma_z(N,config,spin_back)-h))
    return p_top/p_bot
end

function mixing(N, beta,h)
    dim = (2)^N
    M = zeros(dim,dim)
    for ket in (0:dim-1)
        ket_binary = bitstring(ket)
        p_sum = 0
        for SpinIndex in (0:N-1)
            bit = Int(2)^(SpinIndex)
            bra = ket ⊻ bit
            # Sigma_x term
            M[bra+1,ket+1] = (1/N)*accept_prob(N, beta, bra,SpinIndex+1,h)
            p_sum += (1/N)*accept_prob(N, beta, ket,SpinIndex+1,h)
        end
        # Diagonal term
        M[ket+1,ket+1] += 1-p_sum
    end
    return M |> sparse
end


# beta = 6
# h=0
# N_values = (4:15)
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
# save_object("Data/GlaubLoc/GlaubLocOBCgapBeta6", gap_all)

N = 12
h = 0
temp = 10 .^ (range(-2.5,stop=2.5,length=50))
beta_values = 1 ./ temp
gap_all = zeros(length(beta_values))
e1 = zeros(length(beta_values))
e2 = zeros(length(beta_values))
e3 = zeros(length(beta_values))
for j in (1:length(beta_values))
    beta = beta_values[j]
    println(" Working on beta = ",beta)
    M = mixing(N,beta,h)
    e,v  = eigs(M, nev = 3, which=:LM)
    e1[j] = abs(e[1])
    e2[j] = abs(e[2])
    e3[j] = abs(e[3])
    # print("\n largest ", e[1])
    # print(" second ", e[2])
    # print(" third ",e[3])
    gap_all[j] = abs(1-e[2])
end
save_object("Data/GlaubLoc/GlaubLocOBCgapN12", gap_all)
# save_object("Data/GlaubLoc/GlaubLocF0.1OBCgapN10e1", e1)
# save_object("Data/GlaubLoc/GlaubLocF0.1OBCgapN10e2", e2)
# save_object("Data/GlaubLoc/GlaubLocF0.1OBCgapN10e3", e3)