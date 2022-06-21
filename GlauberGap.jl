using LinearAlgebra
using SparseArrays
using Arpack
using JLD2



function sigma_z(N, config, spin)
    # config is in [0,2^N] and spin in [1,N]
    spin_index = spin -1
    return 2*((config>>spin_index)&1)-1
end

function accept_prob(N, beta, config,spin)
    spin_forward = spin + 1
    spin_back = spin - 1
    if spin == N
        spin_forward = 1
    end
    if spin == 1
        spin_back = N
    end
    p_top = exp(-beta*sigma_z(N,config,spin)*(sigma_z(N,config,spin_forward)+sigma_z(N,config,spin_back)))
    p_bot = exp(-beta*(sigma_z(N,config,spin_forward)+sigma_z(N,config,spin_back)))+exp(beta*(sigma_z(N,config,spin_forward)+sigma_z(N,config,spin_back)))
    return p_top/(N*p_bot)
end

function mixing(N, beta)
    dim = (2)^N
    M = zeros(dim,dim)
    for ket in (0:dim-1)
        ket_binary = bitstring(ket)
        Diagonal = Int64(0)
        p_sum = 0
        for SpinIndex in (0:N-1)
            bit = Int(2)^(SpinIndex)
            bra = ket ⊻ bit
            M[bra+1,ket+1] += accept_prob(N, beta, bra,SpinIndex+1)
            p_sum += accept_prob(N, beta, ket,SpinIndex+1)
        end
        M[ket+1,ket+1] += 1-p_sum
    end
    return M |> sparse
end


beta = 10
N_values = (2:15)
N_max = last(N_values)
gap_all = zeros(length(N_values))
for j in (1:length(N_values))
    N = N_values[j]
    println(" Working on N = ",N)
    M = mixing(N,beta)
    e,v  = eigs(M, nev = 2, which=:LR)
    gap_all[j] = abs(e[1]-e[2])
end
save_object("Data/gapBeta10", gap_all)