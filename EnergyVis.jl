using LinearAlgebra
using SparseArrays
using Arpack
using JLD2
using PyPlot
using OrderedCollections
include("qHMC.jl")

# OBC
function ising_energy(N, config)
    # config is in [0,2^N] and spin in [1,N]
    eng = 0
    for SpinIndex in (0:N-1)
        if SpinIndex == N-1
            break
        end
        Si = 2*((config>>SpinIndex)&1)-1
        Si_next = 2*((config>>(SpinIndex+1))&1)-1
        eng += -Si*Si_next
    end
    return eng
end

function sort_mixing(N,Q)
    configs = (0:2^N-1)
    eng = [ising_energy(N, i) for i in configs]
    permvec = sortperm(eng)
    sortedM = zeros(2^N,2^N)
    for j in 1:2^N
        sortedM[:,j] = Q[:,j][permvec]
    end
    for j in 1:2^N
        sortedM[j,:] = sortedM[j,:][permvec]
    end
    return sortedM
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

function average_prob(alpha, eta, N,couplings, h, frac)
    H = ham(alpha, eta, N,couplings, h)
    U = exp(-1im*H)
    prob = (abs.(U)).^2
    B = sort_mixing(N,prob)
    for i in (1:length(B[:,1]))
        B[i,i] = 0
    end
    cutoff = trunc(Int,length(B[:,1])/frac)
    let
        prob_cutoff = 0 
        counter = 0
        for i in (1:length(B[:,1]))
            for j in range(1,cutoff)
                if (i+j)<= length(B[:,1])
                    prob_cutoff += B[i,i+j]
                    counter += 1
                end
            end
        end
    return cutoff/counter
    end
end


# N = 8
# eta = 6
# alpha = 3
# couplings = ones(N)
# couplings[end] = 0 
# h=0
# H = ham(alpha, eta, N,couplings, 0)
# U = exp(-1im*H)
# prob = (abs.(U)).^2
# B = sort_mixing(N,prob)
# for i in (1:length(B[:,1]))
#     B[i,i] = 0
# end

# plt.title(L"Proposal Probability $N = $"*string(N)*L" $\eta =$ "*string(eta)*L" $\alpha = $ "*string(alpha))
# plt.imshow(B, origin="lower",cmap="viridis")
# bar = plt.colorbar()
# plt.xlabel(L"configurations by increasing energy $\rightarrow$")
# plt.ylabel(L"configurations by increasing energy $\rightarrow$")

# bar.set_label("Probability")
# plt.xticks([])
# plt.yticks([])
# name = "Figures/Ising/ProposalProb/N"*string(N)*"eta"*string(eta)*"alpha"*string(alpha)*".png"
# plt.savefig(name)
# # plt.clf()
# plt.show()



N = 8
num_values = 100
beta = 6
eta_values = range(0,30, length=num_values)
alpha_values = range(0,30, length=num_values)
couplings = ones(N)
couplings[end] = 0 
prob_all = zeros(num_values,num_values)

for alpha_i in (1:num_values)
    for eta_i in (1:num_values)
        print("  alpha: ", alpha_i)
        println("  eta: ", eta_i)
        prob_all[alpha_i,eta_i] = average_prob(alpha_values[alpha_i], eta_values[eta_i], N,couplings, 0, 10)
    end
end

name = "Data/GridSearch/alphaEtaParam/PropProp"*string(num_values)*"N"*string(N)*"beta"*string(beta)
save_object(name, prob_all)

