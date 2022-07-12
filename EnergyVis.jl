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

N = 8
eta = 20.7
alpha = 6.45
couplings = ones(N)
couplings[end] = 0 
H = ham(alpha, eta, N,couplings, 0)
U = exp(-1im*H)
prob = (abs.(U)).^2
B = sort_mixing(N,prob)
for i in (1:length(B[:,1]))
    B[i,i] = 0
end
plt.title(L"Proposal Probability $N = $"*string(N)*L" $\eta =$ "*string(eta)*L" $\alpha = $ "*string(alpha))
plt.imshow(B, origin="lower",cmap="viridis")
bar = plt.colorbar()
plt.xlabel(L"configurations by increasing energy $\rightarrow$")
plt.ylabel(L"configurations by increasing energy $\rightarrow$")

bar.set_label("Probability")
plt.xticks([])
plt.yticks([])
# name = "Figures/Ising/ProposalProb/N"*string(N)*"eta"*string(eta)*"alpha"*string(alpha)*".png"
# plt.savefig(name)
# plt.clf()
plt.show()


