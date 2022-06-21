using LinearAlgebra
using SparseArrays
using Arpack
using JLD2

function ising_ham(N)
    dim = (2)^N
    H = zeros(dim,dim)
    for ket in (0:dim-1)
        ket_binary = bitstring(ket)
        Diagonal = Int64(0)
        for SpinIndex in (0:N-1)
            if SpinIndex == N-1
                Si = 2*((ket>>SpinIndex)&1)-1
                Si_next = 2*((ket>>(0))&1)-1
                Bond = Si*Si_next
                H[ket+1,ket+1] += -1 * Bond
                break
            end
            Si = 2*((ket>>SpinIndex)&1)-1
            Si_next = 2*((ket>>(SpinIndex+1))&1)-1
            Bond = Si*Si_next
            Diagonal += Bond
        end
        H[ket+1,ket+1] += -1 * Diagonal
    end
    return H #|> sparse
end

function mixing_ham(N)
    return 1/sqrt(2^N)*ones(2^N,2^N)
end


function ham(gamma,N)
    H = (1-gamma)*ising_ham(N) + gamma*mixing_ham(N)
    return H
end

N = 2
H = ham(0.5,N)
e,v = eigen(H)
display(H)
display(e)
display(v)