using LinearAlgebra
using SparseArrays
using Arpack
using JLD2

function ising_energy(N, config)
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
            # Si = 2*((config>>SpinIndex)&1)-1
            # eng += h*Si
            break
        end
        Si = 2*((config>>SpinIndex)&1)-1
        # eng += h*Si
        Si_next = 2*((config>>(SpinIndex+1))&1)-1
        eng += -Si*Si_next
    end
    return eng
end

function ising_ham(N)
    dim = (2)^N
    H = zeros(dim,dim)
    for ket in (0:dim-1)
        ket_binary = bitstring(ket)
        Diagonal = Int64(0)
        for SpinIndex in (0:N-1)
            # PBC
            # if SpinIndex == N-1
            #     Si = 2*((ket>>SpinIndex)&1)-1
            #     Diagonal += h*Si
            #     Si_next = 2*((ket>>(1))&1)-1
            #     H[ket+1,ket+1] += -1 * Si*Si_next
            #     break
            # end
            # OBC
            if SpinIndex == N-1
                # Si = 2*((ket>>SpinIndex)&1)-1
                # Diagonal += h*Si
                break
            end
            Si = 2*((ket>>SpinIndex)&1)-1
            # Diagonal += h*Si
            Si_next = 2*((ket>>(SpinIndex+1))&1)-1
            Diagonal += -1*Si*Si_next
        end
        H[ket+1,ket+1] += Diagonal
    end
    return H #|> sparse
end

function mixing_ham_dense(N)
    return 1/sqrt(2^N)*ones(2^N,2^N)
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

function ham(gamma,N)
    H = (1-gamma)*ising_ham(N) + gamma*mixing_ham(N)
    return H
end

function mixing_matrix(N,beta,gamma,t)
    H = ham(gamma,N)
    U = exp(-1im*t*H)
    prob = (abs.(U)).^2
    dim = (2)^N
    M = zeros(dim,dim)
    for col in (0:dim-1)
        E_sp = ising_energy(N,col)
        diag = 0 
        for row in (0:dim-1)
            if row != col
                E_s = ising_energy(N,row)
                m = prob[col+1,row+1]*min(1,exp(-beta*(E_sp-E_s)))
                M[col+1,row+1] = m
                diag += m
            end
        M[col+1,col+1] = 1-diag
        end
    end
    return M 
end

function mixing_fixed(N,beta,t,e,v)
    U = v*diagm(exp.(-1im*t*e))*inv(v)
    prob = (abs.(U)).^2
    dim = (2)^N
    M = zeros(dim,dim)
    for col in (0:dim-1)
        E_sp = ising_energy(N,col)
        diag = 0 
        for row in (0:dim-1)
            if row != col
                E_s = ising_energy(N,row)
                m = prob[col+1,row+1]*min(1,exp(-beta*(E_sp-E_s)))
                M[col+1,row+1] = m
                diag += m
            end
        M[col+1,col+1] = 1-diag
        end
    end
    return M 
end

function mixing_fixed_glaub(N,beta,t,e,v)
    U = v*diagm(exp.(-1im*t*e))*inv(v)
    prob = (abs.(U)).^2
    dim = (2)^N
    M = zeros(dim,dim)
    for col in (0:dim-1)
        E_sp = ising_energy(N,col)
        diag = 0 
        for row in (0:dim-1)
            if row != col
                E_s = ising_energy(N,row)
                m = prob[col+1,row+1]*(1/(1+exp(-beta*(E_s-E_sp))))
                M[col+1,row+1] = m
                diag += m
            end
        M[col+1,col+1] = 1-diag
        end
    end
    return M 
end

# N = 10
# gamma = 0.5
# t = 20
# temp = 10 .^ (range(-2.5,stop=2.5,length=50))
# beta_values = 1 ./ temp
# gap_all = zeros(length(beta_values))
# for j in (1:length(beta_values))
#     beta = beta_values[j]
#     println(" Working on beta = ",beta)
#     M = mixing_matrix(N,beta,gamma,t)
#     e,v  = eigs(M, nev = 2, which=:LR)
#     gap_all[j] = abs(1-e[2])
# end
# save_object("Data/qHMC/SigmaX/gamma0.5t20gapN10", gap_all)


# beta = 300
# gamma = 0.5
# t = 20
# N_values = (2:12)
# gap_all = zeros(length(N_values))
# for j in (1:length(N_values))
#     N = N_values[j]
#     println(" Working on N = ",N)
#     M = mixing_matrix(N,beta,gamma,t)
#     e,v  = eigs(M, nev = 2, which=:LR)
#     gap_all[j] = abs(1-e[2])
#     print("Gap: ",abs(1-e[2]))
# end
# save_object("Data/qHMC/SigmaX/gamma0.5t20gapBeta300", gap_all)

