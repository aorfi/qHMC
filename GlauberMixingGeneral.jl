using LinearAlgebra
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
    for col in (0:dim-1)
        E_sp = ising_energy(N,couplings,h,col)
        diag = 0 
        for row in (0:dim-1)
            if row != col
                E_s = ising_energy(N,couplings,h,row)
                m = (1/2^N)*(1/(1+exp(-beta*(E_s-E_sp))))
                M[col+1,row+1] = m
                diag += m
            end
        M[col+1,col+1] = 1-diag
        end
    end
    return M
end

beta = 6
h=0
N = 3
couplings = ones(N)
couplings[end] = 0
M = mixing(N,couplings,h,beta)

# e,v  = eigs(M, nev = 3, which=:LM)

display(M)
display(sum(M[1,:]))
# display(e)

# beta = 6
# N_values = (5:15)
# h=0
# gap_all = zeros(length(N_values))
# for j in (1:length(N_values))
#     N = N_values[j]
#     println(" Working on N = ",N)
#     M = mixing(N,beta,h)
#     # e,v  = eigen(M)
#     # gap_all[j] = abs(1-e[end-1])
#     e,v  = eigs(M, nev = 3, which=:LM)
#     gap_all[j] = 1-abs(e[2])
#     print("Gap: ",1-abs(e[2]))
# end
# save_object("Data/Glaub/CompareOBCGlaubGengapBeta6", gap_all)

# N = 12
# h=0
# temp = 10 .^ (range(-2.5,stop=2.5,length=50))
# beta_values = 1 ./ temp
# gap_all = zeros(length(beta_values))
# for j in (1:length(beta_values))
#     beta = beta_values[j]
#     println(" Working on beta = ",beta)
#     M = mixing(N,beta,h)
#     e,v  = eigs(M, nev = 3, which=:LM)
#     gap_all[j] = abs(1-e[2])
#     # gap_all[j] = abs(1-e[3])
# end
# save_object("Data/Glaub/GlaubGapOBCN12", gap_all)

# N = 10
# h=0.1
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
#     print("\n largest ", e[1])
#     print(" second ", e[2])
#     print(" third ",e[3])
#     gap_all[j] = abs(1-e[2])
# end
# save_object("Data/Glaub/GlaubF0.1OBCgapN10e1", e1)
# save_object("Data/Glaub/GlaubF0.1OBCgapN10e1", e1)
# save_object("Data/Glaub/GlaubF0.1OBCgapN10e2", e2)
# save_object("Data/Glaub/GlaubF0.1OBCgapN10e3", e3)