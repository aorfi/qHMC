using Arpack
using JLD2
include("qHMC.jl")



# Get gaps for different beta
# N = 5
# num_values = 100
# gamma_values = range(0,1, length=num_values)
# t_values = range(0,25, length=num_values)
# beta_values = exp10.(-2:2)
# for beta in beta_values
#     gap_all = zeros(num_values,num_values)
#     for g_i in (1:num_values)
#         # H = ham(gamma_values[g_i],N)
#         # e,v = eigen(H)
#         for t_i in (1:num_values)
#             print("\nworking on gamma: ", g_i)
#             print("  t: ", t_i)
#             # M = mixing_fixed(N,beta,t_values[t_i],e,v)
#             M = mixing_matrix(N,beta,gamma_values[g_i],t_values[t_i])
#             eM,v  = eigen(M)
#             gap_all[1,t_i] = abs(1-eM[end-1])
#         end
#     end
#     name = "Data/GridSearch/SigmaX/GlaubGrid100N"*string(N)*"beta"*string(beta)
#     save_object(name, gap_all)
# end

# Get gaps for different N
beta = 6
num_values = 100
alpha_values = range(0,30, length=num_values)
eta_values = range(0,30, length=num_values)
N_values = (5:5)
for N in N_values
    gap_all = zeros(num_values,num_values)
    for alpha_i in (1:num_values)
        for eta_i in (1:num_values)
            print("\nworking on N: ", N)
            print("\nworking on alpha: ", alpha_i)
            print("  eta: ", eta_i)
            M = mixing_matrix(N,beta,alpha_values[alpha_i], eta_values[eta_i])
            # e,v  = eigen(M)
            # gap_all[alpha_i,eta_i]= abs(1-e[end-1])
            e,v  = eigs(M, nev = 2, which=:LM)
            gap_all[alpha_i,eta_i] = 1-abs(e[2])
        end
    end
    name = "Data/GridSearch/alphaEtaParam/PBC"*string(num_values)*"N"*string(N)*"beta"*string(beta)
    save_object(name, gap_all)
end

