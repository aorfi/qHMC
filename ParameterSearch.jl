using Arpack
using PyPlot
using JLD2
using LaTeXStrings
include("1DIsingQHMC.jl")


# Get gaps for different beta
# N = 2
# num_values = 1000
# gamma_values = range(0,1, length=num_values)
# t_values = range(0,25, length=num_values)
# beta_values = exp10.(-2:2)
# for beta in beta_values
#     gap_all = zeros(num_values,num_values)
#     for g_i in (1:num_values)
#         H = ham(gamma_values[g_i],N)
#         e,v = eigen(H)
#         for t_i in (1:num_values)
#             print("\nworking on gamma: ", g_i)
#             print("  t: ", t_i)
#             # M = mixing_fixed(N,beta,t_values[t_i],e,v)
#             M = mixing_fixed_glaub(N,beta,t_values[t_i],e,v)
#             eM,vM  = eigs(M, nev = 2, which=:LR)
#             gap_all[g_i,t_i] = abs(eM[1]-eM[2])
#         end
#     end
#     name = "Data/GridSearch/GlaubGrid1000N"*string(N)*"beta"*string(beta)
#     save_object(name, gap_all)
# end

# Get gaps for different N
# beta = 100
# num_values = 20
# gamma_values = range(0,1, length=num_values)
# t_values = range(0,25, length=num_values)
# N_values = (10:10)
# for N in N_values
#     gap_all = zeros(num_values,num_values)
#     for g_i in (1:num_values)
#         H = ham(gamma_values[g_i],N)
#         e,v = eigen(H)
#         for t_i in (1:num_values)
#             print("\nworking on N: ", N)
#             print("\nworking on gamma: ", g_i)
#             print("  t: ", t_i)
#             # M = mixing_fixed(N,beta,t_values[t_i],e,v)
#             M = mixing_fixed_glaub(N,beta,t_values[t_i],e,v)
#             eM,vM  = eigs(M, nev = 2, which=:LR)
#             gap_all[g_i,t_i] = abs(1-eM[2])
#         end
#     end
#     name = "Data/GridSearch/OBC/OBCGlaubGrid20N"*string(N)*"beta"*string(beta)
#     save_object(name, gap_all)
# end

# N=10
# beta = 100
# name = "Data/GridSearch/OBC/OBCGlaubGrid20N"*string(N)*"beta"*string(beta)
# gap_all= load_object(name)
# plt.title(L"qHMC Spectral Gap  OBC $\beta= $ "*string(beta)*L" $N = $"*string(N))
# plt.imshow(gap_all,extent = [0, 25, 0 , 1], aspect="auto")
# plt.xlabel(L"$t$")
# plt.ylabel(L"$\gamma$")
# bar = plt.colorbar()
# bar.set_label(L"$\delta$")
# name = "Figures/GridSearch/OBC/OBCgridG20N"*string(N)*"beta"*string(beta)*".png"
# plt.savefig(name)
# plt.show()