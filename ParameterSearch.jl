using LaTeXStrings
using Arpack
using PyPlot
using JLD2
using LsqFit
include("1DIsing.jl")


N = 10
beta = 100
num_values = 100
gamma_values = range(0,1, length=num_values)
t_values = range(0,25, length=num_values)

# Get gaps
# gap_all = zeros(num_values,num_values)
# for g_i in (1:num_values)
#     H = ham(gamma_values[g_i],N)
#     e,v = eigen(H)
#     for t_i in (1:num_values)
#         print("\nworking on gamma: ", g_i)
#         print("  t: ", t_i)
#         M = mixing_fixed(N,beta,t_values[t_i],e,v)
#         eM,vM  = eigs(M, nev = 2, which=:LR)
#         gap_all[g_i,t_i] = abs(eM[1]-eM[2])
#     end
# end
# save_object("Data/gammaTimeGrid50", gap_all)


gap_all= load_object("Data/gammaTimeGrid100")
plt.title(L"qHMC Spectral Gap  $\beta= 100$ $N=10$")
plt.imshow(log.(gap_all),extent = [0, 25, 0 , 1], aspect="auto")
plt.xlabel(L"$t$")
plt.ylabel(L"$\gamma$")
bar = plt.colorbar()
bar.set_label(L"$\delta$")
# plt.savefig("Figures/grid100.png")
plt.show()