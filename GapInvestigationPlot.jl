using LinearAlgebra
using SparseArrays
using Arpack
using JLD2
using PyPlot
using PlotUtils
using OrderedCollections
include("qHMC.jl")
include("GlauberMixingGeneral.jl")

N = 3
alpha = 0
couplings = ones(N)
couplings[end] = 0 
h=0
frac = 4
beta = 6
num_values = 500
eta_values = range(0,10, length=num_values)
U1 = zeros(2^N,num_values)im
for i in (1:num_values)
    println("Working on eta ", i)
    H = ham(alpha, eta_values[i] , N,couplings, 0)
    U = exp(-1im*H)
    U1[:,i]= U[:,1]
end
# name = "Data/GridSearch/alphaEtaParam/AlphaZero/Alpha1.5U1"*string(num_values)*"N"*string(N)*"beta"*string(beta)
# save_object(name, U1)

name = "Data/GridSearch/alphaEtaParam/AlphaZero/"*string(num_values)*"N"*string(N)*"beta"*string(beta)*"e"
eigens = load_object(name)
M = mixing_glab(N,couplings,h,beta)
e,v  = eigen(M)
gapG = 1-abs(e[end-1])
# name = "Data/GridSearch/alphaEtaParam/AlphaZero/Alpha1.5U1"*string(num_values)*"N"*string(N)*"beta"*string(beta)
# U1 = load_object(name)

# for i in (2:2^N)
#     plt.plot(eta_values,abs.(U1[i,:]).^2,linestyle= "dashed")
# end
plt.plot(eta_values, gapG*ones(length(eta_values)), linestyle = "dashdot", color = "red")
plt.plot(eta_values, eigens[end-1,:], label = "1-e2", color = "black")
# plt.plot(eta_values, eigens[end-2,:], label = "1-e3")
# plt.plot(eta_values, eigens[end-3,:], label = "1-e4")
# plt.plot(eta_values, eigens[end-4,:], label = "1-e4")
# plt.plot(eta_values, eigens[end-5,:], label = "1-e4")
# plt.plot(eta_values, eigens[end-6,:], label = "1-e4")
# plt.plot(eta_values,abs.(U1[1,:]).^2, label ="1", linestyle= "dashed")
# plt.plot(eta_values,1/4*cos.(3*eta_values)+3/4*cos.(eta_values), linestyle= "dashed")
plt.plot(eta_values,abs.(U1[2,:]).^2, label ="Prob 1 -> 2/3/5", linestyle= "dashed")
# plt.plot(eta_values,-1/4*sin.(3*eta_values)-1/4*sin.(eta_values), linestyle= "dashed")
plt.plot(eta_values,abs.(U1[4,:]).^2, label ="Prob 1 -> 4/6/7", linestyle= "dashed")
# plt.plot(eta_values,1/4*cos.(3*eta_values)-1/4*cos.(eta_values), linestyle= "dashed")
plt.plot(eta_values,abs.(U1[8,:]).^2, label ="Prob 1 -> 8", linestyle= "dashed")
# plt.plot(eta_values,-1/4*sin.(3*eta_values)+3/4*sin.(eta_values), linestyle= "dashed")


plt.title(L"$N=$ "*string(N)*L" $\alpha=$ "*string(alpha))
plt.xlabel(L"$\eta$")
plt.legend()
plt.grid()
# plt.yscale("log")
# plt.ylim(0,0.1)
# name = "Figures/Ising/GapInvestigation/Alpha8AllStatesN"*string(N)*".png"
# plt.savefig(name)
plt.show()




# alpha = 0
# h=0
# beta = 6
# num_values = 500
# eta_values = range(0,10, length=num_values)
# gapG = load_object("Data/Glaub/OBCBeta6N4-13")
# N_values_big = (4:13)
# N_values = (4:8)

# for i in (1:length(N_values))
#     N = N_values[i]
#     name = "Data/GridSearch/alphaEtaParam/AlphaZero/"*string(num_values)*"N"*string(N)*"beta"*string(beta)
#     gap = load_object(name)
#     c = (i%2,(i+1)%2,i%3/2,1)
#     plt.plot(eta_values, gap, label = string(N), color = c)
#     plt.plot(eta_values,gapG[i]*ones(num_values), color = c, linestyle="dashed")
# end
# plt.title(L"Gap Scaling $\alpha=$ "*string(alpha))
# plt.xlabel(L"$\eta$")
# plt.ylabel("gap")
# plt.legend()
# plt.grid()
# plt.show()