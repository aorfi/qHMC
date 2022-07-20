using LinearAlgebra
using SparseArrays
using Arpack
using JLD2
using PyPlot
using OrderedCollections
include("qHMC.jl")


N = 8
alpha = 1.5
couplings = ones(N)
couplings[end] = 0 
h=0
frac = 4
beta = 6
num_values = 500
eta_values = range(0,10, length=num_values)
# prob = zeros(num_values)
U1 = zeros(2^N,num_values)im
for i in (1:num_values)
    println("Working on eta ", i)
    H = ham(alpha, eta_values[i] , N,couplings, 0)
    U = exp(-1im*H)
    U1[:,i]= U[:,1]
    # prob[i] = average_prob(alpha, eta_values[i], N,couplings, 0, frac)
end
name = "Data/GridSearch/alphaEtaParam/AlphaZero/Alpha1.5U1"*string(num_values)*"N"*string(N)*"beta"*string(beta)
save_object(name, U1)

name = "Data/GridSearch/alphaEtaParam/Alpha1.5"*string(num_values)*"N"*string(N)*"beta"*string(beta)
gap = load_object(name)
name = "Data/GridSearch/alphaEtaParam/AlphaZero/Alpha1.5U1"*string(num_values)*"N"*string(N)*"beta"*string(beta)

U1 = load_object(name)

for i in (2:2^N)
    plt.plot(eta_values,abs.(U1[i,:]).^2,linestyle= "dashed")
end

plt.plot(eta_values, gap, label = "Gap", color = "black")
# plt.plot(eta_values,abs.(U1[1,:]).^2, label ="1", linestyle= "dashed")
# plt.plot(eta_values,1/4*cos.(3*eta_values)+3/4*cos.(eta_values), linestyle= "dashed")
# plt.plot(eta_values,abs.(U1[2,:]).^2, label ="Prob 1 -> 2/3/5", linestyle= "dashed")
# plt.plot(eta_values,-1/4*sin.(3*eta_values)-1/4*sin.(eta_values), linestyle= "dashed")
# plt.plot(eta_values,abs.(U1[4,:]).^2, label ="Prob 1 -> 4/6/7", linestyle= "dashed")
# plt.plot(eta_values,1/4*cos.(3*eta_values)-1/4*cos.(eta_values), linestyle= "dashed")
# plt.plot(eta_values,abs.(U1[8,:]).^2, label ="Prob 1 -> 8", linestyle= "dashed")
# plt.plot(eta_values,-1/4*sin.(3*eta_values)+3/4*sin.(eta_values), linestyle= "dashed")


plt.title(L"$N=$ "*string(N)*L" $\alpha=$ "*string(alpha))
plt.xlabel(L"$\eta$")
plt.legend()
plt.grid()
# plt.ylim(0,0.1)
# name = "Figures/Ising/GapInvestigation/Alpha8AllStatesN"*string(N)*".png"
# plt.savefig(name)
plt.show()




# alpha = 0
# h=0
# beta = 6
# num_values = 500
# eta_values = range(0,10, length=num_values)
# N_values = (4:8)
# for N in N_values
#     name = "Data/GridSearch/alphaEtaParam/AlphaZero/"*string(num_values)*"N"*string(N)*"beta"*string(beta)
#     gap = load_object(name)
#     plt.plot(eta_values, gap, label = string(N))
# end
# plt.title(L"Gap Scaling $\alpha=$ "*string(alpha))
# plt.xlabel(L"$\eta$")
# plt.ylabel("gap")
# plt.legend()
# plt.grid()
# plt.show()