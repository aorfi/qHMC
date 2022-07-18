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
    sort_prob = sort_mixing(N,prob)
    for i in (1:2^N)
        sort_prob[i,i] = 0
    end
    cutoff = 2^N-trunc(Int,2^N/frac)
    # zero everything not below the cutoff
    for i in (1:2^N)
        for j in range(2^N-cutoff,2^N)
            if (i+j)<= 2^N
                sort_prob[i,i+j] = 0
                sort_prob[i+j,i] = 0
            end
        end
    end
    return sum(sort_prob)/2^N
end


N = 8
alpha = 0
couplings = ones(N)
couplings[end] = 0 
h=0
frac = 4
beta = 6
# H = ham(alpha, 1, N,couplings, 0)
# e,v = eigen(H)
# display(e)
# display(inv(v))
num_values = 500
eta_values = range(0,10, length=num_values)
# prob = zeros(num_values)
# U1 = zeros(2^N,num_values)im
# for i in (1:num_values)
#     println("Working on eta ", i)
#     H = ham(alpha, eta_values[i] , N,couplings, 0)
#     U = exp(-1im*H)
#     U1[:,i]= U[:,1]
#     prob[i] = average_prob(alpha, eta_values[i], N,couplings, 0, frac)
# end
# name = "Data/GridSearch/alphaEtaParam/AlphaZero/U1"*string(num_values)*"N"*string(N)*"beta"*string(beta)
# save_object(name, U1)

name = "Data/GridSearch/alphaEtaParam/AlphaZero/"*string(num_values)*"N"*string(N)*"beta"*string(beta)
gap = load_object(name)
name = "Data/GridSearch/alphaEtaParam/AlphaZero/U1"*string(num_values)*"N"*string(N)*"beta"*string(beta)
U1 = load_object(name)


# plt.plot(eta_values, prob, label = "Proposition prob")

# display(U1[:,5])

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


plt.title(L"$N=8$ $\alpha=0$ ")
plt.xlabel(L"$\eta$")
plt.legend()
plt.grid()
plt.ylim(0,0.1)
# name = "Figures/Ising/GapInvestigation/AlphaZeroAllStatesN"*string(N)*".png"
# plt.savefig(name)
plt.show()


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



# N = 3
# num_values = 200
# beta = 6
# eta_values = range(0,30, length=num_values)
# alpha_values = range(0,30, length=num_values)
# couplings = ones(N)
# couplings[end] = 0 
# frac= 4
# prob_all = zeros(num_values,num_values)

# for alpha_i in (1:num_values)
#     for eta_i in (1:num_values)
#         print("  alpha: ", alpha_i)
#         println("  eta: ", eta_i)
#         prob_all[alpha_i,eta_i] = average_prob(alpha_values[alpha_i], eta_values[eta_i], N,couplings, 0, frac)
#     end
# end

# name = "Data/GridSearch/alphaEtaParam/PropProp"*string(num_values)*"N"*string(N)*"beta"*string(beta)*"frac"*string(frac)
# save_object(name, prob_all)



# N=3
# beta = 6
# num_values = 200
# max_param = 30
# frac= 4
# name = "Data/GridSearch/alphaEtaParam/PropProp"*string(num_values)*"N"*string(N)*"beta"*string(beta)*"frac"*string(frac)
# gap_all= load_object(name)
# max,cord = findmax(gap_all)

# plt.title("Nearest 1/"*string(frac)*L" Configurations $\beta= $ "*string(beta)*L" $N = $"*string(N))
# # plt.imshow(gap_all,extent = [0, 30, 0 , 30], aspect="auto")
# plt.imshow(gap_all, origin="lower",extent = [0, max_param, 0 , max_param])
# # plt.clim(0, 0.01) 
# bar = plt.colorbar()
# plt.scatter((cord[2]-1)/num_values*max_param, (cord[1]-1)/num_values*max_param, color="red", marker=".")
# plt.xlabel(L"$\eta$")
# plt.ylabel(L"$\alpha$")

# bar.set_label("Proposal Probability")
# name = "Figures/Ising/GridSearch/alphaEtaParam/PropProp"*string(num_values)*"N"*string(N)*"beta"*string(beta)*"frac"*string(frac)*".png"
# plt.savefig(name)
# # plt.clf()
# plt.show()


# N = 3
# num_values = 500
# beta = 6
# eta_values = range(0,10, length=num_values)
# alpha=0
# couplings = ones(N)
# couplings[end] = 0 
# frac= 2
# prob_all = zeros(num_values)
# for eta_i in (1:num_values)
#     # print("  alpha: ", alpha_i)
#     println("  eta: ", eta_i)
#     prob_all[eta_i] = average_prob(alpha, eta_values[eta_i], N,couplings, 0, frac)
# end
# name = "Data/GridSearch/alphaEtaParam/AlphaZero/PropProp"*string(num_values)*"N"*string(N)*"beta"*string(beta)*"frac"*string(frac)
# save_object(name, prob_all)

# N = 3
# num_values = 500
# beta = 6
# eta_values = range(0,30, length=num_values)
# alpha=0
# couplings = ones(N)
# couplings[end] = 0 
# frac= 2
# name5 = "Data/GridSearch/alphaEtaParam/AlphaZero/PropProp"*string(num_values)*"N"*string(N)*"beta"*string(beta)*"frac"*string(frac)
# prob5 = load_object(name5)
# # frac= 8
# # name8 = "Data/GridSearch/alphaEtaParam/AlphaZero/PropProp"*string(num_values)*"N"*string(N)*"beta"*string(beta)*"frac"*string(frac)
# # prob8 = load_object(name8)
# name =  "Data/GridSearch/alphaEtaParam/AlphaZero/"*string(num_values)*"N"*string(N)*"beta"*string(beta)
# gap = load_object(name)

# eta_values = range(0,30, length=num_values)

# fig, ax1 = plt.subplots()
# ax2 = ax1.twinx()
# ax1.plot(eta_values,gap,label = "Mixing matrix gap", color = "black")
# ax2.plot(eta_values,prob5,label = "Nearest 1/2", color = "tab:blue")
# # ax2.plot(eta_values,prob8,label = "Nearest 1/8", color = "tab:red")
# plt.title(L"qHMC at $\alpha=0$ for $\beta= $ "*string(beta)*L" $N = $"*string(N))

# plt.xlabel(L"$\eta$")
# ax1.set_ylabel(L"$\delta$", color = "black")
# ax2.set_ylabel("Sum of Proposal Probaility", color = "tab:blue")
# # name = "Figures/Ising/GridSearch/alphaEtaParam/AlphaZeroN"*string(N)*"beta"*string(beta)*".png"
# # plt.savefig(name)
# plt.legend()
# plt.show()