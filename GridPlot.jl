using PyPlot
using LaTeXStrings
using JLD2
include("qHMC.jl")

# N_values = (5:10)
# beta = 6
# max_param = 30
# frac = zeros(length(N_values))
# for i in (1:length(N_values))
#     N = N_values[i]
#     num_values = 200
#     if N == 10
#         num_values = 100
#     end
#     name =  "Data/GridSearch/alphaEtaParam/"*string(num_values)*"N"*string(N)*"beta"*string(beta)
#     gap_all= load_object(name)
#     gapG = load_object("Data/Glaub/OBCBeta6")[N-4]
    
#     cutoff = zeros(num_values, num_values)
#     counter = 0
#     all = 0
#     for i in (1:num_values)
#         for j in (1:num_values)
#             if i<=j
#                 all += 1
#                 if gap_all[i,j]>=gapG
#                     cutoff[i,j] = 1
#                     counter += 1
#                 end
#             end
#         end
#     end

#     frac[i] = counter/all

#     # plt.title(L"qHMC Spectral Gap $\beta= $ "*string(beta)*L" $N = $"*string(N))
#     # plt.imshow(gap_all,extent = [0, 30, 0 , 30], aspect="auto")
#     # plt.imshow(cutoff, origin="lower",extent = [0, max_param, 0 , max_param])
#     # # plt.clim(0, 0.01) 
#     # # bar = plt.colorbar()
#     # # plt.scatter((cord[2]-1)/num_values*max_param, (cord[1]-1)/num_values*max_param, color="red", marker=".")
#     # plt.xlabel(L"$\eta$")
#     # plt.ylabel(L"$\alpha$")
#     # bar.set_label(L"$\delta$")
#     # # name = "Figures/Ising/GridSearch/alphaEtaParam/Cutoff/"*string(num_values)*"N"*string(N)*"beta"*string(beta)*".png"
#     # # plt.savefig(name)
#     # # plt.clf()
#     # plt.show()
    
    
# end
# println(frac)
# save_object("Data/qHMC/alphaEtaParam/fracScalingBeta6Triangle",frac)


N_values = (5:10)
frac = load_object("Data/qHMC/alphaEtaParam/fracScalingBeta6")
fracT = load_object("Data/qHMC/alphaEtaParam/fracScalingBeta6Triangle")
plt.scatter(N_values, frac, label = L"$\alpha,\eta\in[0,30]$")
plt.scatter(N_values, fracT, label = L"$\alpha\in[0,30], \eta<\alpha$")
plt.title(L"Fraction of $\alpha,\eta >$ Uniform Gap")
plt.grid()
plt.xlabel("N")
plt.ylabel("Fraction")
plt.legend()
name = "Figures/Ising/GridSearch/alphaEtaParam/Cutoff/Scaling.png"
plt.savefig(name)
plt.show()



# N=8
# beta = 6
# num_values = 1000
# # name =  "Data/GridSearch/alphaEtaParam/SmallAlpha"*string(num_values)*"N"*string(N)*"beta"*string(beta)
# name =  "Data/GridSearch/alphaEtaParam/AlphaZero"*string(num_values)*"N"*string(N)*"beta"*string(beta)
# gap_all= load_object(name)
# maxval = maximum(gap_all)
# positions = [i for (i, x) in enumerate(gap_all) if abs(x-maxval)<0.00002]
# pos = (positions .- 1)*30/num_values
# println(pos)
# # avGap = 0
# # for i in range(1, 9)
# #     global avGap += 1/(9)*(pos[2*i]-pos[2*i-1])
# #     println(pos[2*i]-pos[2*i-1])
# # end
# # print("Average: ", avGap)
# eta_values = range(0,30, length=num_values)
# plt.plot(eta_values,gap_all)
# plt.plot(eta_values, 1 ./ exp.(-1im*eta_values))
# # plt.scatter(pos, gap_all[positions], color= "red", marker=".")
# plt.title(L"qHMC Spectral Gap at $\alpha=0$ for $\beta= $ "*string(beta)*L" $N = $"*string(N))
# plt.grid()
# plt.xlabel(L"$\eta$")
# plt.ylabel(L"$\delta$")
# # name = "Figures/Ising/GridSearch/alphaEtaParam/AlphaZeroN"*string(N)*"beta"*string(beta)*".png"
# # plt.savefig(name)
# plt.show()


