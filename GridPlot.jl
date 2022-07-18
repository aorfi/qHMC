using PyPlot
using LaTeXStrings
using JLD2
include("qHMC.jl")

N=3
beta = 6
num_values = 200
max_param = 30
name =  "Data/GridSearch/alphaEtaParam/"*string(num_values)*"N"*string(N)*"beta"*string(beta)
gap_all= load_object(name)
max,cord = findmax(gap_all)
println((cord[2]-1)/num_values*max_param)
println((cord[1]-1)/num_values*max_param)


plt.title(L"qHMC Spectral Gap $\beta= $ "*string(beta)*L" $N = $"*string(N))
# plt.imshow(gap_all,extent = [0, 30, 0 , 30], aspect="auto")
plt.imshow(gap_all, origin="lower",extent = [0, max_param, 0 , max_param])
# plt.clim(0, 0.01) 
bar = plt.colorbar()
plt.scatter((cord[2]-1)/num_values*max_param, (cord[1]-1)/num_values*max_param, color="red", marker=".")
plt.xlabel(L"$\eta$")
plt.ylabel(L"$\alpha$")

bar.set_label(L"$\delta$")
name = "Figures/Ising/GridSearch/alphaEtaParam/"*string(num_values)*"N"*string(N)*"beta"*string(beta)*".png"
plt.savefig(name)
# plt.clf()
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


