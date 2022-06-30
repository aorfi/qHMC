using LaTeXStrings
using Arpack
using PyPlot
using JLD2
using LsqFit

# # plot gap vs temp for fixed N CLASSICAL
# N = 12
# temp = 10 .^ (range(-2.5,stop=2.5,length=50))
# gapMH = load_object("Data/MH/MHgapOBCN12")
# # gapMH3 = load_object("Data/MH/MHgapOBCN10third")
# gapMHl = load_object("Data/MHLoc/MHLocGapOBCN12")
# # gapMHl3 = load_object("Data/MHLoc/MHLocGapOBCN10third")
# gapG = load_object("Data/Glaub/GlaubGapOBCN12")
# gapGl = load_object("Data/GlaubLoc/GlaubLocOBCgapN12")
# # gapGl3 = load_object("Data/GlaubLoc/GlaubLocOBCgapN10third")

# plt.scatter(temp, gapMH, label = "MH Uniform")
# plt.scatter(temp, gapMHl, label = "MH Local")
# plt.scatter(temp, gapG, label = "Glaubler Uniform")
# plt.scatter(temp, gapGl, label = "Glaubler Local")

# plt.title(string(L"Classical MCMC with OBC $N= $ ", N))
# plt.ylabel(L"$\delta$")
# plt.xlabel(L"$T$")
# plt.yscale("log")
# plt.xscale("log")
# plt.grid("both","major")
# plt.legend(loc=4)
# # plt.savefig("Figures/ClassicalOBCGapN12.png")
# plt.show()




# plot gap vs temp for fixed N CLASSICAL
N = 10
temp = 10 .^ (range(-2.5,stop=2.5,length=50))
gapG = load_object("Data/Glaub/GlaubGapOBCN10")
gap  = load_object("Data/qHMC/SigmaX/gamma0.5t20gapN10")

plt.scatter(temp, gap, label = L"Glauber qHMC $\gamma = 0.5$ $t=20$ ")
plt.scatter(temp, gapG, label = "Glauber Uniform")


plt.title(string(L"Gap Comparison $N= $ ", N))
plt.ylabel(L"$\delta$")
plt.xlabel(L"$T$")
plt.yscale("log")
plt.xscale("log")
plt.grid("both","major")
plt.legend()
plt.savefig("Figures/qHMCScaling/GapGamma0.5t20N10.png")
plt.show()