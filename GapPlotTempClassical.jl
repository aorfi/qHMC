using LaTeXStrings
using Arpack
using PyPlot
using JLD2
using LsqFit

# # plot gap vs temp for fixed N CLASSICAL
N = 10
h = 1.5
temp = 10 .^ (range(-2.5,stop=2.5,length=50))
gapMH = load_object("Data/MH/PBCN10h"*string(h))
gapMHl = load_object("Data/MHLoc/PBCN10h"*string(h))
gapG = load_object("Data/Glaub/PBCN10h"*string(h))
gapGl = load_object("Data/GlaubLoc/PBCN10h"*string(h))

plt.scatter(temp, gapMH, label = "MH Uniform")
plt.scatter(temp, gapMHl, label = "MH Local")
plt.scatter(temp, gapG, label = "Glaubler Uniform")
plt.scatter(temp, gapGl, label = "Glaubler Local")

plt.title(L"Classical MCMC $N= $ "*string(N)*" h = "*string(h))
plt.ylabel(L"$\delta$")
plt.xlabel(L"$T$")
plt.yscale("log")
plt.xscale("log")
plt.grid("both","major")
plt.legend(loc=4)
plt.savefig("Figures/Ising/Classical/PBCGapN10h"*string(h)*".png")
plt.show()




# N = 10
# temp = 10 .^ (range(-2.5,stop=2.5,length=50))
# gapG = load_object("Data/Glaub/Old/GlaubGapOBCN10")
# gap  = load_object("Data/qHMC/SigmaX/gamma0.5t20gapN10")

# plt.scatter(temp, gap, label = L"Glauber qHMC $\gamma = 0.5$ $t=20$ ")
# plt.scatter(temp, gapG, label = "Glauber Uniform")


# plt.title(string(L"Gap Comparison $N= $ ", N))
# plt.ylabel(L"$\delta$")
# plt.xlabel(L"$T$")
# plt.yscale("log")
# plt.xscale("log")
# plt.grid("both","major")
# plt.legend()
# # plt.savefig("Figures/qHMCScaling/GapGamma0.5t20N10.png")
# plt.show()