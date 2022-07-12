using LaTeXStrings
using Arpack
using PyPlot
using JLD2
using LsqFit

# # plot gap vs temp for fixed N CLASSICAL
N = 10
h = 0
temp = 10 .^ (range(-2.5,stop=2.5,length=50))
gapMH = load_object("Data/SK/Classical/avMHN10")
gapMHl = load_object("Data/SK/Classical/avMHLocN10")
gapG = load_object("Data/SK/Classical/avGlaubN10")
gapGl = load_object("Data/SK/Classical/avGlaubLocN10")
gap = load_object("Data/SK/qHMC/alpha1.5eta2.5N10")

plt.scatter(temp, gapMH, label = "MH Uniform")
plt.scatter(temp, gapMHl, label = "MH Local")
plt.scatter(temp, gapG, label = "Glaubler Uniform")
plt.scatter(temp, gapGl, label = "Glaubler Local")
plt.scatter(temp, gap, label = L"qHMC $\alpha = 1.5$ $\eta = 2.5$")

plt.title(L"MCMC SK Model $N= $ "*string(N))
plt.ylabel(L"Average $\delta$")
plt.xlabel(L"$T$")
plt.yscale("log")
plt.xscale("log")
plt.grid("both","major")
plt.legend(loc=4)
plt.savefig("Figures/SK/avGapN10.png")
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