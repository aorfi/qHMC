using LaTeXStrings
using Arpack
using PyPlot
using JLD2
using LsqFit

# plot gap vs temp for fixed N
N = 10
temp = 10 .^ (range(-2.5,stop=2.5,length=50))
gapMH = load_object("Data/MH/MHgapOBCN10")

gapG = load_object("Data/MH/MHgapOBCN10third")
gapGl = load_object("Data/GlaubLoc/GlaubLocOBCgapN10")
# gap = load_object("Data/Ggamma0.75t8gapN10")
# gap5 = load_object("Data/Ggamma0.5t8gapN10")
# plt.scatter(temp, gapG, label = "Glauber Uniform")
plt.scatter(temp, gapG, label = L"MH Uniform $1-e[3]$")
plt.scatter(temp, gapMH, label = L"MH Uniform $1-e[2]$")
# plt.scatter(temp, gap, label = L"qHMC $\gamma = 0.75$ $t=8$")
# plt.scatter(temp, gap5, label = L"qHMC $\gamma = 0.5$ $t=8$")
plt.scatter(temp, gapGl, label = "Glauber Local")
plt.title(string(L"Spectral Gap  $N= $ ", 10))
plt.ylabel(L"$\delta$")
plt.xlabel(L"$T$")
plt.yscale("log")
plt.xscale("log")
plt.grid("both","major")
plt.legend(loc=4)
# plt.savefig("Figures/ParamGapN10.png")
plt.show()