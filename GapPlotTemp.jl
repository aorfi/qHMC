using LaTeXStrings
using Arpack
using PyPlot
using JLD2
using LsqFit

# plot gap vs temp for fixed N
N = 10
temp = 10 .^ (range(-2.5,stop=2.5,length=50))
gapMH = load_object("Data/MH/MHgapOBCN10")
gapMH3 = load_object("Data/MH/MHgapOBCN10third")
gapMHl = load_object("Data/MHLoc/MHLocGapOBCN10")
gapMHl3 = load_object("Data/MHLoc/MHLocGapOBCN10third")
gapGl = load_object("Data/GlaubLoc/GlaubLocOBCgapN10")
gapGl3 = load_object("Data/GlaubLoc/GlaubLocOBCgapN10third")
# gapGle1 = load_object("Data/GlaubLoc/GlaubLocF0.1OBCgapN10e1")
# gapGle2 = load_object("Data/GlaubLoc/GlaubLocF0.1OBCgapN10e2")
# gapGle3 = load_object("Data/GlaubLoc/GlaubLocF0.1OBCgapN10e3")
# gapGle1 = load_object("Data/Glaub/GlaubF0.1OBCgapN10e1")
# gapGle2 = load_object("Data/Glaub/GlaubF0.1OBCgapN10e2")
# gapGle3 = load_object("Data/Glaub/GlaubF0.1OBCgapN10e3")
gapGle1 = load_object("Data/MHLoc/MHLocLocF0.1OBCgapN10e1")
gapGle2 = load_object("Data/MHLoc/MHLocLocF0.1OBCgapN10e2")
gapGle3 = load_object("Data/MHLoc/MHLocLocF0.1OBCgapN10e3")

# plt.scatter(temp, gapGl3, label = L"Glaubler Local $1-e[3]$")
# # plt.scatter(temp, gapGl, label = L"Glaubler Local $1-e[2]$")
# plt.scatter(temp, gapMH, label = L"MH Uniform $1-e[2]$")
# plt.scatter(temp, gapMH3, label = L"MH Uniform $1-e[3]$", alpha = 0.5)
# # plt.scatter(temp, gapMHl, label = L"MH Local $1-e[2]$")
# plt.scatter(temp, gapMHl3, label = L"MH Local $1-e[3]$", alpha = 0.5)

plt.scatter(temp, 1 .- abs.(gapGle1), label = L"$e[1]$")
plt.scatter(temp, 1 .- abs.(gapGle2), label = L"$e[2]$")
plt.scatter(temp, 1 .- abs.(gapGle3), label = L"$e[3]$")


# plt.scatter(temp, gap, label = L"qHMC $\gamma = 0.75$ $t=8$")
# plt.scatter(temp, gap5, label = L"qHMC $\gamma = 0.5$ $t=8$")
plt.title(string(L"MH Local Eigenvalues OBC  with $h=0.1$ $N= $ ", 10))
plt.ylabel("1-Eigenvalue")
plt.xlabel(L"$T$")
plt.yscale("log")
plt.xscale("log")
plt.grid("both","major")
plt.legend(loc=4)
# plt.savefig("Figures/ParamGapN10.png")
plt.savefig("Figures/MHLocal/MHLocEigenvaluesF0.1OBC.png")
# plt.show()