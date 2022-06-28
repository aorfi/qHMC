using LaTeXStrings
using Arpack
using PyPlot
using JLD2
using LsqFit

# plot gap vs temp for fixed N
N = 10
temp = 10 .^ (range(-2.5,stop=2.5,length=50))
gapMH = load_object("Data/MH/MHgapOBCN10")
# gapMH3 = load_object("Data/MH/MHgapOBCN10third")
gapMHl = load_object("Data/MHLoc/MHLocGapOBCN10")
# gapMHl3 = load_object("Data/MHLoc/MHLocGapOBCN10third")
gapG = load_object("Data/Glaub/GlaubGapOBCN10")
gapGl = load_object("Data/GlaubLoc/GlaubLocOBCgapN10")
# gapGl3 = load_object("Data/GlaubLoc/GlaubLocOBCgapN10third")


plt.scatter(temp, gapMH, label = "MH Uniform")
plt.scatter(temp, gapMHl, label = "MH Local")
plt.scatter(temp, gapG, label = "Glaubler Uniform")
plt.scatter(temp, gapGl, label = "Glaubler Local")


plt.title(string(L"Classical MCMC with OBC $N= $ ", 10))
plt.ylabel(L"$\delta$")
plt.xlabel(L"$T$")
plt.yscale("log")
plt.xscale("log")
plt.grid("both","major")
plt.legend(loc=4)
plt.savefig("Figures/ClassicalOBCGapN10.png")
# plt.savefig("Figures/MHLocal/MHLocEigenvaluesF0.1OBC.png")
plt.show()