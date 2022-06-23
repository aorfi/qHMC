using LaTeXStrings
using Arpack
using PyPlot
using JLD2
using LsqFit

# plot gap vs temp for fixed N
N = 10
temp = 10 .^ (range(-2.5,stop=2.5,length=50))
gapMH = load_object("Data/MHgapN10")
gapG = load_object("Data/GlaubGenGapN10")
gap = load_object("Data/gamma0.75t8gapN10")
plt.scatter(temp, gapG, label = "Glauber")
plt.scatter(temp, gapMH, label = "MH")
plt.scatter(temp, gap, label = L"qHMC $\gamma = 0.75$ $t=8$")
plt.title(string(L"Spectral Gap  $N= $ ", 10))
plt.ylabel(L"$\delta$")
plt.xlabel(L"$T$")
plt.yscale("log")
plt.xscale("log")
plt.grid("both","major")
plt.legend()
plt.savefig("Figures/gapN10.png")
plt.show()