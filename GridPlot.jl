using PyPlot
using LaTeXStrings

N=5
beta = 6
num_values = 10
name = name = "Data/GridSearch/alphaEtaParam/"*string(num_values)*"N"*string(N)*"beta"*string(beta)
gap_all= load_object(name)

plt.title(L"qHMC Spectral Gap $\beta= $ "*string(beta)*L" $N = $"*string(N))
plt.imshow(gap_all,extent = [0, 30, 0 , 30], aspect="auto")
plt.xlabel(L"$\eta$")
plt.ylabel(L"$\alpha$")
bar = plt.colorbar()
# plt.clim(0, 0.05) 
bar.set_label(L"$\delta$")
name = "Figures/GridSearch/alphaEtaParam/"*string(num_values)*"N"*string(N)*"beta"*string(beta)*".png"
plt.savefig(name)
plt.show()