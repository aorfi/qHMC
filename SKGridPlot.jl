using PyPlot
using LaTeXStrings
using JLD2

N=5
beta = 5
num_values = 10
name =  "Data/SK/GridSearch/100av"*string(num_values)*"N"*string(N)*"beta"*string(beta)
gap_all= load_object(name)
max,cord = findmax(gap_all)
max_param = 50
plt.title(L"qHMC Spectral Gap $\beta= $ "*string(beta)*L" $N = $"*string(N))
# plt.imshow(gap_all,extent = [0, 30, 0 , 30], aspect="auto")
plt.imshow(gap_all, origin="lower",extent = [0, max_param, 0 , max_param])
bar = plt.colorbar()
plt.scatter(cord[2]/num_values*max_param, cord[1]/num_values*max_param, color="red", marker=".")
plt.xlabel(L"$\eta$")
plt.ylabel(L"$\alpha$")

# plt.clim(0, 0.05) 
bar.set_label(L"$\delta$")
# name = "Figures/SK/GridSearch/"*string(num_values)*"N"*string(N)*"beta"*string(beta)*"i.png"
# plt.savefig(name)
# plt.clf()
plt.show()

