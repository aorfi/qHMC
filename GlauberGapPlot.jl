using LaTeXStrings
using Arpack
using PyPlot
using JLD2
using LsqFit


# plot gap vs N for fixed beta
N_values = (2:15)
gap = load_object("Data/gapBeta6")

@. model(x,p) = p[1]*x+p[2]
p0 = [-1.0,4.0]
fit = curve_fit(model,log.(N_values),log.(gap),p0)
param = fit.param
sigma = stderror(fit)
x = range(2,15, length= 1000)
print("Params: \n ", param)
print("\nError: \n ", sigma)
scale = round(param[1], digits=3)

plt.scatter(N_values, gap)
plt.plot(x,exp.(model(log.(x),param)),linestyle = "dashed", label = string("Scaling = ", scale))
plt.title(string(L"Glauber Dynamics Spectral Gap  $\beta= $ ", 6))
plt.ylabel(L"$\delta$")
plt.xlabel(L"$N$")
plt.yscale("log")
plt.xscale("log")
plt.grid("both","both")
plt.legend()
# plt.savefig("Figures/gapBeta6.png")
plt.show()




# # plot gap vs temp for fixed N
# N = 10
# temp = 10 .^ (range(-2.5,stop=2.5,length=50))
# gap = load_object("Data/gapN10")
# plt.scatter(temp, gap)
# plt.title(string(L"Glauber Dynamics Spectral Gap  $N= $ ", 10))
# plt.ylabel(L"$\delta$")
# plt.xlabel(L"$T$")
# plt.yscale("log")
# plt.xscale("log")
# plt.grid("both","major")
# # plt.savefig("Figures/gapN10.png")
# plt.show()