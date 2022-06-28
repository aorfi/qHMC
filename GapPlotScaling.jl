using LaTeXStrings
using Arpack
using PyPlot
using JLD2
using LsqFit


# plot gap vs N for fixed beta
N_values = (2:12)
N_values13 = (2:13)
gap = load_object("Data/qHMC/OBC/gamma0.75t5gapBeta300")
gapMH = load_object("Data/qHMC/OBC/gamma0.75t5gapBeta300NotNormalizedHmix")


# gapB6t = load_object("Data/GlaubLoc/GlaubLocOBCgapBeta6third")
# gapB300 = load_object("Data/GlaubLoc/GlaubLocOBCgapBeta300")
# gapB300t = load_object("Data/GlaubLoc/GlaubLocOBCgapBeta300third")
# gapB01 = load_object("Data/GlaubLoc/GlaubLocOBCgapBeta0.1")
# gapB01t = load_object("Data/GlaubLoc/GlaubLocOBCgapBeta0.1third")

@. model(x,p) = p[1]*x+p[2]
p0 = [-0.2,-12.0]
fit = curve_fit(model,N_values,log.(gap),p0)
print("\nParams: \n ", fit.param)
print("\nError: \n ", stderror(fit))
param = fit.param
scale = round(param[1]/log(2), digits=4)

fitMH = curve_fit(model,N_values,log.(gapMH),p0)
paramMH = fitMH.param
scaleMH = round(paramMH[1]/log(2), digits=4)

x = range(2,13, length= 1000)
plt.scatter(N_values, gap, label = L"$\gamma = 0.75$ $t=5$ ")
plt.plot(x,exp.(model(x,param)),linestyle = "dashed", label = string(L"$2^{kN}$ with $k=$ ",scale))
plt.scatter(N_values, gapMH, label = "MH Uniform")
plt.plot(x,exp.(model(x,paramMH)),linestyle = "dashed", label = string(L"$2^{kN}$ with $k=$ ",scaleMH))
# plt.scatter(N_values, gapB6t, alpha = 0.5)
# plt.plot(x,exp.(model(log.(x),paramB6t)),linestyle = "dashed", label = string(L"$\beta= 6$ 1-e[3], scaling = ", scaleB6t))


# plt.scatter(N_values, gapB01, alpha = 0.5)
# plt.plot(x,exp.(model(log.(x),paramB01)),linestyle = "dashed", label = string(L"$\beta= 0.1$ 1-e[2], scaling = ", scaleB01))
# plt.scatter(N_values, gapB01t, label=L"$\beta= 0.1$ 1-e[3]", marker = "^")
# plt.scatter(N_values, gapB300, label=L"$\beta= 300$ 1-e[2]", alpha = 0.5)
# plt.scatter(N_values, gapB300t, label=L"$\beta= 300$ 1-e[3]", marker = "^")


plt.title(L"qHMC Spectral Gap Scaling $\beta = 300 $")
plt.ylabel(L"$\delta$")
plt.xlabel(L"$N$")
plt.yscale("log")
# plt.xscale("log")
plt.grid("both","both")
plt.legend()
# plt.savefig("Figures/qHMCScaling/gapScalingBeta2.png")
plt.show()




