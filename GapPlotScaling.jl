using LaTeXStrings
using Arpack
using PyPlot
using JLD2
using LsqFit


# plot gap vs N for fixed beta
# N_values = (2:13)
# gapMH = load_object("Data/MHgapBeta300")
# gapG = load_object("Data/GlaubGengapBeta300")
# @. model(x,p) = p[1]*x+p[2]
# p0 = [-1.0,4.0]
# N_values_even = N_values[1:2:end]
# gapG_even = gapG[1:2:end]
# gapMH_even = gapMH[1:2:end]
# fit = curve_fit(model,N_values_even,log.(gapG_even),p0)
# print("\nParams: \n ", fit.param)
# print("\nError: \n ", stderror(fit))
# paramG = fit.param
# fit = curve_fit(model,N_values_even,log.(gapMH_even),p0)
# print("\nParams: \n ", fit.param)
# print("\nError: \n ", stderror(fit))
# paramMH = fit.param

# scaleG = round(paramG[1]/log(2), digits=4)
# scaleMH = round(paramMH[1]/log(2), digits=4)

# x = range(2,15, length= 1000)
# plt.scatter(N_values, gapG)
# plt.scatter(N_values, gapMH)
# plt.plot(x,exp.(model(x,paramG)),linestyle = "dashed", label = string("Glauber Scaling k= ", scaleG))
# plt.plot(x,exp.(model(x,paramMH)),linestyle = "dashed", label = string("MH Scaling k= ", scaleMH))
# plt.title(string(L"Spectral Gap  $\beta= $ ", 300))
# plt.ylabel(L"$\delta$")
# plt.xlabel(L"$N$")
# plt.yscale("log")
# # plt.xscale("log")
# plt.grid("both","both")
# plt.legend()
# plt.savefig("Figures/gapBeta300fit.png")
# plt.show()

N_values = (2:13)
gap = load_object("Data/gamma0.75t8gapBeta300")
gap5 = load_object("Data/gamma0.5t8gapBeta300")
# @. model(x,p) = p[1]*x+p[2]
# p0 = [-1.0,4.0]
# N_values_even = N_values[1:2:end]
# gapG_even = gapG[1:2:end]
# gapMH_even = gapMH[1:2:end]
# fit = curve_fit(model,N_values_even,log.(gapG_even),p0)
# print("\nParams: \n ", fit.param)
# print("\nError: \n ", stderror(fit))
# paramG = fit.param
# fit = curve_fit(model,N_values_even,log.(gapMH_even),p0)
# print("\nParams: \n ", fit.param)
# print("\nError: \n ", stderror(fit))
# paramMH = fit.param

# scaleG = round(paramG[1]/log(2), digits=4)
# scaleMH = round(paramMH[1]/log(2), digits=4)

x = range(2,15, length= 1000)
plt.scatter(N_values, gap)
plt.scatter(N_values, gap5)
# plt.plot(x,exp.(model(x,paramG)),linestyle = "dashed", label = string("Glauber Scaling k= ", scaleG))
# plt.plot(x,exp.(model(x,paramMH)),linestyle = "dashed", label = string("MH Scaling k= ", scaleMH))
plt.title(string(L"Spectral Gap  $\beta= $ ", 300))
plt.ylabel(L"$\delta$")
plt.xlabel(L"$N$")
plt.yscale("log")
# plt.xscale("log")
plt.grid("both","both")
plt.legend()
plt.savefig("Figures/gapBeta300fit.png")
plt.show()




