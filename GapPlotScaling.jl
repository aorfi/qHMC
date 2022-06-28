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
# gapB6 = load_object("Data/GlaubLoc/GlaubLocOBCgapBeta6")
# gap = load_object("Data/GlaubLoc/GlaubLocPBCgapBeta6")
gapB6 = load_object("Data/Glaub/OBCGlaubGengapBeta300")
gap = load_object("Data/MH/OBCMHgapBeta300")


# gapB6t = load_object("Data/GlaubLoc/GlaubLocOBCgapBeta6third")
# gapB300 = load_object("Data/GlaubLoc/GlaubLocOBCgapBeta300")
# gapB300t = load_object("Data/GlaubLoc/GlaubLocOBCgapBeta300third")
# gapB01 = load_object("Data/GlaubLoc/GlaubLocOBCgapBeta0.1")
# gapB01t = load_object("Data/GlaubLoc/GlaubLocOBCgapBeta0.1third")

@. model(x,p) = p[1]*x+p[2]
p0 = [-0.2,-12.0]
fit = curve_fit(model,N_values,log.(gapB6),p0)
print("\nParams: \n ", fit.param)
print("\nError: \n ", stderror(fit))
param = fit.param


fitp = curve_fit(model,N_values,log.(gap),p0)
paramp = fitp.param


scale = round(param[1]/log(2), digits=4)
scalep = round(paramp[1]/log(2), digits=4)

# fitB6t = curve_fit(model,log.(N_values),log.(gapB6t),p0)
# paramB6t = fitB6t.param
# scaleB6t = round(paramB6t[1], digits=4)

# fitB01 = curve_fit(model,log.(N_values),log.(gapB01),p0)
# paramB01 = fitB01.param
# scaleB01 = round(paramB01[1], digits=4)

x = range(2,13, length= 1000)
plt.scatter(N_values, gapB6)
plt.plot(x,exp.(model(x,param)),linestyle = "dashed", label = L"Glauber Uniform $2^{-N}$ ")
plt.scatter(N_values, gap)
plt.plot(x,exp.(model(x,paramp)),linestyle = "dashed", label = L"MH Uniform $2^{-N}$ ")
# plt.scatter(N_values, gapB6t, alpha = 0.5)
# plt.plot(x,exp.(model(log.(x),paramB6t)),linestyle = "dashed", label = string(L"$\beta= 6$ 1-e[3], scaling = ", scaleB6t))


# plt.scatter(N_values, gapB01, alpha = 0.5)
# plt.plot(x,exp.(model(log.(x),paramB01)),linestyle = "dashed", label = string(L"$\beta= 0.1$ 1-e[2], scaling = ", scaleB01))
# plt.scatter(N_values, gapB01t, label=L"$\beta= 0.1$ 1-e[3]", marker = "^")
# plt.scatter(N_values, gapB300, label=L"$\beta= 300$ 1-e[2]", alpha = 0.5)
# plt.scatter(N_values, gapB300t, label=L"$\beta= 300$ 1-e[3]", marker = "^")


plt.title(L"Uniform Spectral Gap Scaling OBC $\beta = 300 $")
plt.ylabel(L"$\delta$")
plt.xlabel(L"$N$")
plt.yscale("log")
# plt.xscale("log")
plt.grid("both","both")
plt.legend()
# plt.savefig("Figures/MHLocal/gapBCscaling.png")
plt.show()




