using LaTeXStrings
using Arpack
using PyPlot
using JLD2
using LsqFit

beta = 5
@. model(x,p) = p[1]*x+p[2]
p0 = [-0.2,-12.0]


gapG = load_object("Data/SK/Classical/av100GlaubBeta5")
gap = load_object("Data/SK/qHMC/alpha1.5eta2.5beta5")
gapAv = load_object("Data/SK/qHMC/20random0-10Beta5")
gapT = load_object("Data/SK/qHMC/20random0-10Beta5Triangle")

N_values = (5:12)
N_values_short = (5:10)

fitG = curve_fit(model,N_values,log.(gapG),p0)
paramG = fitG.param
scaleG = -round(paramG[1]/log(2), digits=4)

fit = curve_fit(model,N_values,log.(gap),p0)
param = fit.param
scale = -round(param[1]/log(2), digits=4)

fitAv = curve_fit(model,N_values_short[1:end-1],log.(gapAv[1:end-1]),p0)
paramAv = fitAv.param
scaleAv = -round(paramAv[1]/log(2), digits=4)

fitT = curve_fit(model,N_values_short,log.(gapT),p0)
paramT = fitT.param
scaleT = -round(paramT[1]/log(2), digits=4)

x = range(5,15, length= 1000)


label_scatter = L"qHMC $\alpha = 1.5$ $\eta=2.5$ $2^{-kN}$ k= "*string(scale)
plt.scatter(N_values, gap, color = "tab:green")
plt.plot(x,exp.(model(x,param)),linestyle = "dashed", label = label_scatter,color = "tab:green")

label_scatterG = L"Uniform Glauber $2^{-kN}$ k= "*string(scaleG)
plt.scatter(N_values , gapG, color = "tab:blue")
plt.plot(x,exp.(model(x,paramG)),linestyle = "dashed", label = label_scatterG,color = "tab:blue")

label_scatter = L"qHMC $\alpha,\eta\in[0,10]$ $2^{-kN}$ k= "*string(scaleAv)
plt.scatter(N_values_short, gapAv, color = "tab:orange")
plt.plot(x,exp.(model(x,paramAv)),linestyle = "dashed", label = label_scatter,color = "tab:orange")

label_scatter = L"qHMC $\alpha\in[0,10], \eta<\alpha$ $2^{-kN}$ k= "*string(scaleT)
plt.scatter(N_values_short, gapT, color = "tab:purple")
plt.plot(x,exp.(model(x,paramT)),linestyle = "dashed", label = label_scatter,color = "tab:purple")

plt.title(L"SK qHMC Gap Scaling Comparison $\beta = $"*string(beta))
plt.ylabel(L"$\delta$")
plt.xlabel(L"$N$")
plt.yscale("log")
# plt.xscale("log")
plt.grid("both","both")
plt.legend(loc=3)
name = "Figures/SK/avGapBeta5Random.png"
plt.savefig(name)
plt.show()