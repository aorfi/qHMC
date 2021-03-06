using LaTeXStrings
using Arpack
using PyPlot
using JLD2
using LsqFit

# # # plot gap vs N for fixed beta
N_values = (5:12)
beta = 5
gapMH = load_object("Data/SK/Classical/av100MHBeta"*string(beta))
gapMHl = load_object("Data/SK/Classical/av100MHLocBeta"*string(beta))
gapG = load_object("Data/SK/Classical/av100GlaubBeta"*string(beta))
gapGl = load_object("Data/SK/Classical/av100GlaubLocBeta"*string(beta))
gap = load_object("Data/SK/qHMC/alpha1.5eta2.5beta5")
@. model(x,p) = p[1]*x+p[2]
p0 = [-0.2,-12.0]

fit = curve_fit(model,N_values,log.(gap),p0)
param = fit.param
scale = -round(param[1]/log(2), digits=4)

fitG = curve_fit(model,N_values,log.(gapG),p0)
paramG = fitG.param
scaleG = -round(paramG[1]/log(2), digits=4)

fitMH = curve_fit(model,N_values,log.(gapMH),p0)
paramMH = fitMH.param
scaleMH = -round(paramMH[1]/log(2), digits=4)

# fitGl = curve_fit(model,log.(N_values),log.(gapGl),p0)
# paramGl = fitGl.param
# scaleGl = round(paramGl[1], digits=4)

# fitMHl = curve_fit(model,log.(N_values),log.(gapMHl),p0)
# paramMHl = fitMHl.param
# scaleMHl = round(paramMHl[1], digits=4)

x = range(5,13, length= 1000)
plt.scatter(N_values, gapG)
labelG = L"Glauber Uniform with $2^{-kN}$ with $k=$"*string(scaleG)
plt.plot(x,exp.(model(x,paramG)),linestyle = "dashed", label = labelG)

plt.scatter(N_values, gapMH)
labelMH = L"MH Uniform with $2^{-kN}$ with $k=$"*string(scaleMH)
plt.plot(x,exp.(model(x,paramMH)),linestyle = "dashed", label = labelMH)
# plt.scatter(N_values, gapMHl, label = "MH Local")
# labelMHl = L"MH Local with $N^{b}$ b= "*string(scaleMHl)
# plt.plot(x,exp.(model(log.(x),paramMHl)),linestyle = "dashed", label = labelMHl)
# plt.scatter(N_values, gapGl, label = "Glaub Local")
# labelGl = L"Glauber Local with $N^{b}$ b= "*string(scaleGl)
# plt.plot(x,exp.(model(log.(x),paramGl)),linestyle = "dashed", label = labelGl)
plt.scatter(N_values, gap)
label = L"qHMC $\alpha=1.5$ $\eta=2.5$ $2^{-kN}$ with $k=$"*string(scale)
plt.plot(x,exp.(model(x,param)),linestyle = "dashed", label = label)

plt.title(L"SK Gap Scaling Comparison $\beta = $"*string(beta))
plt.ylabel(L"Average $\delta$")
plt.xlabel(L"$N$")
plt.yscale("log")
# plt.xscale("log")
plt.grid("both","both")
plt.legend(loc=3)
plt.savefig("Figures/SK/avGapBeta"*string(beta)*".png")
plt.show()
