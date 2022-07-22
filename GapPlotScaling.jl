using LaTeXStrings
using Arpack
using PyPlot
using JLD2
using LsqFit

# Get best gap values
# best_gaps = zeros(6)
# beta = 6
# num_values = 200
# for N in (5:9)
#     name =  "Data/GridSearch/alphaEtaParam/"*string(num_values)*"N"*string(N)*"beta"*string(beta)
#     gap_all= load_object(name)
#     max,cord = findmax(gap_all)
#     best_gaps[N-4] = max
# end
# name =  "Data/GridSearch/alphaEtaParam/100N10beta"*string(beta)
# gap_all= load_object(name)
# max,cord = findmax(gap_all)
# best_gaps[6] = max
# name = "Data/qHMC/alphaEtaParam/bestGapBeta"*string(beta)
# save_object(name, best_gaps)



# Get Scaling
# beta = 6
# N_values = (4:10)
# alpha = 0
# num_values = 500
# eta_values = range(0,10, length=num_values)
# gap_all = load_object("Data/qHMC/alphaEtaParam/alpha"*string(alpha)*"scalingBeta6")

# @. model(x,p) = p[1]*x+p[2]
# p0 = [-0.2,-12.0]

# scaleEx =  Union{Float64, Nothing}[nothing for _ in 1:num_values]
# scalePol = Union{Float64, Nothing}[nothing for _ in 1:num_values]
# scaleExAll = zeros(num_values)
# scalePolAll = zeros(num_values)
# for i in (2:num_values)
#     gap = gap_all[3:end,i]
#     fitEx = curve_fit(model,N_values,log.(gap),p0)
#     paramEx = fitEx.param
#     fitPol = curve_fit(model,log.(N_values),log.(gap),p0)
#     paramPol = fitPol.param
#     errorEx = stderror(fitEx)[1]
#     errorPol = stderror(fitPol)[1]
#     println("Ex Error", errorEx)
#     println("Pol Error", errorPol)
#     scalePolAll[i] = -paramPol[1]
#     scaleExAll[i] = -paramEx[1]/log(2)
#     if errorPol < 0.1
#         scalePol[i] = -paramPol[1]
#     end
#     if errorEx < 0.1
#         scaleEx[i] = -paramEx[1]/log(2)
#     end
# end
# name = "Data/GridSearch/alphaEtaParam/AlphaZero/"*string(num_values)*"N8beta"*string(beta)
# gapN8 = load_object(name)
# fig, axs = plt.subplots(2)
# fig.suptitle(L"Scaling $\alpha=0$ $\beta=6$")


# axs[1].plot(eta_values,scaleExAll, color = "tab:blue", linestyle = "dashed")
# axs[1].plot(eta_values,scaleEx, color = "tab:blue")
# axs[2].plot(eta_values,scalePolAll, color = "tab:green", linestyle = "dashed")
# axs[2].plot(eta_values,scalePol, color = "tab:green")
# # ax2.plot(eta_values,gapN8, color = "tab:orange", linestyle= "dashed")
# axs[2].set_ylabel("Polynominal Scaling", color = "tab:green")
# axs[2].tick_params(axis="y", colors="tab:green")
# axs[1].set_ylabel("Exponential Scaling", color = "tab:blue")
# axs[1].tick_params(axis="y", colors="tab:blue")
# axs[2].set_xlabel(L"$\eta$")
# axs[1].grid()
# axs[2].grid()
# name = "Figures/Ising/GapInvestigation/AlphaZeroAllScalingBoth.png"
# plt.savefig(name)
# plt.show()






beta = 6
@. model(x,p) = p[1]*x+p[2]
p0 = [-0.2,-12.0]

gapGl = load_object("Data/GlaubLoc/OBCBeta6")
gapG = load_object("Data/Glaub/OBCBeta6")
gap = load_object("Data/qHMC/alphaEtaParam/alpha0etaPI+0.01beta6e2")
gap1 = load_object("Data/qHMC/alphaEtaParam/50random0-30TriangleBeta6e2")
gap2 = load_object("Data/qHMC/alphaEtaParam/50random0-30TriangleBeta6e3")

N_values = (5:10)
N_values_big = (5:13)

# fit = curve_fit(model,N_values,log.(gap),p0)
# param = fit.param
# scale = -round(param[1]/log(2), digits=4)

# fit = curve_fit(model,log.(N_values),log.(gap),p0)
# param = fit.param
# scale = round(param[1], digits=4)

# fit1 = curve_fit(model,log.(N_values),log.(gap1),p0)
# param1 = fit1.param
# scale1 = round(param1[1], digits=4)


# fit2 = curve_fit(model,log.(N_values),log.(gap2),p0)
# param2 = fit2.param
# scale2 = round(param2[1], digits=4)

fit1 = curve_fit(model,N_values[1:end-1],log.(gap1[1:end-1]),p0)
param1 = fit1.param
scale1 = -round(param1[1]/log(2), digits=4)

fit2 = curve_fit(model,N_values,log.(gap2),p0)
param2 = fit2.param
scale2 = -round(param2[1]/log(2), digits=4)

fitG = curve_fit(model,N_values_big,log.(gapG),p0)
paramG = fitG.param
scaleG = -round(paramG[1]/log(2), digits=4)

fitGl = curve_fit(model,log.(N_values_big),log.(gapGl),p0)
paramGl = fitGl.param
scaleGl = round(paramGl[1], digits=4)


x = range(5,15, length= 1000)

# label_scatter1 = L"qHMC $\alpha,\eta\in[0,10]$ $2^{-kN}$ k= "*string(scale1)
# plt.scatter(N_values, gap1, color = "tab:green")
# plt.plot(x,exp.(model(x,param1 )),linestyle = "dashed", label = label_scatter1 ,color = "tab:green")


# label_scatterPer1 = L"qHMC $\alpha = 1$ $\eta=0.1$ $N^{b}$ b= "*string(scale1)
# plt.scatter(N_values, gapPer1)
# plt.plot(x,exp.(model(log.(x),param1 )),linestyle = "dashed", label = label_scatterPer1 )


# label_scatter = L"$\alpha = 0$ $\eta = \pi + 0.01$ e[2] $N^{b}$ b= "*string(scale)
# plt.scatter(N_values, gap, color = "tab:blue")
# plt.plot(x,exp.(model(log.(x),param)),linestyle = "dashed", label = label_scatter,color = "tab:blue")

# label_scatter = L"$\alpha = 0$ $\eta = \pi + 0.01$ e[3] $N^{b}$ b= "*string(scale1)
# plt.scatter(N_values, gap1, color = "tab:orange")
# plt.plot(x,exp.(model(log.(x),param1)),linestyle = "dashed", label = label_scatter,color = "tab:orange")

# label_scatter1 = L"$\alpha = 0$ $\eta = \pi/2 + 0.01$ $N^{b}$ b= "*string(scale1)
# plt.scatter(N_values, gap1, color = "tab:red")
# plt.plot(x,exp.(model(log.(x),param1)),linestyle = "dashed", label = label_scatter1,color = "tab:red")


# label_scatter2 = L"$\alpha = 0$ $\eta = \pi/2 + 0.001$ $N^{b}$ b= "*string(scale2)
# plt.scatter(N_values, gap2, color = "tab:purple")
# plt.plot(x,exp.(model(log.(x),param2)),linestyle = "dashed", label = label_scatter2,color = "tab:purple")

label_scatterG = L"Uniform Glauber $2^{-kN}$ k= "*string(scaleG)
plt.scatter(N_values_big , gapG, color = "tab:orange")
plt.plot(x,exp.(model(x,paramG)),linestyle = "dashdot", label = label_scatterG,color = "tab:orange")

label_scatter1 = L"qHMC $\alpha\in[0,10]$ $\eta<\alpha$ e[2] $2^{-kN}$ k= "*string(scale1)
plt.scatter(N_values, gap1, color = "tab:green")
plt.plot(x,exp.(model(x,param1 )),linestyle = "dashed", label = label_scatter1 ,color = "tab:green")

label_scatter2 = L"qHMC $\alpha\in[0,10]$ $\eta<\alpha$ e[3] $2^{-kN}$ k= "*string(scale2)
plt.scatter(N_values, gap2, color = "tab:red")
plt.plot(x,exp.(model(x,param2 )),linestyle = "dashed", label = label_scatter2 ,color = "tab:red")


# label_scatterGl = L"Glauber Local $N^{b}$ b= "*string(scaleGl)
# plt.scatter(N_values_big ,  gapGl, color = "tab:green")
# plt.plot(x,exp.(model(log.(x),paramGl)),linestyle = "dashed", label = label_scatterGl,color = "tab:green")




plt.title(L"Gap Scaling Comparison $\beta = $"*string(beta))
plt.ylabel(L"$\delta$")
plt.xlabel(L"$N$")
plt.yscale("log")
# plt.xscale("log")
plt.grid("both","both")
plt.legend()
# name = "Figures/Ising/qHMCScaling/alphaEtaParam/alpha0etaPIe.png"
# plt.savefig(name)
plt.show()