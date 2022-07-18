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
beta = 6
N_values = (4:8)
alpha = 0
num_values = 500
eta_values = range(0,10, length=num_values)
gap_all = load_object("Data/qHMC/alphaEtaParam/alpha"*string(alpha)*"scalingBeta6")

@. model(x,p) = p[1]*x+p[2]
p0 = [-0.2,-12.0]

scale = zeros(num_values)
for i in (2:num_values)
    gap = gap_all[:,i]
    fit = curve_fit(model,N_values,log.(gap),p0)
    param = fit.param
    scale[i] = -param[1]/log(2)
    # print("\nParams: \n ", fit.param)
    error = stderror(fit)[1]
    if error >= 0.5
        println( "Warning fit error of ", error)
        println(" i: ", i)
    end
end
name = "Data/GridSearch/alphaEtaParam/AlphaZero/"*string(num_values)*"N8beta"*string(beta)
gapN8 = load_object(name)


fig, ax1 = plt.subplots()
plt.title(L"Gap and Scaling $\alpha=0$")
plt.xlabel(L"$\eta$")
ax2 = ax1.twinx()
ax1.plot(eta_values[2:end],scale[2:end], color = "tab:blue")
ax2.plot(eta_values,gapN8, color = "tab:orange", linestyle= "dashed")
ax2.set_ylabel("Gap N8", color = "tab:orange")
ax2.tick_params(axis="y", colors="tab:orange")
ax1.set_ylabel("Scaling Exponent", color = "tab:blue")
ax1.tick_params(axis="y", colors="tab:blue")
plt.grid()
plt.show()






# beta = 6
# @. model(x,p) = p[1]*x+p[2]
# p0 = [-0.2,-12.0]

# gapGl = load_object("Data/GlaubLoc/Old/CmpareGlaubLocOBCgapBeta6")
# gapG = load_object("Data/Glaub/Old/CompareOBCGlaubGengapBeta6")
# gap = load_object("Data/qHMC/alphaEtaParam/alpha0eta3.14beta6")
# # gap1 = load_object("Data/qHMC/alphaEtaParam/alpha0eta1.63beta6")
# # gap2 = load_object("Data/qHMC/alphaEtaParam/alpha0eta2.1beta6")
# display(gap)
# N_values = (5:12)
# N_values_big = (5:15)

# fit = curve_fit(model,N_values,log.(gap),p0)
# param = fit.param
# scale = -round(param[1]/log(2), digits=4)

# # fit1 = curve_fit(model,N_values,log.(gap1),p0)
# # param1 = fit1.param
# # scale1 = -round(param1[1]/log(2), digits=4)

# # fit2 = curve_fit(model,N_values,log.(gap2),p0)
# # param2 = fit2.param
# # scale2 = -round(param2[1]/log(2), digits=4)

# fitG = curve_fit(model,N_values_big,log.(gapG),p0)
# paramG = fitG.param
# println(paramG)
# scaleG = -round(paramG[1]/log(2), digits=4)

# x = range(5,15, length= 1000)

# label_scatter = L"qHMC $\alpha = 0$ $\eta\sim\pi$ $2^{-kN}$ k= "*string(scale)
# plt.scatter(N_values, gap, color = "tab:red")
# plt.plot(x,exp.(model(x,param )),linestyle = "dashed", label = label_scatter ,color = "tab:red")

# # label_scatter1 = L"qHMC $\alpha = 0$ $\eta=1.63$ $2^{-kN}$ k= "*string(scale1)
# # plt.scatter(N_values, gap1, color = "tab:blue")
# # plt.plot(x,exp.(model(x,param1 )),linestyle = "dashed", label = label_scatter1 ,color = "tab:blue")


# # label_scatterPer1 = L"qHMC $\alpha = 1$ $\eta=0.1$ $N^{b}$ b= "*string(scale1)
# # plt.scatter(N_values, gapPer1)
# # plt.plot(x,exp.(model(log.(x),param1 )),linestyle = "dashed", label = label_scatterPer1 )

# # plt.scatter(N_values, gapPer4, label = L"qHMC $\alpha = 0$ $\eta=1$")
# # label_scatter = L"Glauber qHMC Best $2^{-kN}$ k= "*string(scale)
# # plt.scatter(N_values, gap, color = "tab:blue")
# # plt.plot(x,exp.(model(x,param)),linestyle = "dashed", label = label_scatter,color = "tab:blue")

# label_scatterG = L"Uniform Glauber $2^{-kN}$ k= "*string(scaleG)
# plt.scatter(N_values_big , gapG, color = "tab:orange")
# plt.plot(x,exp.(model(x,paramG)),linestyle = "dashdot", label = label_scatterG,color = "tab:orange")

# # label_scatter2 = L"qHMC $\alpha = 0$ $\eta=2.1$ $2^{-kN}$ k= "*string(scale2)
# # plt.scatter(N_values, gap2, color = "tab:green")
# # plt.plot(x,exp.(model(x,param2 )),linestyle = "dashed", label = label_scatter2 ,color = "tab:green")


# # label_scatterGl = L"Local Glauber $N^{b}$ b= "*string(scaleGl)
# # plt.scatter(N_values_big , gapGl, color = "tab:green")
# # plt.plot(x,exp.(model(log.(x),paramGl)),linestyle = "dashed", label = label_scatterGl,color = "tab:green")



# plt.title(L"qHMC Gap Scaling Comparison $\beta = $"*string(beta))
# plt.ylabel(L"$\delta$")
# plt.xlabel(L"$N$")
# plt.yscale("log")
# # plt.xscale("log")
# plt.grid("both","both")
# plt.legend(loc=1)
# name = "Figures/Ising/GapInvestigation/AlphaZeroEtaPiScaling.png"
# plt.savefig(name)
# plt.show()