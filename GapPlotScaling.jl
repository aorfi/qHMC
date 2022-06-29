# using LaTeXStrings
# using Arpack
# using PyPlot
# using JLD2
# using LsqFit


# # plot gap vs N for fixed beta
# N_values = (2:12)
# N_values13 = (2:13)
# gap = load_object("Data/qHMC/SigmaX/gamma0.5t10Beta300")
# gap3 = load_object("Data/qHMC/SigmaX/gamma0.3t10Beta300")
# gap1 = load_object("Data/qHMC/SigmaX/gamma0.1t10Beta300")
# gap6 = load_object("Data/qHMC/SigmaX/gamma0.6t10Beta300")
# gapG = load_object("Data/Glaub/OBCGlaubGengapBeta300")


# # gapB6t = load_object("Data/GlaubLoc/GlaubLocOBCgapBeta6third")
# # gapB300 = load_object("Data/GlaubLoc/GlaubLocOBCgapBeta300")
# # gapB300t = load_object("Data/GlaubLoc/GlaubLocOBCgapBeta300third")
# # gapB01 = load_object("Data/GlaubLoc/GlaubLocOBCgapBeta0.1")
# # gapB01t = load_object("Data/GlaubLoc/GlaubLocOBCgapBeta0.1third")

# @. model(x,p) = p[1]*x+p[2]
# p0 = [-0.2,-12.0]
# fit = curve_fit(model,N_values,log.(gap),p0)
# # print("\nParams: \n ", fit.param)
# # print("\nError: \n ", stderror(fit))
# param = fit.param
# scale = round(param[1]/log(2), digits=4)

# fitG = curve_fit(model,N_values13,log.(gapG),p0)
# paramG = fitG.param
# scaleG = round(paramG[1]/log(2), digits=4)

# fit3 = curve_fit(model,N_values,log.(gap3),p0)
# param3 = fit3.param
# scale3 = round(param3[1]/log(2), digits=4)

# fit1 = curve_fit(model,N_values,log.(gap1),p0)
# param1 = fit1.param
# scale1 = round(param1[1]/log(2), digits=4)
# fit6 = curve_fit(model,N_values,log.(gap1),p0)
# param6 = fit6.param
# scale6 = round(param6[1]/log(2), digits=4)

# x = range(2,13, length= 1000)
# plt.scatter(N_values, gap, label = L"Glauber qHMC $\gamma = 0.5$ $t=10$ ")
# plt.plot(x,exp.(model(x,param)),linestyle = "dashed", label = string(L"$2^{kN}$ with $k=$ ",scale))
# # plt.scatter(N_values13, gapG, label = "Glauber Uniform")
# # plt.plot(x,exp.(model(x,paramG)),linestyle = "dashed", label = string(L"$2^{kN}$ with $k=$ ",scaleG))

# plt.scatter(N_values, gap3, label = L"Glauber qHMC $\gamma = 0.3$ $t=10$ ")
# plt.plot(x,exp.(model(x,param3)),linestyle = "dashed", label = string(L"$2^{kN}$ with $k=$ ",scale3))
# plt.scatter(N_values, gap1, label = L"Glauber qHMC $\gamma = 0.1$ $t=10$ ")
# plt.plot(x,exp.(model(x,param1)),linestyle = "dashed", label = string(L"$2^{kN}$ with $k=$ ",scale1))
# # plt.scatter(N_values, gapB6t, alpha = 0.5)
# # plt.plot(x,exp.(model(log.(x),paramB6t)),linestyle = "dashed", label = string(L"$\beta= 6$ 1-e[3], scaling = ", scaleB6t))


# # plt.scatter(N_values, gapB01, alpha = 0.5)
# # plt.plot(x,exp.(model(log.(x),paramB01)),linestyle = "dashed", label = string(L"$\beta= 0.1$ 1-e[2], scaling = ", scaleB01))
# # plt.scatter(N_values, gapB01t, label=L"$\beta= 0.1$ 1-e[3]", marker = "^")
# # plt.scatter(N_values, gapB300, label=L"$\beta= 300$ 1-e[2]", alpha = 0.5)
# # plt.scatter(N_values, gapB300t, label=L"$\beta= 300$ 1-e[3]", marker = "^")


# plt.title(L"qHMC Gap Scaling Comparison $\beta = 300 $")
# plt.ylabel(L"$\delta$")
# plt.xlabel(L"$N$")
# plt.yscale("log")
# # plt.xscale("log")
# plt.grid("both","both")
# plt.legend()
# # plt.savefig("Figures/qHMCScaling/SigmaXgapScalingBeta300CompareT10.png")
# plt.show()



gamma_values = [0.2,0.3,0.5]
t = 20
beta = 300
@. model(x,p) = p[1]*x+p[2]
p0 = [-0.2,-12.0]

for gamma in gamma_values
    N_values = (2:12)
    name = "Data/qHMC/SigmaX/gamma"*string(gamma)*"t"*string(t)*"Beta"*string(beta)
    gap = load_object(name)

    fit = curve_fit(model,N_values,log.(gap),p0)
    param = fit.param
    scale = -round(param[1]/log(2), digits=4)
    x = range(2,13, length= 1000)
    label_scatter = L"Glauber qHMC $\gamma$ = "*string(gamma)*" t = "*string(t)*" k= "*string(scale)
    plt.scatter(N_values, gap)
    plt.plot(x,exp.(model(x,param)),linestyle = "dashed", label = label_scatter)
end
plt.title(L"qHMC Gap Scaling Comparison $\beta = 300 $")
plt.ylabel(L"$\delta$")
plt.xlabel(L"$N$")
plt.yscale("log")
# plt.xscale("log")
plt.grid("both","both")
plt.legend()
name = "Figures/qHMCScaling/SigmaXgapScalingBeta"*string(beta)*"CompareT"*string(t)*".png"
plt.savefig(name)
plt.show()