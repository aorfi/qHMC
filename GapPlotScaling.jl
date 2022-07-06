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





beta = 6
@. model(x,p) = p[1]*x+p[2]
p0 = [-0.2,-12.0]

gapGl = load_object("Data/GlaubLoc/CmpareGlaubLocOBCgapBeta6")
gapG = load_object("Data/Glaub/CompareOBCGlaubGengapBeta6")
gapPer = load_object("Data/qHMC/alphaEtaParam/alpha1eta1beta6")

name = "Data/qHMC/alphaEtaParam/bestGapBeta"*string(beta)
gap = load_object(name)

N_values = (5:12)
N_values_big = (5:15)

# fit = curve_fit(model,N_values,log.(gap),p0)
# param = fit.param
# println(param)
# scale = -round(param[1]/log(2), digits=4)


# fitG = curve_fit(model,N_values_big,log.(gapG),p0)
# paramG = fitG.param
# println(paramG)
# scaleG = -round(paramG[1]/log(2), digits=4)

fitGl = curve_fit(model,log.(N_values_big),log.(gapGl),p0)
paramGl = fitGl.param
println(paramGl)
scaleGl = round(paramGl[1], digits=4)

fitPer = curve_fit(model,N_values,log.(gapPer),p0)
paramPer = fitPer.param
scalePer = round(paramPer[1]/log(2), digits=4)

label_scatterPer = L"qHMC $\alpha = 1$ $\eta=1$ $2^{-kN}$ k= "*string(scalePer)
plt.scatter(N_values, gapPer, color = "tab:red")
plt.plot(x,exp.(model(x,paramPer )),linestyle = "dashed", label = label_scatterPer ,color = "tab:red")

x = range(5,15, length= 1000)
# label_scatter = L"Glauber qHMC Best $2^{-kN}$ k= "*string(scale)
# plt.scatter(N_values, gap, color = "tab:blue")
# plt.plot(x,exp.(model(x,param)),linestyle = "dashed", label = label_scatter,color = "tab:blue")

# label_scatterG = L"Uniform Glauber $2^{-kN}$ k= "*string(scaleG)
# plt.scatter(N_values_big , gapG, color = "tab:orange")
# plt.plot(x,exp.(model(x,paramG)),linestyle = "dashed", label = label_scatterG,color = "tab:orange")


label_scatterGl = L"Local Glauber $N^{b}$ b= "*string(scaleGl)
plt.scatter(N_values_big , gapGl, color = "tab:green")
plt.plot(x,exp.(model(log.(x),paramGl)),linestyle = "dashed", label = label_scatterGl,color = "tab:green")



plt.title(L"qHMC Gap Scaling Comparison $\beta = $"*string(beta))
plt.ylabel(L"$\delta$")
plt.xlabel(L"$N$")
plt.yscale("log")
# plt.xscale("log")
plt.grid("both","both")
plt.legend()
name = "Figures/qHMCScaling/alphaEtaParam/alpha1eta1beta"*string(beta)*".png"
plt.savefig(name)
plt.show()