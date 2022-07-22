using PyPlot
using LaTeXStrings
using JLD2

N_values = (5:8)
beta = 5
num_values = 100
frac = zeros(length(N_values))
for i in (1:length(N_values))
    N = N_values[i]
    name =  "Data/SK/GridSearch/100av"*string(num_values)*"N"*string(N)*"beta"*string(beta)
    gap_all= load_object(name)
    gapG = load_object("Data/SK/Classical/av100GlaubBeta5")[N-4]
    # max,cord = findmax(gap_all)
    # println((cord[2]-1)/num_values*max_param)
    # println((cord[1]-1)/num_values*max_param)
    max_param = 30

    cutoff = zeros(num_values, num_values)
    counter = 0
    all = 0
    for i in (1:num_values)
        for j in (1:num_values)

            all += 1
            if gap_all[i,j]>=gapG
                cutoff[i,j] = 1
                counter += 1
            end

        end
    end

    frac[i] = counter/all

    plt.title(L"qHMC Spectral Gap $\beta= $ "*string(beta)*L" $N = $"*string(N))
    plt.imshow(cutoff, origin="lower",extent = [0, max_param, 0 , max_param])
    # bar = plt.colorbar()
    # plt.scatter((cord[2]-1)/num_values*max_param, (cord[1]-1)/num_values*max_param, color="red", marker=".")
    plt.xlabel(L"$\eta$")
    plt.ylabel(L"$\alpha$")

    # plt.clim(0, 0.05) 
    # bar.set_label(L"$\delta$")
    name = "Figures/SK/GridSearch/Cutoff/"*string(num_values)*"N"*string(N)*"beta"*string(beta)*".png"
    plt.savefig(name)
    plt.clf()
    # plt.show()
end
print(frac)
save_object("Data/SK/qHMC/fracScalingBeta6",frac)



# N_values = (5:8)
# frac = load_object("Data/SK/qHMC/fracScalingBeta6")
# fracT = load_object("Data/SK/qHMC/fracScalingBeta6Triangle")
# plt.scatter(N_values, frac, label = L"$\alpha,\eta\in[0,30]$")
# plt.scatter(N_values, fracT, label = L"$\alpha\in[0,30], \eta<\alpha$")
# plt.title(L"Fraction of $\alpha,\eta >$ Uniform Gap")
# plt.grid()
# plt.xlabel("N")
# plt.ylabel("Fraction")
# plt.ylim(0.7,1)
# plt.legend()
# name = "Figures/SK/GridSearch/Cutoff/Scaling.png"
# plt.savefig(name)
# plt.show()
