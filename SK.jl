using LinearAlgebra
using SparseArrays
using Arpack
using JLD2
using Random, Distributions
Random.seed!(123)
include("GlauberMixingGeneral.jl")
include("GlauberMixingLocal.jl")
include("MHMixingLocal.jl")
include("MHMixing.jl")


# N = 10
# h=0
# runs = 100

# temp = 10 .^ (range(-2.5,stop=2.5,length=50))
# beta_values = 1 ./ temp
# gap_all_glab = zeros(length(beta_values))
# gap_all_glab_loc = zeros(length(beta_values))
# gap_all_MH = zeros(length(beta_values))
# gap_all_MH_loc = zeros(length(beta_values))
# for j in (1:length(beta_values))
#     beta = beta_values[j]
#     println(" Working on beta = ",beta)
#     gap_glab = zeros(runs)
#     gap_glab_loc = zeros(runs)
#     gap_MH = zeros(runs)
#     gap_MH_loc = zeros(runs)
#     for i in (1:runs)
#         println(" Working on run = ",i)
#         couplings = rand(Normal(0,1),N)
#         M_glab = mixing_glab(N,couplings,h,beta)
#         try
#             e,v  = eigs(M_glab, nev = 2, which=:LM)
#             gap_glab[i] = 1-abs(e[2])
#         catch
#             println("Arpack method out of iteration")
#             e,v = eigen(M_glab)
#             gap_glab[i] = 1-abs(e[end-1])
#         end
#         M_glab_loc = mixing_glab_loc(N,couplings,h,beta)
#         try
#             e,v  = eigs(M_glab_loc, nev = 2, which=:LM)
#             gap_all_glab_loc[j] = 1-abs(e[2])
#         catch
#             println("Arpack method out of iteration")
#             e,v = eigen(M_glab_loc)
#             gap_glab_loc[i] = 1-abs(e[end-1])
#         end
#         M_MH = mixing_MH(N,couplings,h,beta)
#         try
#             e,v  = eigs(M_MH, nev = 2, which=:LM)
#             gap_MH[i] = 1-abs(e[2])
#         catch
#             println("Arpack method out of iteration")
#             e,v = eigen(M_MH)
#             gap_MH[i] = 1-abs(e[end-1])
#         end
#         M_MH_loc = mixing_MH_loc(N,couplings,h,beta)
#         try
#             e,v  = eigs(M_MH_loc , nev = 2, which=:LM)
#             gap_all_MH_loc[j] = 1-abs(e[2])
#         catch
#             println("Arpack method out of iteration")
#             e,v = eigen(M_MH_loc)
#             gap_MH_loc[i] = 1-abs(e[end-1])
#         end
#     end
#     gap_all_glab[j] = (1/runs)*sum(gap_glab)
#     gap_all_glab_loc[j] = (1/runs)*sum(gap_glab_loc)
#     gap_all_MH[j] = (1/runs)*sum(gap_MH)
#     gap_all_MH_loc[j] = (1/runs)*sum(gap_MH_loc)
#     save_object("Data/SK/Classical/avGlaubN10", gap_all_glab)
#     save_object("Data/SK/Classical/avGlaubLocN10", gap_all_glab_loc)
#     save_object("Data/SK/Classical/avMHN10", gap_all_MH)
#     save_object("Data/SK/Classical/avMHLocN10", gap_all_MH_loc)
# end


beta = 5
N_values = (5:12)
h=0
runs = 100
gap_all_glab = zeros(length(N_values))
gap_all_glab_loc = zeros(length(N_values))
gap_all_MH = zeros(length(N_values))
gap_all_MH_loc = zeros(length(N_values))


for j in (1:length(N_values))
    N = N_values[j]
    println(" Working on N = ",N)
    gap_glab = zeros(runs)
    gap_glab_loc = zeros(runs)
    gap_MH = zeros(runs)
    gap_MH_loc = zeros(runs)
    for i in (1:runs)
        println(" Working on run = ",i)
        couplings = rand(Normal(0,1),N)
        M_glab = mixing_glab(N,couplings,h,beta)
        try
            e,v  = eigs(M_glab, nev = 2, which=:LM)
            gap_glab[i] = abs(1-abs(e[2]))
        catch
            println("Arpack method out of iteration")
            e,v = eigen(M_glab)
            gap_glab[i] = abs(1-abs(e[end-1]))
        end
        M_glab_loc = mixing_glab_loc(N,couplings,h,beta)
        try
            e,v  = eigs(M_glab_loc, nev = 2, which=:LM)
            gap_glab_loc[i] = abs(1-abs(e[2]))
        catch
            println("Arpack method out of iteration")
            e,v = eigen(M_glab_loc)
            gap_glab_loc[i] = abs(1-abs(e[end-1]))
        end
        M_MH = mixing_MH(N,couplings,h,beta)
        try
            e,v  = eigs(M_MH, nev = 2, which=:LM)
            gap_MH[i] = abs(1-abs(e[2]))
        catch
            println("Arpack method out of iteration")
            e,v = eigen(M_MH)
            gap_MH[i] = abs(1-abs(e[end-1]))
        end
        M_MH_loc = mixing_MH_loc(N,couplings,h,beta)
        try
            e,v  = eigs(M_MH_loc , nev = 2, which=:LM)
            gap_MH_loc[i] = abs(1-abs(e[2]))
        catch
            println("Arpack method out of iteration")
            e,v = eigen(M_MH_loc)
            gap_MH_loc[i] = abs(1-abs(e[end-1]))
        end
    end
    gap_all_glab[j] = (1/runs)*sum(gap_glab)
    gap_all_glab_loc[j] = (1/runs)*sum(gap_glab_loc)
    gap_all_MH[j] = (1/runs)*sum(gap_MH)
    gap_all_MH_loc[j] = (1/runs)*sum(gap_MH_loc)
    save_object("Data/SK/Classical/av100GlaubBeta5", gap_all_glab)
    save_object("Data/SK/Classical/av100GlaubLocBeta5", gap_all_glab_loc)
    save_object("Data/SK/Classical/av100MHBeta5", gap_all_MH)
    save_object("Data/SK/Classical/av100MHLocBeta5", gap_all_MH_loc)
end
