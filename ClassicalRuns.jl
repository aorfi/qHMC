using LinearAlgebra
using SparseArrays
using Arpack
using JLD2
using Random, Distributions
include("GlauberMixingGeneral.jl")
include("GlauberMixingLocal.jl")
include("MHMixingLocal.jl")
include("MHMixing.jl")


N = 10
h=1.5
couplings = ones(N)
# couplings[end] = 0
temp = 10 .^ (range(-2.5,stop=2.5,length=50))
beta_values = 1 ./ temp
gap_glab = zeros(length(beta_values))
gap_glab_loc = zeros(length(beta_values))
gap_MH = zeros(length(beta_values))
gap_MH_loc = zeros(length(beta_values))
for j in (1:length(beta_values))
    beta = beta_values[j]
    println(" Working on beta = ",beta)

    M_glab = mixing_glab(N,couplings,h,beta)
    try
        e,v  = eigs(M_glab, nev = 2, which=:LM)
        gap_glab[j] = abs(1-abs(e[2]))
    catch
        println("Arpack method out of iteration")
        e,v = eigen(M_glab)
        gap_glab[j] = abs(1-abs(e[end-1]))
    end
    M_glab_loc = mixing_glab_loc(N,couplings,h,beta)
    try
        e,v  = eigs(M_glab_loc, nev = 2, which=:LM)
        gap_glab_loc[j] = abs(1-abs(e[2]))
    catch
        println("Arpack method out of iteration")
        e,v = eigen(M_glab_loc)
        gap_glab_loc[j] = abs(1-abs(e[end-1]))
    end
    M_MH = mixing_MH(N,couplings,h,beta)
    try
        e,v  = eigs(M_MH, nev = 2, which=:LM)
        gap_MH[j] = abs(1-abs(e[2]))
    catch
        println("Arpack method out of iteration")
        e,v = eigen(M_MH)
        gap_MH[j] = abs(1-abs(e[end-1]))
    end
    M_MH_loc = mixing_MH_loc(N,couplings,h,beta)
    try
        e,v  = eigs(M_MH_loc , nev = 2, which=:LM)
        gap_MH_loc[j] = abs(1-abs(e[2]))
    catch
        println("Arpack method out of iteration")
        e,v = eigen(M_MH_loc)
        gap_MH_loc[j] = abs(1-abs(e[end-1]))
    end
end
save_object("Data/Glaub/PBCN10h"*string(h), gap_glab)
save_object("Data/GlaubLoc/PBCN10h"*string(h), gap_glab_loc)
save_object("Data/MH/PBCN10h"*string(h), gap_MH)
save_object("Data/MHLoc/PBCN10h"*string(h), gap_MH_loc)
println("SAVED")