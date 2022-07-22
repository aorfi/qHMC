using LinearAlgebra
using SparseArrays
using Arpack
using JLD2
using PyPlot
using OrderedCollections
include("qHMC.jl")

function mixing_matrix_random(N,couplings,h,beta,alpha_start,alpha_end, eta_start, eta_end,num) 
    dim = 2^N
    eta_values = range(eta_start, eta_end, num)
    alpha_values = range(alpha_start, alpha_end, num)
    prob_avg = zeros(dim,dim)
    counter = 0
    for i in (1:num)
        println("working on i = ", i)
        for j in (1:num)
                counter +=1
                H = ham(alpha_values[j], eta_values[i], N,couplings, h)
                U = exp(-1im*H)
                prob = (abs.(U)).^2
                prob_avg += prob
        end
    end
    prob_avg = 1/counter*prob_avg
    M = zeros(dim,dim)
    for col in (0:dim-1)
        E_sp = ising_energy(N, couplings,h , col)
        diag = 0 
        for row in (0:dim-1)
            if row != col
                E_s = ising_energy(N, couplings,h,row)
                # MH
                # m = prob[col+1,row+1]*min(1,exp(-beta*(E_sp-E_s)))
                # Glaub
                m = prob_avg[col+1,row+1]*(1/(1+exp(-beta*(E_s-E_sp))))
                M[col+1,row+1] = m
                diag += m
            end
        M[col+1,col+1] = 1-diag
        end
    end
    return M 
end


beta = 5
N_values = (5:10)
h=0
num = 20
alpha_start = 0
alpha_end = 10
eta_start = 0
eta_end = 10
gap_all = zeros(length(N_values))
runs = 50

for j in (1:length(N_values))
    N = N_values[j]
    println("  N = ",N)
    gap = zeros(runs)
    for i in (1:runs)
        couplings = ones(N)
        couplings[end] = 0 
        M = mixing_matrix_random(N,couplings,h,beta,alpha_start,alpha_end, eta_start, eta_end,num) 
        try
            e,v  = eigs(M, nev = 2, which=:LM)
            gap[i] = 1-abs(e[2])
        catch
            println("Arpack method out of iteration")
            e,v  = eigen(M)
            gap[i] = 1-abs(e[end-1])
        end
        gap_all[j] = (1/runs)*sum(gap)
    end
    save_object("Data/SK/qHMC/20random0-10Beta5", gap_all)
end


