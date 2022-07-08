using Arpack
using JLD2
include("qHMC.jl")


# Get gaps for different N
beta = 5
h=0

runs = 100
num_values = 10
alpha_values = range(0,30, length=num_values)
eta_values = range(0,30, length=num_values)

N_values = (5:5)
for N in N_values
    gap_all = zeros(num_values,num_values)
    print("\nworking on N: ", N)
    for alpha_i in (1:num_values)
        for eta_i in (1:num_values)
            print("\nworking on alpha: ", alpha_i)
            print("  eta: ", eta_i)
            gap = zeros(runs)
            for i in (1:runs)
                println(" Working on run = ",i)
                couplings = rand(Normal(0,1),N)
                M = mixing_matrix(N,couplings,h,beta,alpha_values[alpha_i], eta_values[eta_i]) 
                # e,v  = eigen(M)
                # gap_all[alpha_i,eta_i]= abs(1-e[end-1])
                e,v  = eigs(M, nev = 2, which=:LM)
                gap[i] = 1-abs(e[2])
            end
            gap_all[alpha_i,eta_i] = (1/runs)*sum(gap)
        end
    end
    name = "Data/SK/GridSearch/100av"*string(num_values)*"N"*string(N)*"beta"*string(beta)
    save_object(name, gap_all)
end

