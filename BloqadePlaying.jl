using Bloqade
using PythonCall # Use matplotlib to generate plots
matplotlib = pyimport("matplotlib")
plt = pyimport("matplotlib.pyplot")

nsites = 3
Ω = 2π * 4
C6 = 2π * 862690
Rb = (C6/Ω)^(1/6)
a = Rb/sqrt(2)

atoms = generate_sites(ChainLattice(),nsites,scale=a) 
h = rydberg_h(atoms; Ω = Ω, Δ = C6/(2*a^6))
# H = mat(h)
# display(H)
# init_state = zero_state(nsites)



# prob = SchrodingerProblem(init_state, 2π, h)
# integrator = init(prob, Vern8())
# emulate!(prob)

# rydberg_populations = map(1:nsites) do i
#     rydberg_density(prob.init_state, i)
# end

# Tmax = 6.0
# nsteps = 2001
# times = LinRange(0, Tmax, nsteps)
# dt = Tmax / (nsteps - 1)

# prob = SchrodingerProblem(init_state, Tmax, h, dt = dt, adaptive = false);
# integrator = init(prob, Vern8());

# densities = [] # Time evolve the system in the full space
# for _ in TimeChoiceIterator(integrator, 0.0:dt:Tmax)
#     push!(densities, rydberg_density(init_state, 1))
# end

# ax = plt.subplot(1, 1, 1)
# plt.plot(times, real(densities), "k", label = "Full space")
# plt.xlabel("Time (us)")
# plt.ylabel("Rydberg density")
# plt.tight_layout()
# plt.legend()
# plt.show()