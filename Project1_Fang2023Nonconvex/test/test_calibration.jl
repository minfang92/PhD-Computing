pwd()
include("../src/FHANK.jl")


p = Param(;α = 0.25, ν = 0.60, S = 0.33, ξbar = 0.45, ϕ_k = 2.4, σz = 0.05)
ss = SteadyState(p;U=0,w=1.2)
@show wss = @time find_SteadyState(p,ss;howard_on=false)


N = 1000
T = 500
simul=SS_Simulation(p,N,T)


x = stochastic_SS_simulation(p,ss,simul)

simul.m1[1]
simul.m2[1]
simul.m3[1]
simul.m4[1]
simul.m5[1]
simul.m6[1]
simul.m7[1]
