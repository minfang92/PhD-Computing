pwd()
using Revise
includet("../src/FHANK.jl")

#####################################################################################################
# Section 1: Solve steadystates under different ACs [Moments]
#####################################################################################################
"""
Average investment rate                10.4%
Standard Deviation                     0.16
Inaction rate                          23.7%
Spike rate                             14.4%
Negative rate                          10% ??
Serial correlation of investment rates 0.40
PE elasticity to real interest rate    7
"""

#####################################################################################################
# Section 2: Let's Run the Tests
#####################################################################################################

#----------------------------------------------------------------------------
#Test 1: SteadyState Moments
#----------------------------------------------------------------------------
ξbar = 0.6
p = Param(;α = 0.25, ν = 0.60, S = 0.0, μξ = ξbar/2, σξ = ξbar/sqrt(12), ϕ_k = 0.0, a = 0.001, σz = 0.05)
ss = SteadyState(p;U=0,w=1.2878046517828643)
@show wss = @time find_SteadyState(p,ss;updating = 0.2)

N = 50000
T = 750

simul = SS_Simulation(p,N,T)
stochastic_SS_simulation(p,ss,simul)

simul.m1[1]
simul.m2[1]
simul.m3[1]
simul.m4[1]
simul.m5[1]
simul.m6[1]
simul.m7[1]



#----------------------------------------------------------------------------
#Test 2: Transitional Moments
#----------------------------------------------------------------------------
T = 200

#PE transition
tr = Transition(p,ss,T)
Rshock = 0.0025
for t in 2:2
    tr.TΛ[t] =  tr.TΛ[1] + Rshock * 0.5^(t-2)
end

backward(p, tr, tr.TΛ, tr.Tpw, tr.Tw, tr.TU)
forward(p,tr,tr.Tpw,tr.Tw,tr.TU)

RI =  (tr.TI[2] - tr.TI[1])/tr.TI[1]
@show els = -RI/Rshock
