using Revise
includet("../src/FHANK.jl") # ".." to go upper layer directory




N = 11
ρ = 0.95
σ = 0.05
μ = -σ^2/2

sigma = [1,2,3,4]

sigma[1:2] == [1,2]


x = rouwenhorst(N, ρ, σ, 0)


collect(x.state_values)

x.p


y = tauchen(11, 0.95, 0.05)

z = tauchen(11, 0.95, 0.2)

y.state_values

z.state_values


x = stochastic_volatility_tauchen(11, 0.95, [0.05,0.2])






x[1].state_values

x[2].state_values



x[1].p
y.p

x[2].p
z.p

TR = zeros(11,11)

TR[2,1]

zgrid = collect(stochastic_volatility_tauchen(Nz, ρz, [0.05,0.2])[1].state_values)

zgrid[1]

zgrid[11]

TR = tail_risk(N, 0.02, 0, zgrid)

TR + x[1].p

x[1].p
