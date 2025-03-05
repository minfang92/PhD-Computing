using Revise
includet("../src/FHANK.jl")

p = Param()

p.kgrid


u = Uncertainty(p,0)
u.Δz[1]
u = Uncertainty(p,1)
u.Δz[1]
u = Uncertainty(p,2)
u.Δz[1]
u = Uncertainty(p,3)
u.Δz[1]


p.zgrid


u.Δz

u.Πz

p.T

p.mc_z
