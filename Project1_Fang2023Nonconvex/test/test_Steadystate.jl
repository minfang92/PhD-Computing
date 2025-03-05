using Revise
includet("../src/FHANK.jl")
pwd()

#Step1: Initialize steady state
p = Param()
ss = SteadyState(p;w=1.0)

ss.U
ss.q
@show ss


wss = @time find_SteadyState(p,ss)

@show ss.Qss




p.zgrid[2,1]

p.kpolicy[:,4] - p.kgrid

p.ξstar




#Step2: Distributional Targets
p.vf -ss.vfss
p.kpolicy - ss.kpolicyss

sum(ss.Dist_kz, dims=1)
sum(ss.Dist_kz, dims=2)
sum(sum(ss.Dist_kz, dims=2))
plot(p.kgrid, p.kpolicy[:,9])

#Step3: Aggregate Targets
@show ss.Pss
@show ss.Qss

ss.vfss - p.vf
ss.EVss - expectation(p,ss.vfss,p.Πz)


EV_func = expectation(p,ss.vfss,p.Πz)
vf, kpolicy = bellman(p, ss.Pss[4], ss.Pss[5], ss.Pss[3], EV_func, ss.profitss)
vf - ss.vfss #this does NOT LOOK RIGHT !!!
