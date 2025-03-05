using Revise
includet("../src/FHANK.jl")

#Step1: Initialize steady state

pw = 1.0
w = 1.4
Λ = 0.99
U = 0
q = 1.0
p = Param()
Firm_VFI(p, U, pw, w, Λ, q)
p.kpolicyA

p.Dist_kz

x = Young(p,U,p.kpolicyA,p.kpolicyNA,p.ξstar,p.Dist_kz)

x = Young(p,U,p.kpolicyA,p.kpolicyNA,p.ξstar,x)


#Step3: Test find stationary equilibrium
find_stationary_distributions(p,U)
p.Dist_kz


#Step4: Test aggregation over distribution
@time N, Y, I, IA, INA, K, AC_K = calculate_aggregation(p,U,p.Dist_kz,p.ξstar,p.lpolicy,p.kpolicyA,p.kpolicyNA)
