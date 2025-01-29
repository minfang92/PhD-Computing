pwd()
using Revise
includet("../src/FHANK.jl")

#----------------------------------------------------------------------------
# BENCHMARK MODEL
#----------------------------------------------------------------------------
ξbar = 0.6
p = Param(;μξ = 0.0,
           σξ = ξbar/sqrt(12),
           σz = 0.05,)
ss = SteadyState(p; w = 1.33)
@show wss = @time find_SteadyState(p,ss;updating = 0.1)

###Save to Plot
tr = DataFrame()
ρ = 0.8

#----------------------------------------------------------------------------
# UP 5
#----------------------------------------------------------------------------
Ashock = 5/100
T = 100
tr1 = Transition(p,ss,T)
Tπ = zeros(T)
for t in 2:T-1
    tr1.TA[t] =  1 + Ashock * ρ^(t-2)
end

### GE Transition
@time find_transition_RBC(p,tr1; tol = 1e-6,update = 0.2)
@show ((tr1.TI .- tr1.TI[1])./tr1.TI[1])[1:5]
@show ((tr1.TY .- tr1.TY[1])./tr1.TY[1])[1:5]

tr.Iup5  = (tr1.TI .- tr1.TI[1])./tr1.TI[1] *100


#----------------------------------------------------------------------------
# DN 5
#----------------------------------------------------------------------------
Ashock = -5/100
T = 100
tr1 = Transition(p,ss,T)
Tπ = zeros(T)
for t in 2:T-1
    tr1.TA[t] =  1 + Ashock * ρ^(t-2)
end

### GE Transition
@time find_transition_RBC(p,tr1; tol = 1e-6,update = 0.2)
@show ((tr1.TI .- tr1.TI[1])./tr1.TI[1])[1:5]
@show ((tr1.TY .- tr1.TY[1])./tr1.TY[1])[1:5]

tr.Idn5  = -(tr1.TI .- tr1.TI[1])./tr1.TI[1] *100

#----------------------------------------------------------------------------
# UP 10
#----------------------------------------------------------------------------
Ashock = 10/100
T = 100
tr1 = Transition(p,ss,T)
Tπ = zeros(T)
for t in 2:T-1
    tr1.TA[t] =  1 + Ashock * ρ^(t-2)
end

### GE Transition
@time find_transition_RBC(p,tr1; tol = 1e-6,update = 0.2)
@show ((tr1.TI .- tr1.TI[1])./tr1.TI[1])[1:5]
@show ((tr1.TY .- tr1.TY[1])./tr1.TY[1])[1:5]

tr.Iup10  = (tr1.TI .- tr1.TI[1])./tr1.TI[1] *100


#----------------------------------------------------------------------------
# DN 10
#----------------------------------------------------------------------------
Ashock = -10/100
T = 100
tr1 = Transition(p,ss,T)
Tπ = zeros(T)
for t in 2:T-1
    tr1.TA[t] =  1 + Ashock * ρ^(t-2)
end

### GE Transition
@time find_transition_RBC(p,tr1; tol = 1e-6,update = 0.2)
@show ((tr1.TI .- tr1.TI[1])./tr1.TI[1])[1:5]
@show ((tr1.TY .- tr1.TY[1])./tr1.TY[1])[1:5]

tr.Idn10  = -(tr1.TI .- tr1.TI[1])./tr1.TI[1] *100

### Save Results
CSV.write("3dynamics/transition_zero_mu.csv", tr)
