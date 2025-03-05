
pwd()
using Revise
includet("../src/FHANK.jl")

### Load parameters and steadystate
p = Param()
ss = SteadyState(p;U=0,w=1.2669395837406647)
@show wss = @time find_SteadyState(p,ss)

ss.ξstarss
ss.vfss


###Setup and Check Prices
T = 100
tr = Transition(p,ss,T)

tr.TvfA[100][1]



find_transition_RBC(p, tr)



for t in 2:2
    tr.Tq[t] = 0.9
end



@show tr.Tq[1:5]




###Test Backward: Should have reasonable Monotonicities
    #Also, the errors along the transition should be small
@time backward(p, tr, tr.TΛ, tr.Tpw, tr.Tw, tr.TU;Tq = tr.Tq)
@time forward(p,tr,tr.Tpw,tr.Tw,tr.TU)

@show tr.TI[1:5]


















##########################################################################
### Test for NO SHOCK: The stability of Transition which always stay at SS
##########################################################################
###Setup and Check Prices
T = 200
tr = Transition(p,ss,T)
@show tr.Tϵm[1:10]
@show tr.Tπ[1:10]
@show tr.Tr[1:10]
@show tr.TΛ[1:10]
@show tr.Tpw[1:10]
@show tr.Tw[1:10]

###Test Terminal Conditions
#terminal Expected Value Function should be same [by defination of Transition]
sum(abs, tr.TEV[T][1] - ss.EVss)
#one time backward update should be the same [if expectation function works]
tr.TEV[T-1][1] = expectation(p,tr.Tvf[T][1],p.Πz)
sum(abs, tr.TEV[T-1][1] - tr.TEV[T][1])
#Bellman does not converge ???
tr.Tvf[T-1][1], tr.Tkpolicy[T-1][1]  = bellman(p,tr.Tpw[T-1],tr.Tw[T-1],tr.TΛ[T-1],tr.TEV[T-1][1],tr.Tprofit[T-1][1])
sum(abs, tr.Tvf[T-1][1] - tr.Tvf[T][1])






### Test Aggregation: Should have reasonable Monotonicities
    #Also, the errors along the transition should be small
for t in 2:T-1
    tr.TC[t] =  tr.TY[t] - tr.TI[t] - tr.TAC_K[t]  - p.ψ/2*tr.Tπ[t]^2*tr.TY[t]
    tr.Tw[t] =  p.θ*tr.TC[t]^p.η*tr.TN[t]^(p.χ-1)
end
@show tr.TC
@show tr.Tw









##########################################################################
### Test for Small MP SHOCK: The stability of Transition not far from SS
##########################################################################
### Test LBJ method [One-Time small shock]
# Create Initial Transition: Monetary Policy Shock Only and SS
T = 200
tr = Transition(p,ss,T)
#Get shock and an initial guess
     # !!!! Shock comes at t=2 !!!!
mshock = -0.0624/100
for t in 2:2
    tr.Tϵm[t] = mshock * 0.5^(t-2)
end

#------------ First Round -------------#

### Step 1: Run LBJ [Recursive Block]
@time LBJ(p,tr;maxit=100, Dampening=0.25, DampeningThresh=1e-4)
# Results
    # By Calculation on Paper, this One-Time Shock should almost only change
    # 1. π[2] = -ϵm[2]/ϕ_π = 0.0624/1.25 = 0.04992%
    # 2. pw[2] = - ψ/(γ-1)ϕ_π * ϵm[2] = 8 * 0.0624 = 0.4992%
@show tr.TC[1:10]
@show tr.Tr[1:10]
@show tr.TΛ[1:10]
@show tr.Tw[1:10]
@show tr.Tπ[1:10]
@show tr.Tpw[1:10]

### Step 2: Run Backward [HA Block]
@time backward(p, tr, tr.TΛ, tr.Tpw, tr.Tw, tr.TS)
#Changes in HA Block
    # Labor supply policy go up
    # Investment policy go up
    # VF go up
gap_lpolicy = zeros(T)
for t in 1:T
    gap_lpolicy[t] = maximum(abs,tr.Tlpolicy[1][1]-tr.Tlpolicy[t][1])
end
@show gap_lpolicy[1:10]
# Monotonicity of Tprofit
gap_profit = zeros(T)
for t in 1:T
    gap_profit[t] = maximum(abs,tr.Tprofit[1][1]-tr.Tprofit[t][1])
end
@show gap_profit[1:10]
# Monotonicity of TEV
gap_EV= zeros(T)
for t in 1:T
    gap_EV[t] = maximum(abs,tr.TEV[1][1]-tr.TEV[t][1])
end
@show gap_EV[1:10]
# Monotonicity of Tvf
gap_vf = zeros(T)
for t in 1:T
    gap_vf[t] = maximum(abs,tr.Tvf[1][1]-tr.Tvf[t][1])
end
@show gap_vf[1:10]
# Monotonicity of Tkpolicy
gap_kpolicy = zeros(T)
for t in 1:T
    gap_kpolicy[t] = maximum(abs,tr.Tkpolicy[1][1]-tr.Tkpolicy[t][1])
end
@show gap_kpolicy[1:10]

### Step 3: Run Forward[HA Block]
@time forward(p,tr,tr.Tpw,tr.Tw, tr.TS)
#Changes in HA Block
    # Labor supply policy go up
    # Investment policy go up
    # VF go up
gap_Y = zeros(T)
for t in 1:T
    gap_Y[t] = maximum(abs,tr.TY[200]-tr.TY[t])
end
@show gap_Y[1:20]
# Monotonicity of TN
gap_N = zeros(T)
for t in 1:T
    gap_N[t] = maximum(abs,tr.TN[200]-tr.TN[t])
end
@show gap_N[1:20]
# Monotonicity of TK
gap_K = zeros(T)
for t in 1:T
    gap_K[t] = maximum(abs,tr.TK[200]-tr.TK[t])
end
@show gap_K[1:20]
# Monotonicity of TI
gap_I = zeros(T)
for t in 1:T
    gap_I[t] = maximum(abs,tr.TI[200]-tr.TI[t])
end
@show gap_I[1:20]


#------------ Second Round -------------#
    # The Nominal Interest is updated at wrong direction !!!
    # how to fix this?

### Step 1: Run LBJ [Recursive Block]
@time LBJ(p,tr;maxit=10000, Dampening=0.25, DampeningThresh=1e-4)
# Results
    # By Calculation on Paper, this One-Time Shock should almost only change
    # 1. π[2] = -ϵm[2]/ϕ_π = 0.0624/1.25 = 0.04992%
    # 2. pw[2] = - ψ/(γ-1)ϕ_π * ϵm[2] = 8 * 0.0624 = 0.4992%
@show tr.TC[1:10]
@show tr.Tr[1:10]
@show tr.TΛ[1:10]
@show tr.Tw[1:10]
@show tr.Tπ[1:10]
@show tr.Tpw[1:10]

### Step 2: Run Backward [HA Block]
@time backward(p, tr, tr.TΛ, tr.Tpw, tr.Tw, tr.TS)
#Changes in HA Block
    # Labor supply policy go up
    # Investment policy go up
    # VF go up
gap_lpolicy = zeros(T)
for t in 1:T
    gap_lpolicy[t] = maximum(abs,tr.Tlpolicy[1][1]-tr.Tlpolicy[t][1])
end
@show gap_lpolicy[1:10]
# Monotonicity of Tprofit
gap_profit = zeros(T)
for t in 1:T
    gap_profit[t] = maximum(abs,tr.Tprofit[1][1]-tr.Tprofit[t][1])
end
@show gap_profit[1:10]
# Monotonicity of TEV
gap_EV= zeros(T)
for t in 1:T
    gap_EV[t] = maximum(abs,tr.TEV[1][1]-tr.TEV[t][1])
end
@show gap_EV[1:10]
# Monotonicity of Tvf
gap_vf = zeros(T)
for t in 1:T
    gap_vf[t] = maximum(abs,tr.Tvf[1][1]-tr.Tvf[t][1])
end
@show gap_vf[1:10]
# Monotonicity of Tkpolicy
gap_kpolicy = zeros(T)
for t in 1:T
    gap_kpolicy[t] = maximum(abs,tr.Tkpolicy[1][1]-tr.Tkpolicy[t][1])
end
@show gap_kpolicy[1:10]

### Step 3: Run Forward[HA Block]
@time forward(p,tr,tr.Tpw,tr.Tw, tr.TS)
#Changes in HA Block
    # Labor supply policy go up
    # Investment policy go up
    # VF go up
gap_Y = zeros(T)
for t in 1:T
    gap_Y[t] = maximum(abs,tr.TY[200]-tr.TY[t])
end
@show gap_Y[1:20]
# Monotonicity of TN
gap_N = zeros(T)
for t in 1:T
    gap_N[t] = maximum(abs,tr.TN[200]-tr.TN[t])
end
@show gap_N[1:20]
# Monotonicity of TK
gap_K = zeros(T)
for t in 1:T
    gap_K[t] = maximum(abs,tr.TK[200]-tr.TK[t])
end
@show gap_K[1:20]
# Monotonicity of TI
gap_I = zeros(T)
for t in 1:T
    gap_I[t] = maximum(abs,tr.TI[200]-tr.TI[t])
end
@show gap_I[1:20]
