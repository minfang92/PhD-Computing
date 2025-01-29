pwd()
using Revise
includet("../src/FHANK.jl")

#####################################################################################################
# Section One: Solve steadystates
#####################################################################################################
ξbar = 0.6
p = Param(;β = 0.99,
           α = 0.25,
           ν = 0.60,
           S = 0.00,
           a = 0.001,
           ϕ_k = 0.00,
           σz = 0.05,
           μξ = ξbar/2,
           σξ = ξbar/sqrt(12)
           )

#-------------- Low uncertianty SS --------------#
ss = SteadyState(p; U = 0, w = 1.33)
@show wss = @time find_SteadyState(p,ss)


#####################################################################################################
# Section Two: Convert Data for R plot
#####################################################################################################

#-------------- [Plot Investment Rate] --------------#
Δk = zeros(p.Nk,p.Nz)
for iz in 1:p.Nz, ik in 1:p.Nk
    Δk[ik,iz] = (ss.kpolicyAss[ik,iz] - (1-p.δ)*p.kgrid[ik])/p.kgrid[ik]
end
Δk #Investment Rate
knots = (p.kgrid, p.zgrid)
Δk_itp = interpolate(knots, Δk, Gridded(Linear()))

kplot = collect(range(p.kgrid[1], stop=p.kgrid[p.Nk], length=50))
zplot = collect(range(p.zgrid[1], stop=p.zgrid[p.Nz], length=50))
Δk_plot = zeros(length(kplot),length(zplot))
for ik in 1:length(kplot)
    for iz in 1:length(zplot)
        Δk_plot[ik,iz] = Δk_itp(kplot[ik],zplot[iz])
    end
end
Δk_plot
data = DataFrame(Δk_plot)
CSV.write("2steadystate/InvestmentRate_bundled.csv", data)


#-------------- [Plot Adjustment Prob.] --------------#
Prob = zeros(p.Nk,p.Nz)
for iz in 1:p.Nz, ik in 1:p.Nk
    Prob[ik,iz] = (ss.ξstarss[ik,iz].-p.ξlowbar)./(p.ξupbar-p.ξlowbar)
end
Prob
knots = (p.kgrid, p.zgrid)
Prob_itp = interpolate(knots, Prob, Gridded(Linear()))

kplot = collect(range(p.kgrid[1], stop=p.kgrid[p.Nk], length=50))
zplot = collect(range(p.zgrid[1], stop=p.zgrid[p.Nz], length=50))
Prob_plot = zeros(length(kplot),length(zplot))
for ik in 1:length(kplot)
    for iz in 1:length(zplot)
        Prob_plot[ik,iz] = Prob_itp(kplot[ik],zplot[iz])
    end
end
Prob_plot
data = DataFrame(Prob_plot)
CSV.write("2steadystate/Prob_of_Adjustment_bundled.csv", data)



#####################################################################################################
# Section Three: Distribution of the Response
#####################################################################################################
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

###Save to Plot
tr_output = DataFrame()
tr.TI

#Investment Counterfactural
TIex = deepcopy(tr.TI)
TIin = deepcopy(tr.TI)
for t in 2:T-1
    #Intensive
    N, Y, I1, IA, INA, K, AC_k = calculate_aggregation(p,0,tr.TDist_kz[t][1],ss.ξstarss,tr.Tlpolicy[t][1],tr.TkpolicyA[t][1],tr.TkpolicyNA[t][1])
    TIin[t] = I1
    #Extensive
    N, Y, I2, IA, INA, K, AC_k = calculate_aggregation(p,0,tr.TDist_kz[t][1],tr.Tξstar[t][1],tr.Tlpolicy[t][1],ss.kpolicyAss,ss.kpolicyNAss)
    TIex[t] = I2
end
tr_output.Iex = (TIex .- TIex[1])./TIex[1] *100
tr_output.Iin = (TIin .- TIin[1])./TIin[1] *100
tr_output.I  = (tr.TI .- tr.TI[1])./tr.TI[1] *100

### Save Results
CSV.write("2steadystate/transition_bundled.csv", tr_output)
# Readback
tr_output = CSV.read("2steadystate/transition_bundled.csv",  copycols = true, typemap = Dict(Int64=>Float64))
