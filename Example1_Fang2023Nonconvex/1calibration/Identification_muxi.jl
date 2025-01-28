pwd()
using Revise
includet("../src/FHANK.jl")


#####################################################################################################
# Section 1: Auto Correlation to identify ξbar
#####################################################################################################
!!! Don't Rerun this part [will erease all previous outputs]
# First Time create Dataframe to hold results
!!!
Identification_df = DataFrame()
!!!
CSV.write("1calibration/02.Identification_results/Identification_df_Auto2.csv", Identification_df)

# Firstload dataframe
Identification_df = CSV.read("1calibration/02.Identification_results/Identification_df_Auto2.csv", copycols = true, typemap = Dict(Int64=>Float64))


# Function
function Identification_df_Auto2(α = 0.25,ν = 0.60; #KEY
                          ϕk_lower = 0.0,
                          ϕk_upper = 2.4,
                          muxi_lower = 0.0,
                          muxi_upper = 0.65,
                          ξ_number = 40,
                          ξbar = 0.65,
                          T = 50,
                          Rshock = 0.0025,
                          w = 1.2)

        for iϕk in 0:1
            for iξ in 1:ξ_number

                    #Params
                    ϕ_k = ϕk_lower + ϕk_upper*iϕk
                    muxi = muxi_lower + (muxi_upper-muxi_lower)/ξ_number*iξ
                    #Params
                    p = Param(α = α, ν = ν; μξ = muxi, σξ = ξbar/sqrt(12), ϕ_k = ϕ_k)
                    ss = SteadyState(p; w = w)

                    #Solve
                    @show wss = @time find_SteadyState(p,ss;updating = 0.2)
                    w = wss
                    tr = Transition(p,ss,T)

                    #PE transition
                    for t in 2:2
                        tr.TΛ[t] =  tr.TΛ[1] + Rshock * 0.5^(t-2)
                    end

                    backward(p, tr, tr.TΛ, tr.Tpw, tr.Tw, tr.TU)
                    forward(p,tr,tr.Tpw,tr.Tw,tr.TU)

                    RI =  (tr.TI[2] - tr.TI[1])/tr.TI[1]
                    @show els = RI/Rshock

                    #Update Parameters and Moments
                    results = (alpha=p.α, nu=p.ν, phi_k=p.ϕ_k, muxi=muxi, a=p.a, wss=wss, els=els)
                    push!(Identification_df,results)
                    CSV.write("1calibration/02.Identification_results/Identification_df_Auto2.csv", Identification_df)

            end
        end

end

Identification_df_Auto2()
