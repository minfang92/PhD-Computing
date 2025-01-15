#This File include several functions to calibrate the SteadyState
#files including:
#  Stochastic Simulation
#


#This File include several functions to solve the stedaystate
"""
Generate Investment distribution Moments
------------------------------------------------------
##### Arguments
- `Yi` : Distribution of SteadyState Investment (Annually)
##### Output
- `Mom` : A Vector of Moments
##### Moments [M1,M2,M3,M4,M5]
------------------------------------------------------
PAPER ONE
Investment Distribution                 & LRD
Investment Inaction Rate                & 8.1%
Fraction of Negative Investment         & 10.4%
Spike Rate: Positive Investment         & 18.6%
Spike Rate: Negative Investment         & 1.8%
Serial correlation of Investment Rates  & 0.058
------------------------------------------------------
PAPER TWO [Prefered]
Average investment rate                10.4%
Standard Error                         0.16
Skewness                                3.6
Inaction rate                          23.7%
Spike rate                             14.4%
Serial correlation of investment rates 0.40
Aggregate investment rate              6.9%  (not matching)
Spike share of aggregate investment    24.4% (not matching)
------------------------------------------------------
"""
function SS_Simulation(p::NamedTuple,N::Int,T::Int;
                       Tdrop::Int=250)

    #Periods
    T_simul = T + Tdrop
    T = T
    N = N

    #Firm Sample
    Ti = zeros(T,N) #investment
    Tk = zeros(T,N) #capital
    Tz_value = zeros(T,N) #productivity grid
    Taction = zeros(Int,T,N) #action
    Tage = zeros(T,N) #age sign last active adjustment
    Tlogkz = zeros(T,N) # log(k/z) Named: capital-productivity-ratio
    Tlogkzgap = zeros(T,N) # log(k/z) - mean(log(k/z)) Named: capital-productivity-ratio gap

    #Annualized Moments
    m1 = [mean(Ti),]
    m2 = [std(Ti),]
    m3 = [count(i->(-0.01<=i<=0.01), Ti)/(N*T),]
    m4 = [count(i->(0.2<=i), Ti)/(N*T),]
    m5 = [count(i->(i<0), Ti)/(N*T),]
    m6 = [mean(autocor(Ti,[1])),]
    m7 = [0.0,] #cov(Tlogkzgap,Tage)
    m8 = [skewness(Ti),]
    m9 = [mean(Tlogkz),]

    ###Output
    simul = (N=N,T=T,T_simul=T_simul,
             Ti=Ti,Tk=Tk,Tz_value=Tz_value,Tage=Tage,Taction=Taction,Tlogkz=Tlogkz,Tlogkzgap=Tlogkzgap,
             m1=m1,m2=m2,m3=m3,m4=m4,m5=m5,m6=m6,m7=m7,m8=m8,m9=m9)

end


"""
Stochastic Simulation to generate investment distribution
------------------------------------------------------
##### Arguments
- `p` : Parameters
- `ss`: SteadyState
- `N` : number of firms
- `T`: Years of simulation
- `s`: Uncertainty State
------------------------------------------------------
"""
function stochastic_SS_simulation(p::NamedTuple,ss::NamedTuple,simul::NamedTuple;
                               Q::Int = 4*simul.T_simul, #quarters
                               T_simul::Int = simul.T_simul,
                               T::Int = simul.T,
                               N::Int = simul.N,
                               Ti::Array = zeros(Float64,Q,N),
                               Yi::Array = zeros(Float64,T,N),
                               Tz::Array = zeros(Int,Q,N),
                               Tz_value::Array = zeros(Q,N),
                               Tk::Array=ones(Float64,Q,N),
                               Tξ::Array=ones(Float64,Q,N))

    ###Uncertainty Process
    U = ss.U
    u = Uncertainty(p,U)

    ###Generate Interpolations
    @time PLFA = [interpolate(p.kgrid, ss.kpolicyAss[:,iz],FritschButlandMonotonicInterpolation()) for iz in 1:p.Nz]
    @time PLFNA = [interpolate(p.kgrid, ss.kpolicyNAss[:,iz],FritschButlandMonotonicInterpolation()) for iz in 1:p.Nz]
    @time PLFξ = [interpolate(p.kgrid, ss.ξstarss[:,iz],FritschButlandMonotonicInterpolation()) for iz in 1:p.Nz]

    ###simulate shocks
    @time for i in 1:1:N
        Tz[:,i] = simulate(u.mc_z,Q) #simulate from QuantEcon
    end
    @time for i in 1:1:N, t in 1:1:Q
        Tξ[t,i] = rand()*(p.ξupbar - p.ξlowbar) + p.ξlowbar
    end

    ###Zgrid values
    for i in 1:1:N
        for t in 1:1:Q
            Tz_value[t,i]=exp(p.zgrid[Tz[t,i]])
        end
    end

    ###simulate capital
    # @time Threads.@threads for i in 1:1:N
    @time for i in 1:1:N
        Tk[1,i] = p.kgrid[Int(ceil(p.Nk/2))]
        for t in 2:Q
            if Tξ[t-1,i] > PLFξ[Tz[t-1,i]](Tk[t-1,i])
                Tk[t,i] =  max(p.kgrid[1],PLFNA[Tz[t,i]](Tk[t-1,i]))
            else
                Tk[t,i] = PLFA[Tz[t,i]](Tk[t-1,i])
            end
        end
    end

    ###calculate investment RATE
    @time for i in 1:N
        for t in 2:Q
            Ti[t-1,i] = (Tk[t,i] - (1-p.δ)*Tk[t-1,i]) / Tk[t-1,i]
        end
    end

    #---------Annually-------------
    Tdrop = T_simul - T

    ###Annualized Distribution Path [drop first 100 years]
    @time for i in 1:N
        for t in 1:T
            q = (t-1+Tdrop)*4 #drop first 100 years
            simul.Ti[t,i] =  Ti[q+1,i] + Ti[q+2,i] + Ti[q+3,i] + Ti[q+4,i]
            if abs(simul.Ti[t,i])>0.01
                simul.Taction[t,i] = 1
            end
            simul.Tk[t,i] =  Tk[q+4,i]
            simul.Tz_value[t,i] = Tz_value[q+4,i]
            simul.Tlogkz[t,i] = log(simul.Tk[t,i]/simul.Tz_value[t,i])
        end
    end

    @time for i in 1:N
        if simul.Taction[1,i]==0
            simul.Tage[1,i] = 1
        else
            simul.Tage[1,i] = 0
        end
        for t in 2:T
            if simul.Taction[t,i]==0
                simul.Tage[t,i] = simul.Tage[t-1,i] + 1
            else
                simul.Tage[t,i] = 0
            end
        end
    end

    #---------Annuallized Moments-------------
    simul.m1[1,] = mean(simul.Ti)
    simul.m2[1,] = std(simul.Ti)
    simul.m3[1,] = count(i->(-0.01<=i<=0.01), simul.Ti)/(N*T)
    simul.m4[1,] = count(i->(0.2<=i), simul.Ti)/(N*T)
    simul.m5[1,] = count(i->(i<0), simul.Ti)/(N*T)
    simul.m6[1,] = mean(StatsBase.autocor(simul.Ti[T-110:T,:],[1]))
    simul.m8[1,] = skewness(simul.Ti)

    #--mean logkz should be weighted by size
    K = sum(simul.Tk)
    simul.m9[1,] = mean(simul.Tlogkz.*simul.Tk./K) # E[Tlogkz]

    @time for i in 1:N
        for t in 1:T
            simul.Tlogkzgap[t,i] = simul.Tlogkz[t,i] - simul.m9[1,]
        end
    end

    cov_gap_age = zeros(N)
    for i in 1:N
        cov_gap_age[i] = cov(simul.Tlogkzgap[:,i],simul.Tage[:,i])
    end
    simul.m7[1,] = mean(cov_gap_age)

end


"""
Stochastic Simulation to generate investment distribution over Transition
"""
function TR_Simulation(p::NamedTuple,ss::NamedTuple,tr::NamedTuple, N::Int, Tss::Int; Ttt::Int = Tss + tr.T) #quarterly

    #scalers
    Ttr = tr.T
    TIQRsg = zeros(Ttt)
    TIQRsr = zeros(Ttt)
    Tipo = zeros(Ttt)

    #Firm Sample
    Ti = zeros(Ttt,N) #investment
    Tirate = zeros(Ttt,N) # investment rate
    Tk = zeros(Ttt,N) #capital
    Tac = zeros(Ttt,N) #all capital adjustment costs
    Tξ = zeros(Ttt,N)
    Tl = zeros(Ttt,N) #labor
    Ty = zeros(Ttt,N) #output
    Td = zeros(Ttt,N)
    Tvf = zeros(Ttt,N)
    Tz = ones(Int,Ttt,N)
    Tz_value = zeros(Ttt,N) #productivity grid
    Tstockreturn = zeros(Ttt,N)
    Tsalesgrowth = zeros(Ttt,N)
    Tpub = zeros(Ttt,N)
    Tpubage = zeros(Ttt,N)
    Tpub25plus = zeros(Ttt,N)

    ###Output
    simul = (N=N,Ttt=Ttt,Tss=Tss,Ttr=Ttr,TIQRsg=TIQRsg,TIQRsr=TIQRsr,Tipo=Tipo,
             Ti=Ti,Tirate=Tirate,Tac=Tac,
             Tk=Tk,Tξ=Tξ,Tl=Tl,Ty=Ty,Tz=Tz,Tz_value=Tz_value,
             Td=Td,Tvf=Tvf,Tstockreturn=Tstockreturn,
             Tsalesgrowth=Tsalesgrowth,Tpub=Tpub,Tpubage=Tpubage,Tpub25plus=Tpub25plus)

end


"""
Stochastic Simulation to generate investment distribution over Transition
"""
function stochastic_TR_simulation(p::NamedTuple,ss::NamedTuple,tr::NamedTuple,simul::NamedTuple;
                               Ttt::Int = simul.Ttt,
                               Tss::Int = simul.Tss,
                               Ttr::Int = simul.Ttr,
                               N::Int = simul.N)

    ####---------------------------------------------------------------------------------------------------------------
    #Part2: SteadyState Simulation [no shocks] [T: 1--->Tss+1]

    ###Uncertainty Process
    Uss = ss.U
    uss = Uncertainty(p,Uss)

    ###Generate Interpolations
    PLFA = [interpolate(p.kgrid, ss.kpolicyAss[:,iz],FritschButlandMonotonicInterpolation()) for iz in 1:p.Nz]
    PLFNA = [interpolate(p.kgrid, ss.kpolicyNAss[:,iz],FritschButlandMonotonicInterpolation()) for iz in 1:p.Nz]
    PLFξ = [interpolate(p.kgrid, ss.ξstarss[:,iz],FritschButlandMonotonicInterpolation()) for iz in 1:p.Nz]
    PLFL = [interpolate(p.kgrid, ss.lpolicyss[:,iz], FritschButlandMonotonicInterpolation()) for iz in 1:p.Nz]
    VF = [interpolate(p.kgrid, ss.vfss[:,iz], FritschButlandMonotonicInterpolation()) for iz in 1:p.Nz]
    PLFprofit = [interpolate(p.kgrid, ss.profitss[:,iz], FritschButlandMonotonicInterpolation()) for iz in 1:p.Nz]

    ###simulate shocks
    for i in 1:1:N
        simul.Tz[1:Tss+1,i] = simulate(uss.mc_z,Tss+1) #simulate from QuantEcon
    end
    for i in 1:1:N, t in 1:1:Tss+1
        simul.Tξ[t,i] = rand()*(p.ξupbar - p.ξlowbar) + p.ξlowbar
    end

    ###Zgrid values
    for i in 1:1:N
        for t in 1:1:Tss+1
            simul.Tz_value[t,i]=exp(p.zgrid[simul.Tz[t,i]]+uss.Δz[simul.Tz[t,i]])
        end
    end

    ###simulate capital
    for i in 1:1:N
        simul.Tk[1,i] = p.kgrid[Int(ceil(p.Nk/2))]
        for t in 2:Tss+1
            if simul.Tξ[t-1,i] > PLFξ[simul.Tz[t-1,i]](simul.Tk[t-1,i])
                simul.Tk[t,i] =  max(p.kgrid[1],PLFNA[simul.Tz[t,i]](simul.Tk[t-1,i]))
            else
                simul.Tk[t,i] = PLFA[simul.Tz[t,i]](simul.Tk[t-1,i])
            end
            simul.Ti[t,i] = simul.Tk[t,i] - simul.Tk[t-1,i]
            simul.Tirate[t,i] = simul.Ti[t,i]/simul.Tk[t,i]
            simul.Tac[t,i] = simul.Tξ[t,i] + p.ϕ_k*abs(simul.Ti[t,i])^2/(2*simul.Tk[t,i]) + p.S*(simul.Tirate[t,i]<0)*simul.Ti[t,i]
            simul.Td[t,i] = PLFprofit[simul.Tz[t,i]](simul.Tk[t,i]) - simul.Ti[t,i] - simul.Tac[t,i]
            simul.Tvf[t,i] = VF[simul.Tz[t,i]](simul.Tk[t,i])
        end
    end

    ###simulate output [sales]
    for i in 1:1:N
        for t in 1:Tss+1
            simul.Tl[t,i] = PLFL[simul.Tz[t,i]](simul.Tk[t,i])
            simul.Ty[t,i] = ss.Pss[4]*simul.Tz_value[t,i] * simul.Tk[t,i]^p.α * simul.Tl[t,i]^p.ν
            #simul.Ty[t,i] = simul.Tz_value[t,i] * simul.Tk[t,i]^p.α * simul.Tl[t,i]^p.ν
        end
    end

    ####---------------------------------------------------------------------------------------------------------------
    #Part2: Transition Simulation [with potential shocks] [T: Tss+1--->Tss+Ttr]
    for ttr in Tss+2:1:Ttt

        t = ttr - Tss

        ###Uncertainty Process
        U = tr.TU[t-1]
        u = Uncertainty(p,U)

        ###Generate Interpolations
        PLFA = [interpolate(p.kgrid, tr.TkpolicyA[t][1][:,iz],FritschButlandMonotonicInterpolation()) for iz in 1:p.Nz]
        PLFNA = [interpolate(p.kgrid, tr.TkpolicyNA[t][1][:,iz],FritschButlandMonotonicInterpolation()) for iz in 1:p.Nz]
        PLFξ = [interpolate(p.kgrid, tr.Tξstar[t][1][:,iz],FritschButlandMonotonicInterpolation()) for iz in 1:p.Nz]
        PLFL = [interpolate(p.kgrid, tr.Tlpolicy[t][1][:,iz], FritschButlandMonotonicInterpolation()) for iz in 1:p.Nz]
        VF = [interpolate(p.kgrid, tr.Tvf[t][1][:,iz], FritschButlandMonotonicInterpolation()) for iz in 1:p.Nz]
        PLFprofit = [interpolate(p.kgrid, tr.Tprofit[t][1][:,iz], FritschButlandMonotonicInterpolation()) for iz in 1:p.Nz]

        ###simulate shocks
        for i in 1:1:N
            simul.Tz[ttr,i] = simulate(u.mc_z,2;init=simul.Tz[ttr-1,i])[2] #simulate just one period ahead for next period which is [2] from QuantEcon
            simul.Tz_value[ttr,i] = exp(p.zgrid[simul.Tz[ttr,i]]+u.Δz[simul.Tz[ttr,i]])
        end
        for i in 1:1:N
            simul.Tξ[ttr,i] = rand()*p.ξbar
        end

        ###simulate capital
        for i in 1:1:N
            if simul.Tξ[ttr-1,i] > PLFξ[simul.Tz[ttr-1,i]](simul.Tk[ttr-1,i])
                simul.Tk[ttr,i] =  max(p.kgrid[1],PLFNA[simul.Tz[ttr,i]](simul.Tk[ttr-1,i]))
            else
                simul.Tk[ttr,i] = PLFA[simul.Tz[ttr,i]](simul.Tk[ttr-1,i])
            end
            simul.Ti[ttr,i] = simul.Tk[ttr,i] - simul.Tk[ttr-1,i]
            simul.Tirate[ttr,i] = simul.Ti[ttr,i]/simul.Tk[ttr,i]
            simul.Tac[ttr,i] = simul.Tξ[ttr,i] + p.ϕ_k*abs(simul.Ti[ttr,i])^2/(2*simul.Tk[ttr,i]) + p.S*(simul.Tirate[ttr,i]<0)*simul.Ti[ttr,i]
            simul.Td[ttr,i] = PLFprofit[simul.Tz[ttr,i]](simul.Tk[ttr,i]) - simul.Ti[ttr,i] - simul.Tac[ttr,i]
            simul.Tvf[ttr,i] = VF[simul.Tz[ttr,i]](simul.Tk[ttr,i])
        end

        ###simulate output [sales]
        for i in 1:1:N
            simul.Tl[ttr,i] = PLFL[simul.Tz[ttr,i]](simul.Tk[ttr,i])
            simul.Ty[ttr,i] = tr.Tpw[t]*simul.Tz_value[ttr,i] * simul.Tk[ttr,i]^p.α * simul.Tl[ttr,i]^p.ν
            #simul.Ty[ttr,i] = simul.Tz_value[ttr,i] * simul.Tk[ttr,i]^p.α * simul.Tl[ttr,i]^p.ν
        end

    end

    #index public firms by cutoff: 10% (17%) largest firms account for 10/22=45% of the total output [2019 data]
    for t in 1:Ttt
        simul.Tipo[t] = quantile(simul.Ty[t,:], 0.90)
        for i in 1:1:N
            if simul.Ty[t,i] > simul.Tipo[t]
                simul.Tpub[t,i] = 1
            end
        end
    end

    #calculate the age of being a public firm
    for i in 1:N
        if simul.Tpub[1,i] == 1
            simul.Tpubage[1,i] = 1
        else
            simul.Tpubage[1,i] = 0
        end
        for t in 2:Ttt
            if simul.Tpub[t-1,i] == 1
                simul.Tpubage[t,i] = simul.Tpubage[t-1,i] + 1
            else
                simul.Tpubage[t,i] = 0
            end
        end
    end

    for i in 1:N
        for t in 1:Ttt
            if simul.Tpubage[t,i] >= 100
                simul.Tpub25plus[t,i] = 1
            end
        end
    end

    #Calculate sales growth
    #Calculate stock return
    for i in 1:N
        for t in 5:Ttt
            simul.Tsalesgrowth[t,i] = 2*(simul.Ty[t,i] - simul.Ty[t-4,i])/(simul.Ty[t,i] + simul.Ty[t-4,i])
            simul.Tstockreturn[t,i] = simul.Tvf[t,i]/(simul.Tvf[t-1,i]-simul.Td[t-1,i])
        end
    end

    for t in 1:1:Ttt
    #    simul.TIQRsg[t] = iqr(simul.Tsalesgrowth[t,:])
         #z only took 5 points, so sales growth are accumulated around -1.0, -0.5, 0, 0.5, 1.0
         #with volatility shock, the IQR measure jumped from 0.07 [-0.035->0.035] to [-0.5->0.5]
         #while SD only jump from 0.16 to 0.45, because the distribution is approximately normal
         #I use iqr=1.35*sd which is more reliable
         simul.TIQRsg[t] = 1.35*std(simul.Tsalesgrowth[t,:])
         simul.TIQRsr[t] = 1.35*std(simul.Tstockreturn[t,:])
    end

end
