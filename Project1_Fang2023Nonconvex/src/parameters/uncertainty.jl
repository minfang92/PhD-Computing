"""
This function holds all parameters and grids for uncertainty process
    ##### Types of Micro Uncertainty Shocks
    - micro volatility
    - micro tail risk
    - micro knight uncertainty (ambiguity aversion)
    [xxx forget this first] - micro noise
"""
function Uncertainty(p::NamedTuple, U::Int;
    #Uncertainty Shocks
    Δσ =  2.6,              #Micro Volatility increase in high uncertainty state
    Δσ4 = 1.7,
    T = 0.20,               #Micro Tail Risk increase in high uncertainty state
    λ_T = 0.0,              #Micro Tail Risk increase in high uncertainty state (distribution)
    a = 0.70                #Micro knight uncertainty increase in high uncertainty state: expected mean prod: 1 ---> 0.95
    )

    #Micro Volatility Shock
    Πz1 = stochastic_volatility_tauchen(p.Nz, p.ρz, [p.σz, p.σz*Δσ])[2].p
    #Δz1 = [ - (p.σz*Δσ)^2 / (2*(1-p.ρz)) for i in 1:p.Nz ]
    Δz1 = [ - (p.σz*Δσ)^2 / (2*(1-p.ρz^2)) for i in 1:p.Nz ]

    Πz4 = stochastic_volatility_tauchen(p.Nz, p.ρz, [p.σz, p.σz*Δσ4])[2].p
    Δz4 = [ - (p.σz*Δσ4)^2 / (2*(1-p.ρz^2)) for i in 1:p.Nz ]

    #Micro Tail Risk Shock
    TR = tail_risk(p.Nz, T, λ_T, p.zgrid)
    Πz2 = p.Πz + TR
    #Δz2 = - TR * p.zgrid
    #Δz2 = [ - (p.σz)^2 / (2*(1-p.ρz)) for i in 1:p.Nz ]
    Δz2 = [ - (p.σz)^2 / (2*(1-p.ρz^2)) for i in 1:p.Nz ]

    #Micro Ambiguity Shock
    Πz3 = p.Πz #
    #Δz3 = [ - (p.σz)^2 / (2*(1-p.ρz)) - a for i in 1:p.Nz ]
    Δz3 = [ - (p.σz)^2 / (2*(1-p.ρz^2)) - a for i in 1:p.Nz ]

    #No Uncertainty Shock
    if U == 0
        Πz = p.Πz
        #Δz = [ - (p.σz)^2 / (2*(1-p.ρz)) for i in 1:p.Nz ]
        Δz = [ - (p.σz)^2 / (2*(1-p.ρz^2)) for i in 1:p.Nz ]
        mc_z = MarkovChain(Πz)
    end

    #Micro Volatility Shock
    if U == 1
        Πz = Πz1
        Δz = Δz1
        mc_z = MarkovChain(Πz)
    end

    #Micro Tail Risk Shock
    if U == 2
        Πz = Πz2
        Δz = Δz2
        mc_z = MarkovChain(Πz)
    end

    #Micro Ambiguity Shock
    if U == 3
        Πz = Πz3
        Δz = Δz3
        mc_z = MarkovChain(Πz)
    end

    #Micro Volatility Shock
    if U == 4
        Πz = Πz4
        Δz = Δz4
        mc_z = MarkovChain(Πz)
    end

    ###Output
    u =  (Πz=Πz,Δz=Δz,mc_z=mc_z)

end
