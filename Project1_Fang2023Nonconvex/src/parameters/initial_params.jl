"""
This function holds all parameters and grids for the model
"""
function Param( ;
    #### Parameters
    #Preferences and Technology
    β = 0.99,               #Annual discount factor of 95%
    η = 1.0,                #Unit elasticity of intertemporal substitution (Khan and Thomas (2008))
    θ = 2.0,                #Leisure preference, households spend 1/3 of time working
    χ = 1.0,                #Infinite Frisch elasticity of labor supply (Khan and Thomas (2008))
    α = 0.25,               #DRS production, imply a total returns to scale of 85%
    ν = 0.60,               #DRS labor share of 2/3, capital share of 1/3 Karabarbounis and Neiman (2013)
    #Adjustment Costs
    δ = 0.025,
    μξ = 0.2,
    σξ = 0.1,
    S = 0.0,                #Resale Costs
    ϕ_k = 2.4,              #Quad Adjustment Costs
    a = 0.001,               #free region
    #Uncertainty Process
    Nz = 5,                 #State of Idiosyncratic prod. process
    ρz = 0.95,              #Quarterly persistence of idiosyncratic productivity (Khan and Thomas (2008))
    σz = 0.05,              #Bloom et al.(2018)
    #New Keynesian Block
    γ = 10,                   #Demand elasticity
    ϕ_π = 1.5,               #Taylor rule coefficient
    ϕ_y = 0.1,
    ψ = 90,                   #Price Adjustment Costs
    γ_r = 0.76 )              #inertial taylor rule

    #param
    ξupbar = μξ + sqrt(3)*σξ
    ξlowbar = μξ - sqrt(3)*σξ

    ### Grids
    #Transiction Matrix and z grids [Tauchen]
    Πz = tauchen(Nz, ρz, σz).p
    zgrid = collect(tauchen(Nz, ρz, σz).state_values)
    mc_z = MarkovChain(Πz)

    #Capital Space
    Nk    = 51
    kmin  = 0.21
    # kmax  = kmin + 0.05*Nk
    # kmax  = exp(log(kmin) - (Nk-1)*log(1.0-δ)) #4.2 for Nk = 50
    kmax = 15.1
    logkmin = log(kmin)
    logkmax = log(kmax)
    logkgrid = collect(range(logkmin, stop=logkmax, length=Nk))
    kgrid = collect(exp.(logkgrid))

    #Individual Solution
    vf = zeros(Nk, Nz)
    vfA = zeros(Nk,Nz)
    vfNA = zeros(Nk,Nz)
    ξstar = zeros(Nk,Nz)
    profit = zeros(Nk, Nz)
    kpolicyA = zeros(Nk, Nz)
    kpolicyNA = zeros(Nk, Nz)
    lpolicy = zeros(Nk, Nz)
    #dividend = zeros(Nk,Nz)

    #Aggregate Distribution[dense grids needed for accurency]
    Nk_dense = Nk * 20
    kgrid_dense = collect(range(kgrid[1], stop=kgrid[Nk], length=Nk_dense))
    kgap_dense = kgrid_dense[2]-kgrid_dense[1]
    Dist_kzA = 0.5*ones(Nk_dense, Nz)/(Nk_dense*Nz)
    Dist_kzNA = 0.5*ones(Nk_dense, Nz)/(Nk_dense*Nz)
    Dist_kz = ones(Nk_dense, Nz)/(Nk_dense*Nz)

    ###Output
    p =  (β=β,η=η,θ=θ,χ=χ,α=α,ν=ν,
          δ=δ,ξupbar=ξupbar,ξlowbar=ξlowbar,S=S,ϕ_k=ϕ_k,a=a,μξ=μξ,σξ=σξ,
          Nz=Nz,ρz=ρz,σz=σz,
          γ=γ,ϕ_π=ϕ_π,ϕ_y=ϕ_y,ψ=ψ,γ_r=γ_r,
          Πz=Πz,zgrid=zgrid,mc_z=mc_z,
          Nk=Nk,kmin=kmin,kmax=kmax,kgrid=kgrid,
          vf=vf,vfA=vfA,vfNA=vfNA,ξstar=ξstar,profit=profit,kpolicyA=kpolicyA,kpolicyNA=kpolicyNA,lpolicy=lpolicy,
          Nk_dense=Nk_dense,kgrid_dense=kgrid_dense,kgap_dense=kgap_dense,
          Dist_kzA=Dist_kzA,Dist_kzNA=Dist_kzNA,Dist_kz=Dist_kz)

end
