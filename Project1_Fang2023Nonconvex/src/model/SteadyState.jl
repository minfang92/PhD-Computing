#This File include several functions to solve the stedaystate
"""
This function holds all aggregate variables and distribution at SteadyState
"""
function SteadyState(p::NamedTuple;

    U = 0,

    #Aggregate Prices[Analytically][initial value equals steadystate]
    w = 1.0,
    π = 0.0,            #inflation [logΠ]
    r = 1.0/p.β - 1.0,  #interest [logR]
    Λ = p.β,            #discount factor
    pw = 1.0,           #wholesale price
    q = 1.0 )           #investment price

    #Aggregate Quantities [Too be pin down]
    C = 1.0
    N = 1.0
    Y = 1.0
    K = 1.0
    I = 0.1
    IA = 0.05
    INA = 0.05
    AC_K = 0.01
    AC_P = 0.00

    #Array Form could be updated by function
    Pss = [π r Λ pw w]
    Qss = [C N Y K I IA INA AC_K AC_P]

    #value Function and Policy Function
    vfss = zeros(p.Nk,p.Nz)
    vfAss = zeros(p.Nk,p.Nz)
    vfNAss = zeros(p.Nk,p.Nz)
    EVss = zeros(p.Nk,p.Nz)
    ξstarss = zeros(p.Nk,p.Nz)
    kpolicyAss = zeros(p.Nk,p.Nz)
    kpolicyNAss = zeros(p.Nk,p.Nz)
    lpolicyss = zeros(p.Nk,p.Nz)
    profitss = zeros(p.Nk,p.Nz)

    #Aggregate Distribution[dense grids needed for accurency]
    Dist_kz = ones(p.Nk_dense, p.Nz)/(p.Nk_dense*p.Nz)

    ###Output
    ss = (U=U,q=q,
          Pss=Pss,Qss=Qss,
          vfss=vfss,vfAss=vfAss,vfNAss=vfNAss,EVss=EVss,
          ξstarss=ξstarss,kpolicyAss=kpolicyAss,kpolicyNAss=kpolicyNAss,
          lpolicyss=lpolicyss,profitss=profitss,
          Dist_kz=Dist_kz)

end




"""
This function finds the steadystate
"""
function find_SteadyState(p::NamedTuple,ss::NamedTuple;
                          updating::TF = 0.1,
                          howard_on::Bool = false,
                          tol::TF = 1e-8,
                          maxit::Int = 100,
                          it::Int = 0,
                          dist::TF = 10.0,
                          w_new::TF = ss.Pss[5],
                          w_upd::TF = ss.Pss[5]
                          ) where TF<:AbstractFloat
    ### Unpack the prices
    pw = ss.Pss[4]
    Λ = ss.Pss[3]
    U = ss.U
    q = ss.q

    ### Main Loop to converge value and policy function
    while dist > tol && it < maxit
        #Updating iteration and value function
        it += 1
        w_upd = (1-updating)*w_upd + updating*w_new
        #VFI given prices
        howard = howard_on
        Firm_VFI(p, U, pw, w_upd, Λ, q; tol=tol, howard_on=howard_on)
        #Find Stationary Equilibrium
        find_stationary_distributions(p,U)
        #Calculate Aggregate Variables
        N, Y, I, IA, INA, K, AC_k = calculate_aggregation(p,U,p.Dist_kz,p.ξstar,p.lpolicy,p.kpolicyA,p.kpolicyNA)
        #Update new wage
        C = Y-I-AC_k
        w_new = p.θ * C^p.η * N^(p.χ-1)
        #Updating Distance and the Progress
        @show dist = maximum(abs, w_upd - w_new)
        #Save Results to SteadyState
        #Pss = [π R Λ pw w]
        ss.Pss[5] = w_new
    end

    ### Update steadystate obejects
    update_steadystate(p,ss)

    ### Return
    return ss.Pss[5]



end


"""
This function updates the steadystate
"""
function update_steadystate(p::NamedTuple,ss::NamedTuple)

    U = ss.U
    u = Uncertainty(p,U)

    #1.policy + value
    for ik in 1:p.Nk
        for iz in 1:p.Nz
            ss.lpolicyss[ik,iz] = p.lpolicy[ik,iz]
            ss.kpolicyAss[ik,iz] = p.kpolicyA[ik,iz]
            ss.kpolicyNAss[ik,iz] = p.kpolicyNA[ik,iz]
            ss.profitss[ik,iz] = p.profit[ik,iz]
            ss.ξstarss[ik,iz] = p.ξstar[ik,iz]
            ss.vfss[ik,iz] = p.vf[ik,iz]
            ss.vfAss[ik,iz] = p.vfA[ik,iz]
            ss.vfNAss[ik,iz] = p.vfNA[ik,iz]
        end
    end

    #2.Expectations
    EV_func = expectation(p,ss.vfss,u.Πz)
    for ik in 1:p.Nk
        for iz in 1:p.Nz
            ss.EVss[ik,iz] = EV_func[ik,iz]
        end
    end

    #3.Distribution
    for ik in 1:p.Nk_dense
        for iz in 1:p.Nz
            ss.Dist_kz[ik,iz] = p.Dist_kz[ik,iz]
         end
    end

    #4.Aggregate Variables
    N, Y, I, IA, INA, K, AC_k = calculate_aggregation(p,U,ss.Dist_kz,ss.ξstarss,ss.lpolicyss,ss.kpolicyAss,ss.kpolicyNAss)
    #Calculate C #Qss = [C N Y K AC_K AC_P]
    C = Y-I-AC_k
    ss.Qss[1] = C
    ss.Qss[2] = N
    ss.Qss[3] = Y
    ss.Qss[4] = K
    ss.Qss[5] = I
    ss.Qss[6] = IA
    ss.Qss[7] = INA
    ss.Qss[8] = AC_k


end
