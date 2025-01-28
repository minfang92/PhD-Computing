#This File include several functions to solve dynamics given aggregate
#prices and uncertainty state S, files including:
#  Transition
#  Backward
#  Forward
#  Newton Method
#  find_transition_path
"""
This function holds all aggregate variables and distribution at Transition
"""
function Transition(p::NamedTuple,ss::NamedTuple,T::Int)

    #Aggregate Prices
    Tπ = ones(T)*0.0
    Tr = ones(T)*(1.0/p.β - 1.0)
    TΛ = ones(T)*p.β
    Tpw = ones(T)*ss.Pss[4]
    Tw = ones(T)*ss.Pss[5]

    #Aggregate Quantities
    TC = ones(T)*ss.Qss[1]
    TN = ones(T)*ss.Qss[2]
    TY = ones(T)*ss.Qss[3]
    TK = ones(T)*ss.Qss[4]
    TI = ones(T)*p.δ*ss.Qss[4]
    TAC_K = ones(T)*ss.Qss[8]
    TAC_P = zeros(T)*ss.Qss[9]

    #Aggregate State
    TU = ones(Int,T) * ss.U
    Tϵm = zeros(T)
    Tq = ones(T)
    TA = ones(T)

    #Choices
    Tvf = [[ss.vfss] for j = 1:T]      #!!![terminal value should be steadystate]
    TvfA = [[ss.vfAss] for j = 1:T]      #!!![terminal value should be steadystate]
    TvfNA = [[ss.vfNAss] for j = 1:T]      #!!![terminal value should be steadystate]
    TEV = [[ss.EVss] for j = 1:T]      #initial and terminal value doesnt matter

    Tξstar = [[ss.ξstarss] for j = 1:T]
    TkpolicyA = [[ss.kpolicyAss] for j = 1:T]
    TkpolicyNA = [[ss.kpolicyNAss] for j = 1:T]
    Tlpolicy = [[ss.lpolicyss] for j = 1:T]
    Tprofit = [[ss.profitss] for j = 1:T]  #initial and terminal value doesnt matter

    #Distribution !!![initial value should be steadystate]
    TDist_kz = [[ss.Dist_kz] for j = 1:T]

    ###Output
    tr = (T = T,
          Tπ = Tπ,  Tr = Tr,  TΛ = TΛ,  Tpw = Tpw,  Tw = Tw,
          TC = TC,  TN = TN,  TY = TY,  TI = TI, TK = TK,  TAC_K = TAC_K, TAC_P = TAC_P,
          TU = TU,  Tϵm = Tϵm, Tq = Tq, TA = TA,
          Tvf = Tvf,  TvfA = TvfA, TvfNA = TvfNA, TEV = TEV,
          Tξstar = Tξstar, TkpolicyA = TkpolicyA, TkpolicyNA = TkpolicyNA, Tlpolicy = Tlpolicy,  Tprofit = Tprofit,
          TDist_kz = TDist_kz
          )

end


"""
Update Policy Function BACKWARD using bellman equation
-------------------------------------------------
##### Arguments
- `p` : Param (using grids of k, z, S, lpolicy)
- `tr` : Transition
-------------------------------------------------
"""
function backward(p::NamedTuple,tr::NamedTuple,TΛ::Array{TF,1},Tpw::Array{TF,1},Tw::Array{TF,1},TU::Array{Int,1};
                  TA::Array{TF,1}=tr.TA, Tq::Array{TF,1}=tr.Tq) where TF<:AbstractFloat
    # Backward calculation from T-1 (T is known as SteadyState)
    for t in (tr.T-1):-1:2
        # Calculate Labor Policy, Profit, and Expectation
        tr.Tlpolicy[t][1] = laborpolicy(p,TU[t-1],Tpw[t],Tw[t];A=TA[t])
        tr.Tprofit[t][1] = profit_func(p,TU[t-1],Tpw[t],Tw[t];A=tr.TA[t])
        u = Uncertainty(p,TU[t-1])
        tr.TEV[t][1] = expectation(p,tr.Tvf[t+1][1],u.Πz)

        # Calculate Value Function and Capital Policy
        tr.Tvf[t][1], tr.TkpolicyA[t][1], tr.TkpolicyNA[t][1], tr.Tξstar[t][1], tr.TvfA[t][1], tr.TvfNA[t][1] = bellman(p,Tpw[t],Tw[t],TΛ[t],tr.Tq[t],tr.TEV[t][1],tr.Tprofit[t][1])

    end
end


"""
Update Distribution FORWARD using policy function
-------------------------------------------------
##### Arguments
- `p` : Param (using grids of k, z, S, lpolicy)
- `tr` : Transition
-------------------------------------------------
"""
function forward(p::NamedTuple,tr::NamedTuple,Tpw::Array{TF,1},Tw::Array{TF,1},TU::Array{Int,1};
                 TA::Array{TF,1}=tr.TA ) where TF<:AbstractFloat

    # Forward simulation from t=2 (t=1 is known as SteadyState)
    for t in 2:1:(tr.T-1)
        tr.TDist_kz[t][1] = Young(p,tr.TU[t-1],tr.TkpolicyA[t-1][1],tr.TkpolicyNA[t-1][1],tr.Tξstar[t-1][1],tr.TDist_kz[t-1][1])
    end

    # Calculate Aggregate Quantities [including initial and terminal]
    for t in 2:1:(tr.T-1)
        # From aggregation of distribution
        N, Y, I, IA, INA, K, AC_K = calculate_aggregation(p,TU[t],tr.TDist_kz[t][1],tr.Tξstar[t][1],
                                                          tr.Tlpolicy[t][1],tr.TkpolicyA[t][1],tr.TkpolicyNA[t][1];
                                                          A=TA[t])
        tr.TN[t] = N
        tr.TY[t] = Y
        tr.TK[t] = K
        tr.TI[t] = I
        tr.TAC_K[t] = AC_K
        tr.TAC_P[t] = p.ψ/2*tr.Tπ[t]^2*tr.TY[t]
        tr.TC[t] = tr.TY[t] - tr.TI[t] - tr.TAC_K[t] - tr.TAC_P[t]
    end

end


"""
Find Transition Path in Real Business Cycle (non-NK)
%  This program
-------------------------------------------------
##### Arguments
- `p` : Param (using grids of k, z, S, lpolicy)
- `tr` : Transition
-------------------------------------------------
"""
function find_transition_RBC(p::NamedTuple, tr::NamedTuple;
                              T::Int = tr.T,
                              tol::TF = 1e-6,
                              update::TF = 0.05,
                              maxit::Int = 1000,
                              it_out::Int = 0,
                              dist_out::TF = 10.0,
                              ExDemand::Array{TF,1} = tr.TY) where TF<:AbstractFloat

    ### Target of Updating
    TC = deepcopy(tr.TC)

    #Tracking the progress
    prog = ProgressThresh(tol, "Solving RBC GE Transition: ")

    ### Loop to Achieve Target
    while it_out < maxit && dist_out > tol
        #Step0: update iteration
        it_out += 1

        #Step1: Backward solving decision rules given aggregate prices
        backward(p, tr, tr.TΛ, tr.Tpw, tr.Tw, tr.TU; TA=tr.TA, Tq = tr.Tq)

        #Step2: Forward solving distribution
        forward(p,tr,tr.Tpw,tr.Tw,tr.TU; TA=tr.TA)

        #Step3: Update TC &
        TCp = tr.TY - tr.TI - tr.TAC_K
        ExDemand = (TCp .- TC)./TC
        dist_out = maximum(abs,ExDemand)
        ProgressMeter.update!(prog, dist_out)

        TC = TC .* (1 .+ update*ExDemand)
        for t in tr.T-1:-1:2
            tr.TΛ[t] = p.β * (TC[t]/TC[t+1])^p.η
            tr.Tw[t] = p.θ * TC[t]^p.η
        end

    end

    ### Update Target
    for t in 2:tr.T-1
        tr.TC[t] = TC[t]
    end

end
