#This File include several functions to solve firms' problem given aggregate
#prices and uncertainty state S, files including:
#   Firm_VFI
#   laborpolicy
#   expectation
#   capitalpolicy
#   howard
"""
Value Function Iteration which solves the firm problem
------------------------------------------------------
##### Arguments
- `p` : Parameters
- `pw`: whole sale price
- `w`: wage
- `Λ`: stochastic discount factor
- `S`: state of volatility
------------------------------------------------------
"""
function Firm_VFI(p::NamedTuple, U::Int, pw::TF, w::TF, Λ::TF, q::TF;
              A::TF=1.0,
              tol::TF = 1e-8,
              maxit::Int = 1500,
              it::Int = 0,
              dist::TF = 10.0,
              howard_on::Bool = false,
              V_upd::Array{TF,2} = similar(p.vf)
              ) where TF<:AbstractFloat
    #Tracking the progress
    prog = ProgressThresh(tol, "Solving individual problem by VFI: ")
    #Compute Labor Policy
    lpolicy = laborpolicy(p, U, pw, w; A)
    for ik in 1:p.Nk, iz in 1:p.Nz
        p.lpolicy[ik,iz]=lpolicy[ik,iz]
    end
    #Compute Profits
    profit = profit_func(p, U, pw, w; A)
    for ik in 1:p.Nk, iz in 1:p.Nz
        p.profit[ik,iz]=profit[ik,iz]
    end
    #Main Loop to converge value and policy function
    while dist > tol && it < maxit
        #Updating iteration and value function
        it += 1
        V_upd = deepcopy(p.vf)
        #Calculate Expectation
        u = Uncertainty(p,U)
        EV = expectation(p, p.vf, u.Πz)
        #Update Value Function p.vf & Policy Function p.kpolicy
        vf, kpolicyA, kpolicyNA, ξstar, vfA, vfNA = bellman(p,pw,w,Λ,q,EV,p.profit)
        for ik in 1:p.Nk, iz in 1:p.Nz
            p.vf[ik,iz]= vf[ik,iz]
            p.kpolicyA[ik,iz]= kpolicyA[ik,iz]
            p.kpolicyNA[ik,iz]= kpolicyNA[ik,iz]
            p.ξstar[ik,iz] = ξstar[ik,iz]
            p.vfA[ik,iz] = vfA[ik,iz]
            p.vfNA[ik,iz] = vfNA[ik,iz]
        end
        #Howard Improvement
        if howard_on
            howard(p,U,pw,w,Λ)
        end
        #Updating Distance and the Progress
        dist = maximum(abs, V_upd - p.vf)
        ProgressMeter.update!(prog, dist)
    end
end




"""
Compute the optimal labor policy given wage (Using FOC of firm's choice)
------------------------------------------------------------------------
##### Arguments
- `p` : Param (using grids of k, z)
- `pw`: wholesale price
- `w`: wage
##### Output
- `lpolicy`: labor policy on an [k,z] matrix
------------------------------------------------------------------------
"""
function laborpolicy(p::NamedTuple, U::Int, pw::TF, w::TF;
                     A::TF=1.0,
                     lpolicy::Array{TF,2}=zeros(p.Nk,p.Nz)) where TF<:AbstractFloat

    u = Uncertainty(p,U)

    for (iz, z) = enumerate(p.zgrid)

        prod = exp(z+u.Δz[iz])

        for (ik, k) = enumerate(p.kgrid)
            # Euler Equation to update l
            lpolicy[ik,iz] =  (p.ν * A * prod * k^p.α * pw / w)^(1/(1-p.ν))
        end
    end
    return lpolicy
end


"""
Compute the profit of each grid (Using FOC of firm's choice)
------------------------------------------------------------------------
##### Arguments
- `p` : Param (using grids of k, z)
------------------------------------------------------------------------
"""
function profit_func(p::NamedTuple, U::Int, pw::TF, w::TF;
                A::TF=1.0,
                profit::Array{TF,2}=zeros(p.Nk,p.Nz)) where TF<:AbstractFloat

    u = Uncertainty(p,U)

    for (iz, z) = enumerate(p.zgrid)

        prod = exp(z+u.Δz[iz])

        for (ik, k) = enumerate(p.kgrid)
            # Euler Equation to update l
            profit[ik,iz] = (pw*A*prod*k^p.α)^(1/(1-p.ν)) * w^(-p.ν/(1-p.ν)) * ( p.ν^(p.ν/(1-p.ν))-p.ν^(1/(1-p.ν)) )
        end

    end

    return profit

end


"""
Compute Expected Value Function for next period
-------------------------------------------------
##### Arguments
- `p` : Parameters (using grids of k, z, S)
- `S` : Uncertainty State [1(low) or 2(high)]
-------------------------------------------------
"""
function expectation(p::NamedTuple, vfprime::Array{TF,2}, Πz::Array{TF,2};
                      EV::Array{TF} = zeros(p.Nk,p.Nz)) where TF<:AbstractFloat
    for ikprime in 1:p.Nk
        for iz in 1:p.Nz
            for izprime in 1:p.Nz
                EV[ikprime,iz] += Πz[iz,izprime] * vfprime[ikprime,izprime]
            end
        end
    end
    return EV
end


"""
Update Value Function once using bellman equation
-------------------------------------------------
##### Arguments
- `p` : Param (using grids of k, z, S, lpolicy)
- `EV` : expected value function for next period
-------------------------------------------------
"""
function bellman(p::NamedTuple, pw::TF, w::TF, Λ::TF, q::TF, EV::Array{TF}, profit::Array{TF};
                 vf::Array{TF,2}=zeros(p.Nk,p.Nz),
                 vfNA::Array{TF,2}=zeros(p.Nk,p.Nz),
                 vfA::Array{TF,2}=zeros(p.Nk,p.Nz),
                 ξstar::Array{TF,2}=zeros(p.Nk,p.Nz),
                 kpolicyA::Array{TF,2}=zeros(p.Nk,p.Nz),
                 kpolicyNA::Array{TF,2}=zeros(p.Nk,p.Nz)) where TF<:AbstractFloat
    for (iz, z) = enumerate(p.zgrid)
        EV_func = interpolate(p.kgrid, EV[:,iz], FritschButlandMonotonicInterpolation())
        #EV_func = Spline1D(p.kgrid, EV[:,iz])
        for (ik, k) = enumerate(p.kgrid)
            πl = profit[ik,iz]
            vf[ik,iz], kpolicyA[ik,iz], kpolicyNA[ik,iz], ξstar[ik,iz], vfA[ik,iz], vfNA[ik,iz] = capitalpolicy(p,w,k,πl,Λ,q,EV_func)
        end
    end
    return vf, kpolicyA, kpolicyNA, ξstar, vfA, vfNA
end

"""
Find Optimal Capital Policy Given states and prices
-------------------------------------------------
##### Arguments
- `p` : Param (using grids of k, z, S, lpolicy)
-
-------------------------------------------------
"""
function capitalpolicy(p::NamedTuple, w::TF, k::TF, πl::TF, Λ::TF, q::TF, EV_func;
                       res::Array{TF}=zeros(2,2)) where TF<:AbstractFloat

    #Step 1: Calculate vfNA
    #Case1: i<0
    function case1NA(kprime)
        i = (kprime - (1-p.δ)*k) * q
        AC = abs(i) * p.S + p.ϕ_k*abs(i)^2/(2*k)
        d = πl - i - AC
        return -(d + Λ*EV_func(kprime))
    end
    sol1 = optimize(case1NA, max(p.kgrid[1],(1-p.δ-p.a)*k), max(p.kgrid[1],(1-p.δ)*k))
    res[1,1] = sol1.minimizer
    res[1,2] = -sol1.minimum
    #Case2: i>0
    function case2NA(kprime)
        i = (kprime - (1-p.δ)*k) * q
        AC = p.ϕ_k*abs(i)^2/(2*k)
        d = πl - i - AC
        return -(d + Λ*EV_func(kprime))
    end
    #make sure of the range
    xmin = max( p.kgrid[1], (1-p.δ)*k ) # make sure xmin is in range
    xmax = max( xmin, min( (1-p.δ+p.a)*k, min(k+πl,p.kgrid[p.Nk]) ) ) # make sure xmax is in range
    sol2 = optimize(case2NA, xmin, xmax)
    res[2,1] = sol2.minimizer
    res[2,2] = -sol2.minimum

    #Return the optimal kprime and value function
    vfNA = findmax(res[:,2])[1]
    kpolicyNA = res[findmax(res[:,2])[2],1]

    #Step 2: Calculate vfA
    #        find kprime range that dividend < 0, this only happens when kprime is too large
    #        therefore, the only capital adjustment cost is F (in other words, i>0)
    #Case1: i<0
    function case1(kprime)
        i = (kprime - (1-p.δ)*k) * q
        AC = abs(i) * p.S + p.ϕ_k*abs(i)^2/(2*k)
        d = πl - i - AC
        return -(d + Λ*EV_func(kprime))
    end
    sol1 = optimize(case1, p.kgrid[1], max(p.kgrid[1],(1-p.δ)*k))
    res[1,1] = sol1.minimizer
    res[1,2] = -sol1.minimum
    #Case2: i>0
    function case2(kprime)
        i = (kprime - (1-p.δ)*k) * q
        AC = p.ϕ_k*abs(i)^2/(2*k)
        d = πl - i - AC
        return -(d + Λ*EV_func(kprime))
    end
    #make sure of the range
    xmin = max( p.kgrid[1], (1-p.δ)*k ) # make sure xmin is in range
    xmax = max( xmin, min(k+πl,p.kgrid[p.Nk])) # make sure xmax is in range
    sol2 = optimize(case2, xmin, xmax)
    res[2,1] = sol2.minimizer
    res[2,2] = -sol2.minimum

    #Return the optimal kprime and value function
    vfA = findmax(res[:,2])[1]
    kpolicyA = res[findmax(res[:,2])[2],1]

    ### Step 3:
    # Calculate ξstar
    ξstar = min( p.ξupbar, max( p.ξlowbar, (vfA-vfNA)/w ) )
    #Calculate vf
    vf = ( (ξstar-p.ξlowbar)/(2*sqrt(3)*p.σξ) )*(vfA - w*(ξstar+p.ξlowbar)/2) + ( 1 - (ξstar-p.ξlowbar)/(2*sqrt(3)*p.σξ) )*vfNA

    return vf, kpolicyA, kpolicyNA, ξstar, vfA, vfNA

end
