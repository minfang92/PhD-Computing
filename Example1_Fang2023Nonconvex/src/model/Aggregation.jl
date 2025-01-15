#This File include several functions to solve the aggregate distribution

"""
Find Stationary Distribution given aggregate variables
------------------------------------------------------
##### Arguments
- `p` : Parameters
- `agg`: Aggregate Variables
------------------------------------------------------
"""
function find_stationary_distributions(p::NamedTuple,
                                       S::Int;
                                       tol::TF = 1e-10,  #limit is 1e-17
                                       dist::TF = 10.0,
                                       it::Int = 1,
                                       maxit::Int = 1000,
                                       Dist_old = similar(p.Dist_kz),
                                       Dist_upd = deepcopy(p.Dist_kz)
                                       ) where TF<:AbstractFloat
    #Tracking the progress
    prog = ProgressThresh(tol, "Find stationary distribution by Young(2010):")

    #Main loop
    while dist > tol && it < maxit
        #Updating iteration and Distribution
        it += 1
        Dist_old = deepcopy(Dist_upd)
        #Non-stochastic Simulation, and update the distribution
        Dist_upd = Young(p,S,p.kpolicyA,p.kpolicyNA,p.ξstar,Dist_old)
        #Updating Distance and the Progress
        dist = maximum(abs, Dist_old - Dist_upd)
        ProgressMeter.update!(prog, dist)
    end

    #final update
    for iz in 1:p.Nz
        for ik in 1:p.Nk_dense
        p.Dist_kz[ik,iz] = Dist_upd[ik,iz]
        end
    end

end

"""
Young's method of Non-Stochastic Simulation to update distribution
Given Initial distribution and policy rule, update!
------------------------------------------------------
##### Arguments
- `p` : Parameters
- `agg`: Aggregate Variables
------------------------------------------------------
"""
function Young(p::NamedTuple,
               S::Int,
               kpolicyA::Array{TF,2},
               kpolicyNA::Array{TF,2},
               ξstar::Array{TF,2},
               Dist_kz::Array{TF,2};
               Dist_upd = zeros(p.Nk_dense, p.Nz))  where TF<:AbstractFloat
    #essential to start with 0, don't use "similar(agg.Nk_dense)*0", it creates "NaN"
    #for each state (ik,iz), find policy function, and distribute the evolution to corresponding grids

    u = Uncertainty(p,S)

    for iz in 1:p.Nz
        #interpolate policy function
        kpolicyA_func = interpolate(p.kgrid, p.kpolicyA[:,iz], FritschButlandMonotonicInterpolation())
        kpolicyNA_func = interpolate(p.kgrid, p.kpolicyNA[:,iz], FritschButlandMonotonicInterpolation())
        #kpolicy_func = Spline1D(p.kgrid, p.kpolicy[:,iz])
        ξstar_func = interpolate(p.kgrid, p.ξstar[:,iz], FritschButlandMonotonicInterpolation())
        #ξstar_func = Spline1D(p.kgrid, p.ξstar[:,iz])

        for (ik, k) = enumerate(p.kgrid_dense)
            #active share
            shareA = (ξstar_func(k)-p.ξlowbar)/(p.ξupbar-p.ξlowbar)
            #active kpolicy
            kpA = min(max(p.kgrid_dense[1],kpolicyA_func(k)),p.kgrid_dense[p.Nk_dense])
            kpindex_plusA = searchsortedfirst(p.kgrid_dense, kpA)
            kpindex_minusA = searchsortedlast(p.kgrid_dense, kpA)
            if kpindex_minusA == 0 print("Index Wrong Case 1") end
            kpweightA = (p.kgrid_dense[kpindex_plusA]-kpA)/p.kgap_dense
            #non-active kpolicy
            kpNA = min(max(p.kgrid_dense[1],kpolicyNA_func(k)),p.kgrid_dense[p.Nk_dense])
            kpindex_plusNA = searchsortedfirst(p.kgrid_dense, kpNA)
            kpindex_minusNA = searchsortedlast(p.kgrid_dense, kpNA)
            if kpindex_minusNA == 0 print("Index Wrong Case 2") end
            kpweightNA = (p.kgrid_dense[kpindex_plusNA]-kpNA)/p.kgap_dense

            #REDISTRIBUTE! #if kp is exactly on grid, it's fine, kpweight would be 0 or 1
            for jz in 1:p.Nz
                Dist_upd[kpindex_minusA,jz]  += u.Πz[iz,jz]*Dist_kz[ik,iz] * shareA     * kpweightA
                Dist_upd[kpindex_plusA,jz]   += u.Πz[iz,jz]*Dist_kz[ik,iz] * shareA     * (1-kpweightA)
                Dist_upd[kpindex_minusNA,jz] += u.Πz[iz,jz]*Dist_kz[ik,iz] * (1-shareA) * kpweightNA
                Dist_upd[kpindex_plusNA,jz]  += u.Πz[iz,jz]*Dist_kz[ik,iz] * (1-shareA) * (1-kpweightNA)
            end

        end

    end

    return Dist_upd

end


"""
From Stationary Distribution to aggregate variables
------------------------------------------------------
##### Arguments
- `p` : Parameters
- `s`: uncertainty state of the economy
- `Dist_kz`: Distribution of firms over the (k,z) space
- `lpolicy`: laborpolicy over the (k,z) space
- `kpolicy`: capitalpolicy over the (k,z) space
##### Output
- `N,Y,I,K,AC_K` : aggregate quantities
------------------------------------------------------
"""
function calculate_aggregation(p::NamedTuple,
                               U::Int,
                               Dist_kz::Array{TF,2},
                               ξstar::Array{TF,2},
                               lpolicy::Array{TF,2},
                               kpolicyA::Array{TF,2},
                               kpolicyNA::Array{TF,2};
                               A::TF=1.0,
                               ydist::Array{TF,2} = zeros(p.Nk_dense, p.Nz),
                               ldist::Array{TF,2} = zeros(p.Nk_dense, p.Nz),
                               idist::Array{TF,2} = zeros(p.Nk_dense, p.Nz),
                               idistA::Array{TF,2} = zeros(p.Nk_dense, p.Nz),
                               idistNA::Array{TF,2} = zeros(p.Nk_dense, p.Nz),
                               ackdist::Array{TF,2} = zeros(p.Nk_dense, p.Nz)
                               )where TF<:AbstractFloat

   u = Uncertainty(p,U)

   #Calculate through the distribution
   for (iz, z) = enumerate(p.zgrid)

       prod = A*exp(z+u.Δz[iz])

       kpolicyA_func = interpolate(p.kgrid, kpolicyA[:,iz], FritschButlandMonotonicInterpolation())
       kpolicyNA_func = interpolate(p.kgrid, kpolicyNA[:,iz], FritschButlandMonotonicInterpolation())
       #kpolicy_func = Spline1D(p.kgrid, kpolicy[:,iz])
       lpolicy_func = interpolate(p.kgrid, lpolicy[:,iz], FritschButlandMonotonicInterpolation())
       #lpolicy_func = Spline1D(p.kgrid, lpolicy[:,iz])
       ξstar_func = interpolate(p.kgrid, ξstar[:,iz], FritschButlandMonotonicInterpolation())
       #ξstar_func = Spline1D(p.kgrid, ξstar[:,iz])

       for (ik, k) = enumerate(p.kgrid_dense)
           shareA = (ξstar_func(k)-p.ξlowbar)/(p.ξupbar-p.ξlowbar)
           ldist[ik,iz] =  lpolicy_func(k)
           ydist[ik,iz] =  prod * k^p.α * ldist[ik,iz]^p.ν
           idistA[ik,iz] =  ( kpolicyA_func(k) - (1-p.δ)*k ) * shareA
           idistNA[ik,iz] =  ( kpolicyNA_func(k) - (1-p.δ)*k ) * (1-shareA)
           idist[ik,iz] =  idistA[ik,iz] + idistNA[ik,iz]
           ackdist[ik,iz] = ( abs(idist[ik,iz])*p.S*(idist[ik,iz]<0) + 0.5*p.ϕ_k*(idist[ik,iz]^2/k) )
       end
   end

   #Aggregation
   N = sum(ldist .* Dist_kz)
   Y = sum(ydist .* Dist_kz)
   IA = sum(idistA .* Dist_kz)
   INA = sum(idistNA .* Dist_kz)
   I = sum(idist .* Dist_kz)
   K = sum(p.kgrid_dense .* sum(Dist_kz,dims=2))
  AC_K = sum(ackdist .* Dist_kz)

  return N, Y, I, IA, INA, K, AC_K

end
