#This File include several functions to generate transition matrix of indiosyncratic shocks
#Transition Matrix normal ΠzL, high volatility ΠzH, tail risk ΠzT
#files including:
#   Firm_VFI
#   laborpolicy
#   expectation
#   capitalpolicy
#   howard


"""
This function defines Volatility Process using Tauchen (1986) adapted to stochastic volatility;
---External Functions Needed
* erfc() from SpecialFunctions
* MarkovChain(Π[state], y) from QuantEcon
---Input
* N: number of discretized space grids
* ρ: persistence of AR(1) shocks
* σ: A vector of volatility
* μ: mean of AR(1) process
---Output
* An array of MarkovChain
"""
function stochastic_volatility_tauchen(N::Int, ρ::Real, σ::Vector, μ::Real=0.0, n_std::Int=3) #n_std=3 is the same as quantecon

   # Define necessary functions
    std_norm_cdf(x::T) where {T <: Real} = 0.5 * erfc(-x/sqrt(2))
    std_norm_cdf(x::Array{T}) where {T <: Real} = 0.5 .* erfc(-x./sqrt(2))

    # Get discretized space grids
    a_bar = n_std * sqrt(σ[1]^2 / (1 - ρ^2))
    y = range(-a_bar, stop=a_bar, length=N)
    d = y[2] - y[1]

    # Get transition probabilitie
    Π = [ones(N, N).*1/N for state = 1:length(σ)]
    mc = [MarkovChain(Π[state], y) for state = 1:length(σ)]
    for state = 1:length(σ)
       for row = 1:N
          # Do end points first
          Π[state][row, 1] = std_norm_cdf((y[1] - ρ*y[row] + d/2) / σ[state])
          Π[state][row, N] = 1 - std_norm_cdf((y[N] - ρ*y[row] - d/2) / σ[state])
          # fill in the middle columns
          for col = 2:N-1
              Π[state][row, col] = (std_norm_cdf((y[col] - ρ*y[row] + d/2) / σ[state]) -
                           std_norm_cdf((y[col] - ρ*y[row] - d/2) / σ[state]))
          end
       end
       # Renormalize. In some test cases the rows sum to something that is 2e-15
       # away from 1.0, which caused problems in the MarkovChain constructor
       yy = y .+ μ / (1 - ρ) # center process around its mean (wbar / (1 - rho)) in new variable
       Π[state] = Π[state]./sum(Π[state], dims = 2)
       mc[state] = MarkovChain(Π[state], yy)
    end

    return mc

end

"""
This function defines Tail Risk using
"""
function tail_risk(N::Integer, T::Real, λ_T::Real, zgrid)

    TR = zeros(N,N)
    for i = 1:N
      TR[i,1] =  T*((zgrid[N]-zgrid[i])/(zgrid[N]-zgrid[1]))^λ_T
      TR[i,i] = -T*((zgrid[N]-zgrid[i])/(zgrid[N]-zgrid[1]))^λ_T
   end
   TR[1,1] = 0.0
   return TR
end
