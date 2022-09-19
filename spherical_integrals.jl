
import CGcoefficient

using Memoize
using QuadGK
using Infinity
using AssociatedLegendrePolynomials

function Y(l, m, θ, ϕ)
    # https://demonstrations.wolfram.com/HydrogenAtomRadialFunctions/
    # println(n, " ", l, " ", η, " ", r)
    #=
    if true
        N, err = quadgk(
            r -> r * (sqrt(factorial(n-1-l) / (2*n*factorial(n+l)^3)) *
            (2 / (n * a0))^(l + 3/2.) *
            η * (η*r)^l * exp(-(η*r)/(n*a0)) *
            laguerrel(n+l, 2*l+1, 2*(η*r)/(n*a0)) )^2,
            0, Infinity.∞, rtol=1e-8
        )
    end
    println(N)
    =#
    return λlm(l, m, cos(θ)) * exp(1im*m*ϕ)
end