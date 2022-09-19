
import CGcoefficient

using Memoize
using QuadGK
using Infinity
using AssociatedLegendrePolynomials

function Y(l, m, θ, ϕ)
    if m >= 0
        return λlm(l, m, cos(θ)) * exp(1im*m*ϕ)
    else
        m *= -1
        return λlm(l, m, cos(θ)) * exp(1im*m*ϕ) * (-1)^m * factorial(l-m) / factorial(l+m)
    end
end