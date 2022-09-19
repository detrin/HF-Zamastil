
import CGcoefficient
using ClassicalOrthogonalPolynomials
using QuadGK
# using Infinity
using Polynomials

include("./spherical_integrals.jl")
include("./radial_integrals.jl")

function state(n, l, m, r, θ, ϕ; a0=1., η=1.)
    return Y(l, m, θ, ϕ) * R(n, l, η, r; a0=a0)
end