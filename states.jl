
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

function get_states_discrete(n_sum)
    states = []
    for n2 in 1:(n_sum-1)
        for n1 in n2:(n_sum-n2)
            for l in 0:(min(n1, n2)-1)
                push!(states, [n2, n1, l, 1])
            end
        end
    end
    return states
end