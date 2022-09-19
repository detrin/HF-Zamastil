
import CGcoefficient
using ClassicalOrthogonalPolynomials
using QuadGK
import Infinity
using Polynomials
include("./fock_elements.jl")
include("./state_overlaps.jl")
include("./spherical_integrals.jl")
include("./states.jl")


l, m, θ, ϕ = 0, 0, 0, 0
println(Y(l, m, θ, ϕ))

n1, l1, m1 = 1, 1, 0
n2, l2, m2 = 1, 1, 0
integral, err = quadgk(
    r -> begin
    integral, err = quadgk(
        θ -> begin
        integral, err = quadgk(
            ϕ -> state(n1, l1, m1, r, θ, ϕ) * adjoint(state(n2, l2, m2, r, θ, ϕ)), 
            0, 2*π, rtol=1e-10
        );
        integral
        end,
        0, π, rtol=1e-10
    );
    integral
    end,
    0, Infinity.∞, rtol=1e-10
)
println(integral)
