
import CGcoefficient
using ClassicalOrthogonalPolynomials
using QuadGK
import Infinity
using Polynomials
include("./fock_elements.jl")
include("./state_overlaps.jl")
include("./spherical_integrals.jl")
include("./states.jl")

n1j, lj, mj, n1k, lk, mk = 2, 1, -1, 2, 1, -1

println(two_states_num_integral(n1j, lj, mj, n1k, lk, mk; η=1.) )