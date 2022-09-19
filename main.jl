
import CGcoefficient
using ClassicalOrthogonalPolynomials
using QuadGK
# using Infinity
using Polynomials
include("./fock_elements.jl")
include("./state_overlaps.jl")
include("./spherical_integrals.jl")

n, l, η, r = 1, 0, 3, 0
println(R(n, l, η, r; a0=1.0))
println()

n1j, n2j, n1k, n2k = 1, 1, 1, 1
lj, lk = 0, 0
xi1j, xi2j, xi1k, xi2k = 1, 1, 1, 1

println(fock_r12(n1j, n2j, n1k, n2k, lj, lk, xi1j, xi2j, xi1k, xi2k; sign_j=+1, sign_k=+1))
println(fock_ham_one(n1j, n2j, n1k, n2k, lj, lk; η=0.9))
println()

coeff = 3.7714 
println("state_overlap")
states = [
    [1, 1, 0, 1],
    [1, 2, 0, 1],
    [1, 3, 0, 1],
    [2, 2, 0, 1],
    [2, 2, 1, -1]
]
S1 = zeros(Float64, length(states), length(states))
S2 = zeros(Float64, length(states), length(states))
for st1_i in 1:length(states)
    st1 = states[st1_i]
    n1j, n2j, lj, sign_j = st1
    for st2_i in 1:length(states)
        st2 = states[st2_i]
        n1k, n2k, lk, sign_k = st2
        # println(st1, " ", st2)
        val = state_overlap_old(n1j, n2j, n1k, n2k, lj, lk, xi1j, xi2j, xi1k, xi2k; sign_j=sign_j, sign_k=sign_k )
        # println(val)
        # println("---")
        # println(val * coeff)
        S1[st1_i, st2_i] = val 
        S2[st1_i, st2_i] = val * coeff
    end
end
println()
println(display(S1))
println(display(S2))
println()

coeff = 4.8550 / 0.625
println("fock_r12")
S1 = zeros(Float64, length(states), length(states))
S2 = zeros(Float64, length(states), length(states))
for st1_i in 1:length(states)
    st1 = states[st1_i]
    n1j, n2j, lj, sign_j = st1
    for st2_i in 1:length(states)
        st2 = states[st2_i]
        n1k, n2k, lk, sign_k = st2
        # println(st1, " ", st2)
        val = fock_r12(n1j, n2j, n1k, n2k, lj, lk, xi1j, xi2j, xi1k, xi2k; sign_j=sign_j, sign_k=sign_k )
        # println(val)
        # println("---")
        # println(val * coeff)
        S1[st1_i, st2_i] = val 
        S2[st1_i, st2_i] = val * coeff
    end
end
println()
println(display(S1))
println(display(S2))
println()

println("check")
println(fock_ham_one(n1j, n2j, n1k, n2k, lj, lk; η=0.9))
println(fock_ham_one(n1j, n2j, n1k, n2k, lj, lk; η=1))
println(fock_ham_one(n1j, n2j, n1k, n2k, lj, lk; η=2))

println("square_brackets")
n1, n2, n3, n4 = 1, 1, 1, 1
l1_, l2_, l3_, l4_ = 0, 0, 0, 0
xi1, xi2, xi3, xi4 = 1, 1, 1, 1
println(square_brackets(n1, n2, n3, n4, l1_, l2_, l3_, l4_, xi1, xi2, xi3, xi4))

n1, n2, n3, n4 = 2, 2, 2, 2
l1_, l2_, l3_, l4_ = 1, 1, 1, 1
xi1, xi2, xi3, xi4 = 1, 1, 1, 1
println(square_brackets(n1, n2, n3, n4, l1_, l2_, l3_, l4_, xi1, xi2, xi3, xi4))

n1, n2, n3, n4 = 2, 2, 2, 2
l1_, l2_, l3_, l4_ = 2, 2, 2, 2
xi1, xi2, xi3, xi4 = 1, 1, 1, 1
println(square_brackets(n1, n2, n3, n4, l1_, l2_, l3_, l4_, xi1, xi2, xi3, xi4))

println("round_brackets")
n1, n2, n3, n4 = 1, 1, 1, 1
l1_, l2_, l3_, l4_ = 0, 0, 0, 0
xi1, xi2, xi3, xi4 = 1, 1, 1, 1
println(round_brackets(n1, n2, n3, n4, l1_, l2_, l3_, l4_, xi1, xi2, xi3, xi4))

n1, n2, n3, n4 = 1, 1, 1, 2
l1_, l2_, l3_, l4_ = 0, 0, 0, 1
xi1, xi2, xi3, xi4 = 1, 1, 1, 1
println("(0, 0, 0, 1), ", round_brackets(n1, n2, n3, n4, l1_, l2_, l3_, l4_, xi1, xi2, xi3, xi4))

n1, n2, n3, n4 = 1, 1, 2, 1
l1_, l2_, l3_, l4_ = 0, 0, 1, 0
xi1, xi2, xi3, xi4 = 1, 1, 1, 1
println("(0, 0, 1, 0), ", round_brackets(n1, n2, n3, n4, l1_, l2_, l3_, l4_, xi1, xi2, xi3, xi4))

n1, n2, n3, n4 = 1, 2, 1, 1
l1_, l2_, l3_, l4_ = 0, 1, 0, 0
xi1, xi2, xi3, xi4 = 1, 1, 1, 1
println("(0, 1, 0, 0), ", round_brackets(n1, n2, n3, n4, l1_, l2_, l3_, l4_, xi1, xi2, xi3, xi4))

n1, n2, n3, n4 = 2, 1, 1, 1
l1_, l2_, l3_, l4_ = 1, 0, 0, 0
xi1, xi2, xi3, xi4 = 1, 1, 1, 1
println("(1, 0, 0, 0), ", round_brackets(n1, n2, n3, n4, l1_, l2_, l3_, l4_, xi1, xi2, xi3, xi4))

n1, n2, n3, n4 = 2, 2, 2, 2
l1_, l2_, l3_, l4_ = 1, 1, 1, 1
xi1, xi2, xi3, xi4 = 1, 1, 1, 1
println(round_brackets(n1, n2, n3, n4, l1_, l2_, l3_, l4_, xi1, xi2, xi3, xi4))

println("integral")
n1, n2, n3, n4 = 2, 2, 2, 2
l1, l2, l3, l4 = 1, 1, 1, 1
xi1, xi2, xi3, xi4 = 1, 1, 1, 1
integral, err = quadgk(
            r -> r^2*R(n1, l1, xi1, r)*R(n2, l2, xi2, r),
            0, Infinity.∞, rtol=1e-10
        )
println(integral)

integral, err = quadgk(
            r -> r^2*R(n3, l3, xi3, r)*R(n4, l4, xi4, r), 
            0, Infinity.∞, rtol=1e-10
        )
println(integral)

integral, err = quadgk(
            r -> r^2*R(n1, l1, xi1, r)*R(n4, l4, xi4, r),
            0, Infinity.∞, rtol=1e-10
        )
println(integral)

integral, err = quadgk(
            r -> r^2*R(n3, l3, xi3, r)*R(n2, l2, xi2, r), 
            0, Infinity.∞, rtol=1e-10
        )
println(integral)

l, m, θ, ϕ = 1, 0, 0, 0
println(Y(l, m, θ, ϕ))