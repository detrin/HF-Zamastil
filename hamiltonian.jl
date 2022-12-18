import CGcoefficient
using ClassicalOrthogonalPolynomials
using QuadGK
# using Infinity
using Polynomials
using LinearAlgebra
include("./fock_elements.jl")
include("./state_overlaps.jl")
include("./spherical_integrals.jl")


function get_operators(states, xi1j, xi2j, xi1k, xi2k, Z, η)
    Ham_0 = zeros(Float64, length(states), length(states))
    S_mat = zeros(Float64, length(states), length(states))
    W_mat = zeros(Float64, length(states), length(states))
    for st1_i in length(states):-1:1
        st1 = states[st1_i]
        n1j, n2j, lj, sign_j = st1
        for st2_i in length(states):-1:1
            st2 = states[st2_i]
            n1k, n2k, lk, sign_k = st2
            # println(st1, " ", st2)
            h0 = 0.
            h0 += fock_ham_one(n1j, n2j, n1k, n2k, lj, lk; η=η)
            h0 += fock_ham_one(n2j, n1j, n2k, n1k, lj, lk; η=η)
            h0 += sign_j*sign_k*fock_ham_one(n1j, n2j, n2k, n1k, lj, lk; η=η)
            h0 += sign_j*sign_k*fock_ham_one(n2j, n1j, n1k, n2k, lj, lk; η=η)
            Ham_0[st1_i, st2_i] = h0 / 4

            S = state_overlap_exact(
                n1j, n2j, n1k, n2k, 
                lj, lk, 
                xi1j, xi2j, xi1k, xi2k; 
                sign_j=sign_j, sign_k=sign_k
            )
            S_mat[st1_i, st2_i] = S

            W = fock_r12(
                n1j, n2j, n1k, n2k, 
                lj, lk, 
                xi1j, xi2j, xi1k, xi2k; 
                sign_j=sign_j, sign_k=sign_k, Z=Z, η=η
            )
            W_mat[st1_i, st2_i] = W
        end
    end
    return Ham_0, S_mat, W_mat
end

function test_hamiltonian()
    states = [
        [1, 1, 0, 1],
        [1, 2, 0, 1],
        [1, 3, 0, 1],
        [2, 2, 0, 1],
        [2, 2, 1, 1]
    ]
    xi1j, xi2j, xi1k, xi2k = 1, 1, 1, 1
    Z = 2
    η = 0.971
    Ham_0, S_mat, W_mat = get_operators(states, xi1j, xi2j, xi1k, xi2k, Z, η)
    println()
    println("Ham_0")
    println(display(Ham_0))
    Ham_0_ref = [ 0.058 -0.020506 0. 0. 0.;
    -0.020506 0.5435 -0.017759 -0.727613 0.;
    0. -0.017759 1.058 0. 0.;
    0. -0.727613 0. 4.116 0.;
    0. 0. 0. 0. 4.116]
    println(display(Ham_0-Ham_0_ref))
    println("S_mat")
    println(display(S_mat))
    S_mat_ref = [ 1.     -0.707107  0.        0.5       0.      ;
    -0.707107  1.25     -0.612372 -1.414214  0.      ;
    0.       -0.612372  1.5       0.866025  0.      ;
    0.5      -1.414214  0.866025  4.        0.      ;
    0.        0.        0.        0.        4.      ]
    println(display(S_mat-S_mat_ref))
    println()
    println("W_mat")
    println(display(W_mat))
    W_mat_ref = [ 0.303437 -0.107281 -0.052557  0.091031 -0.122633;
    -0.107281  0.227578 -0.037163 -0.203835  0.130072;
    -0.052557 -0.037163  0.20482   0.091975  0.007586;
    0.091031 -0.203835  0.091975  0.584117 -0.197088;
    -0.122633  0.130072  0.007586 -0.197088  0.842039]
    println(display(W_mat-W_mat_ref))

end

