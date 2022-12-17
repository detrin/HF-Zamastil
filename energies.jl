using LinearAlgebra
include("./hamiltonian.jl")
include("./states.jl")

function get_energies(states, xi1j, xi2j, xi1k, xi2k, Z, η)
    Ham_0, S_mat, W_mat = get_operators(states, xi1j, xi2j, xi1k, xi2k, Z, η)

    A = (Ham_0 + W_mat - S_mat) * Z^2 / η^2
    results = eigen(A, S_mat)
    return results.values
end

function test_energies()
    n_sum = 4
    states = get_states_discrete(n_sum)
    xi1j, xi2j, xi1k, xi2k = 1, 1, 1, 1
    Z = 2
    η = 0.971
    energies = get_energies(states, xi1j, xi2j, xi1k, xi2k, Z, η)
    println(energies)
    println("ΔE0:", energies[1] - (-2.895444))
end
