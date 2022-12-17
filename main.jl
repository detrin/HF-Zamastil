
using Optim
using ArgParse
include("./hamiltonian.jl")
include("./energies.jl")

s = ArgParseSettings()
@add_arg_table s begin
    "--mode", "-m"
        help = "select from 'referece' and 'optim' modes"
        arg_type = String
        required = true
    "--n12", "-n"
        help = "a value of n_12"
        arg_type = Int64
        required = true
end
parsed_args = parse_args(ARGS, s)
mode = parsed_args["mode"]
n_sum = parsed_args["n12"]
# test_hamiltonian()
# test_energies()

energies_ref = Dict(
    2 => -2.847656,
    3=> -2.847656,
    4=> -2.895444,
    5=> -2.897109,
    6=> -2.900714,
    7=> -2.901452,
    8=> -2.902341,
    9=> -2.902654,
    10=> -2.902975,
    11=> -2.903127,
)
etas_ref = Dict(
    2=> 32 / 27.,
    3=> 32 / 27.,
    4=> 0.971,
    5=> 0.940,
    6=> 0.796,
    7=> 0.760,
    8=> 0.682,
    9=> 0.648,
    10=> 0.595,
    11=> 0.566,
)

if mode == "reference"
    energy_ref = energies_ref[n_sum]
    eta_ref = etas_ref[n_sum]

    states = get_states_discrete(n_sum)
    xi1j, xi2j, xi1k, xi2k = 1, 1, 1, 1
    Z = 2
    η = eta_ref
    println("n_12 = ", n_sum)
    # println("calculating matrices ...")
    Ham_0, S_mat, W_mat = get_operators(states, xi1j, xi2j, xi1k, xi2k, Z, η)
    # println("solving generalized eigenvalue problem ...")
    A = (Ham_0 + W_mat - S_mat) * Z^2 / η^2
    results = eigen(A, S_mat)
    energies = results.values
    
    println("E0 = ", energies[1])
    println("η = ", η)
    println("ΔE0:", energies[1] - energy_ref)
end

function objective(args; n_sum=4)
    η = args[1]
    states = get_states_discrete(n_sum)
    xi1j, xi2j, xi1k, xi2k = 1., 1., 1., 1.
    Z = 2.
    Ham_0, S_mat, W_mat = get_operators(states, xi1j, xi2j, xi1k, xi2k, Z, η)
    A = (Ham_0 + W_mat - S_mat) * Z^2 / η^2
    results = eigen(A, S_mat)
    energies = results.values
    return energies[1]
end

if mode == "optim"
    println("n_12 = ", n_sum)
    objective_(args) = objective(args; n_sum=n_sum)
    results = optimize(objective_, [1.0], Newton())
    println("iterations = ", Optim.iterations(results))
    println("E0 = ", Optim.minimum(results))
    println("η = ", Optim.minimizer(results)[1])
end
