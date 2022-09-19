
import CGcoefficient
include("./utils.jl")
include("./radial_integrals.jl")


function fock_r12(n1j::Int64, n2j::Int64, n1k::Int64, n2k::Int64, lj::Int64, lk::Int64, xi1j, xi2j, xi1k, xi2k; sign_j=+1, sign_k=+1)
    sign = sign_j * sign_k
    ret = 0.
    l = abs(lj - lk)
    while l <= lj + lk
        
        spherical_part = sqrt((2*lj + 1)*(2*lk + 1)) / (2*l + 1) * (-1)^(lj + lk) * 
            float(CGcoefficient.CG(lj, lk, l, 0, 0, 0))^2
        
        ret += spherical_part 
        l += 2
    end
    round_brackets_1_plus = round_brackets(n1j, n1k, n2j, n2k, lj, lk, lj, lk, xi1j, xi1k, xi2j, xi2k; sign=+1, rtol=1e-8)
    round_brackets_1_minus = round_brackets(n1j, n1k, n2j, n2k, lj, lk, lj, lk, xi1j, xi1k, xi2j, xi2k; sign=-1, rtol=1e-8)
    round_brackets_2_plus = round_brackets(n1j, n2k, n2j, n1k, lj, lk, lj, lk, xi1j, xi2k, xi2j, xi1k; sign=+1, rtol=1e-8)
    round_brackets_2_minus = round_brackets(n1j, n2k, n2j, n1k, lj, lk, lj, lk, xi1j, xi2k, xi2j, xi1k; sign=-1, rtol=1e-8)
    # println(round_brackets_1_plus)
    # println(round_brackets_1_minus)
    # println(round_brackets_2_plus)
    # println(round_brackets_2_minus)
    radial_part = round_brackets_1_plus + round_brackets_1_minus +
        sign*round_brackets_2_plus + sign*round_brackets_2_minus
    return ret * radial_part / 2
end

function fock_ham_one(n1j::Int64, n2j::Int64, n1k::Int64, n2k::Int64, lj::Int64, lk::Int64; η=0.9::Float64)
    if lk != lj 
        return 0.
    end
    if n1k != n1j 
        return 0.
    end
    return (n1j - η) * (
        δ(n2k, n2j) * n2j - 
        sqrt((n2j - lj) * (n2j + lj + 1)) / 2 * δ(n2k, n2j+1) - 
        sqrt((n2j + lj) * (n2j - lj - 1)) / 2 * δ(n2k, n2j-1) 
    )
end