
import CGcoefficient
include("./utils.jl")
include("./radial_integrals.jl")


function fock_r12(n1j::Int64, n2j::Int64, n1k::Int64, n2k::Int64, lj::Int64, lk::Int64, xi1j, xi2j, xi1k, xi2k; sign_j=+1, sign_k=+1, η=1.0, Z=2)
    sign = sign_j * sign_k
    l = abs(lj - lk)
    C_part = 0.
    V_part = 0.
    while l <= lj + lk
        
        spherical_part =  sqrt((2*lj + 1)*(2*lk + 1)) / (2*l + 1) * (-1)^(lj + lk) * 
            float(CGcoefficient.CG(lj, lk, l, 0, 0, 0))^2 / Z * η

        C_part += spherical_part * (
            round_brackets([n1j, n1k, n2j, n2k, lj, lk, lj, lk, xi1j, xi1k, xi2j, xi2k, +1, l])
            + round_brackets([n1j, n1k, n2j, n2k, lj, lk, lj, lk, xi1j, xi1k, xi2j, xi2k, -1, l])
        ) 

        V_part += spherical_part * (
            round_brackets([n1j, n2k, n2j, n1k, lj, lk, lj, lk, xi1j, xi2k, xi2j, xi1k, +1, l])
            + round_brackets([n1j, n2k, n2j, n1k, lj, lk, lj, lk, xi1j, xi2k, xi2j, xi1k, -1, l])
        ) 
        l += 2
    end
    
    return (C_part + sign*V_part) / 2 
end

function fock_ham_one(n1j::Int64, n2j::Int64, n1k::Int64, n2k::Int64, lj::Int64, lk::Int64; η=1.0::Float64)
    if lk != lj 
        return 0.
    end
    h1 = δ(n1k, n1j) * (n1j - η) * (
        δ(n2k, n2j) * n2j - 
        sqrt((n2j - lj) * (n2j + lj + 1)) / 2 * δ(n2k, n2j+1) - 
        sqrt((n2j + lj) * (n2j - lj - 1)) / 2 * δ(n2k, n2j-1) 
    )
    h2 = δ(n2k, n2j) * (n2j - η) * (
        δ(n1k, n1j) * n1j - 
        sqrt((n1j - lj) * (n1j + lj + 1)) / 2 * δ(n1k, n1j+1) - 
        sqrt((n1j + lj) * (n1j - lj - 1)) / 2 * δ(n1k, n1j-1) 
    )
    # println("one_particle_integral ", n1j, " ", n2j, " ", lj, " | ", n1k, " ", n2k, " ", lk)
    # println("fock_ham_one h1 ", h1)
    # println("fock_ham_one h2 ", h2)
    # println("fock_ham_one ", h1 + h2)
    return h1 + h2
end