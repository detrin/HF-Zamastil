
using QuadGK
import Infinity
include("./radial_integrals.jl")

function state_overlap_old(n1j::Int64, n2j::Int64, n1k::Int64, n2k::Int64, lj::Int64, lk::Int64, xi1j, xi2j, xi1k, xi2k; sign_j=+1, sign_k=+1, rtol=1e-8)
    println("function state_overlap")
    sign = sign_j * sign_k
    ret = 0.
    l = abs(lj - lk)
    while l <= lj + lk
        println(lj, " ", lk, " ", l)
        spherical_part = sqrt((2*lj + 1)*(2*lk + 1)) / (2*l + 1) * (-1)^(lj + lk) * 
            float(CGcoefficient.CG(lj, lk, l, 0, 0, 0))^2
        integral11, err = quadgk(
            r1 -> r1^2*R(n1j, lj, xi1j, r1)*R(n1k, lk, xi1k, r1),
            0, Infinity.∞, rtol=rtol
            )
        integral12, err = quadgk(
            r2 -> r2^2*R(n2j, lj, xi2j, r2)*R(n2k, lk, xi2k, r2),
            0, Infinity.∞, rtol=rtol
            )
        integral21, err = quadgk(
            r1 -> r1^2*R(n1j, lj, xi1j, r1)*R(n2k, lk, xi2k, r1),
            0, Infinity.∞, rtol=rtol
            )
        integral22, err = quadgk(
            r2 -> r2^2*R(n2j, lj, xi2j, r2)*R(n1k, lk, xi1k, r2),
            0, Infinity.∞, rtol=rtol
            )
        #=
        println("state_overlap")
        println(float(CGcoefficient.CG(lj, lk, l, 0, 0, 0))^2)
        println(sqrt((2*lj + 1)*(2*lk + 1)) / (2*l + 1) * (-1)^(lj + lk))
        println("spherical_part", spherical_part)
        println(integral11)
        println(integral12)
        println(integral21)
        println(integral22)
        
        =#
        ret += spherical_part * (integral11*integral12 + sign*integral21*integral22)
        l += 2
    end
    return ret / 2
end

function state_overlap_old2(n1j::Int64, n2j::Int64, n1k::Int64, n2k::Int64, lj::Int64, lk::Int64, xi1j, xi2j, xi1k, xi2k; sign_j=+1, sign_k=+1, rtol=1e-8)
    println("function state_overlap")
    sign = sign_j * sign_k
    ret = 0.
    # (n1j, lj), (n2j, lj) -> (n1k, lk), (n2k, lk)
    for mj in -lj:lj
        for mk in -lk:lk
            if n1j == n1k && n2j == n2k && lj == lk && mj == mk
                ret += float(CGcoefficient.CG(lj, lj, 0, mj, -mj, 0)*CGcoefficient.CG(lk, lk, 0, mk, -mk, 0))
            end
        end
    end
    # (n2j, lj), (n1j, lj) -> (n1k, lk), (n2k, lk) 
    for mj in -lj:lj
        for mk in -lk:lk
            if n2j == n1k && n1j == n2k && lj == lk && mj == mk
                ret += sign*float(CGcoefficient.CG(lj, lj, 0, mj, -mj, 0)*CGcoefficient.CG(lk, lk, 0, mk, -mk, 0))
            end
        end
    end
    # (n1j, lj), (n2j, lj) -> (n2k, lk), (n1k, lk) 
    for mj in -lj:lj
        for mk in -lk:lk
            if n1j == n2k && n2j == n1k && lj == lk && mj == mk
                ret += sign*float(CGcoefficient.CG(lj, lj, 0, mj, -mj, 0)*CGcoefficient.CG(lk, lk, 0, mk, -mk, 0))
            end
        end
    end
    # (n2j, lj), (n1j, lj) -> (n2k, lk), (n1k, lk)
    for mj in -lj:lj
        for mk in -lk:lk
            if n2j == n1k && n2j == n1k && lj == lk && mj == mk
                ret += float(CGcoefficient.CG(lj, lj, 0, mj, -mj, 0)*CGcoefficient.CG(lk, lk, 0, mk, -mk, 0))
            end
        end
    end
    return ret 
end

function state_overlap(n1j::Int64, n2j::Int64, n1k::Int64, n2k::Int64, lj::Int64, lk::Int64, xi1j, xi2j, xi1k, xi2k; sign_j=+1, sign_k=+1, rtol=1e-8)
    println("function state_overlap")
    sign = sign_j * sign_k
    ret = 0.
    # (n1j, lj), (n2j, lj) -> (n1k, lk), (n2k, lk)
    for mj in -lj:lj
        for mk in -lk:lk
            if n1j == n1k && n2j == n2k && lj == lk && mj == mk
                ret += (-1)^(mj+lj) / sqrt(2*lj + 1) * (-1)^(mk+lk) / sqrt(2*lk + 1)
            end
        end
    end
    # (n2j, lj), (n1j, lj) -> (n1k, lk), (n2k, lk) 
    for mj in -lj:lj
        for mk in -lk:lk
            if n2j == n1k && n1j == n2k && lj == lk && mj == mk
                ret += sign*(-1)^(mj+lj) / sqrt(2*lj + 1) * (-1)^(mk+lk) / sqrt(2*lk + 1)
            end
        end
    end
    # (n1j, lj), (n2j, lj) -> (n2k, lk), (n1k, lk) 
    for mj in -lj:lj
        for mk in -lk:lk
            if n1j == n2k && n2j == n1k && lj == lk && mj == mk
                ret += sign*(-1)^(mj+lj) / sqrt(2*lj + 1) * (-1)^(mk+lk) / sqrt(2*lk + 1)
            end
        end
    end
    # (n2j, lj), (n1j, lj) -> (n2k, lk), (n1k, lk)
    for mj in -lj:lj
        for mk in -lk:lk
            if n2j == n2k && n1j == n1k && lj == lk && mj == mk
                ret += (-1)^(mj+lj) / sqrt(2*lj + 1) * (-1)^(mk+lk) / sqrt(2*lk + 1)
            end
        end
    end
    return ret 
end