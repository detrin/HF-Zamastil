

import CGcoefficient

using Memoize
using QuadGK
using Infinity

Base.factorial(x::Float64) = Float64(Base.factorial(Int64(x)))

# @memoize
function R_normalisation(n, l, η; a0=1.0)
    if n == 1 && l == 0
        return 1.
    end
    if n == 2 && l == 0
        return 2.
    end
    if n == 2 && l == 1
        return 1.3888888888888888
    end
    if n == 3 && l == 0
        return 3.0
    end
    if n == 3 && l == 1
        return 0.03038194444444443
    end
    if n == 3 && l == 2
        return 0.04666666666666665
    end
end

function R(n, l, η, r; a0=1.0)
    # Zamastil 2017 - 4.96
    #=
    if n == 1 && l == 0
        return 2 * sqrt(a0^3) * exp(-a0 * r)
    end
    if n == 2 && l == 0
        c0 = sqrt(2) * a0^(3 / 2) / 2
        c1 = -sqrt(2) * a0^(5 / 2) / 4
        return (c0 + c1 * r) * exp(-a0 * r / 2)
    end
    if n == 2 && l == 1
        return 2/sqrt(a0^5 / factorial(3)) * 2*r * exp(-a0 * r / 2)
    end
    if n == 3 && l == 2
        return 2*sqrt(a0^7 / factorial(5)) * (2*r)^2 * exp(-a0 * r / 3)
    end
    =#
    if n-1 == l && n == 1
        return 2*η / sqrt(factorial(2*l+1)) * (2*η*r)^l * exp(-η*r)
    end
    if n == 2 && l == 0
        return 2*η * sqrt(2) * (1 - η*r) * exp(-η*r)
    end
    # if n == 2 && l == 1
    #     return 2*η * sqrt(2) * (1 - η*r) * exp(-η*r)
    # end
    if n == 3 && l == 0
        return 2*η * sqrt(3) *
            (1 - 2 * η*r + 2 * (η*r)^2 / 3) *
            exp(-η*r)
    end
    # https://demonstrations.wolfram.com/HydrogenAtomRadialFunctions/
    # println(n, " ", l, " ", η, " ", r)
    #=
    if true
        N, err = quadgk(
            r -> r * (sqrt(factorial(n-1-l) / (2*n*factorial(n+l)^3)) *
            (2 / (n * a0))^(l + 3/2.) *
            η * (η*r)^l * exp(-(η*r)/(n*a0)) *
            laguerrel(n+l, 2*l+1, 2*(η*r)/(n*a0)) )^2,
            0, Infinity.∞, rtol=1e-8
        )
    end
    println(N)
    =#
    # println(n, l, η)
    N = R_normalisation(n, l, η; a0=1.0)
    # println(n, l, η, N)
    return sqrt(factorial(n-1-l) / (2*n*factorial(n+l)^3)) *
        (2 / (n * a0))^(l + 3/2.) *
        η * (η*r)^l * exp(-(η*r)/(n*a0)) *
        laguerrel(n+l, 2*l+1, 2*(η*r)/(n*a0)) # / sqrt(N) # / 2
end


function J_integral(args)
    eta, xi, a, b, l = args
    # println("J_integral: ", eta, " ", xi, " ", a, " ", b, " ", l)
    ret = 0.
    for q in 0:(a-l-1)
        # println(q)
        ret += factorial(b+ l + q) / factorial(q) * eta^(-(a-l-q)) * (eta+xi)^(-(b+l+q+1))
    end
    return factorial(a - l - 1) * ret
end 

function I_integral(args)
    eta, xi, a, b, l = args
    return J_integral([eta, xi, a, b, l]) + J_integral([xi, eta, b, a, l])
end

@memoize function square_brakcet_integral_root(args)
    l1, l2, l3, l4 = args
    return (
        1
        / sqrt(
            factorial(2 * l4 + 1)
            * factorial(2 * l3 + 1)
            * factorial(2 * l2 + 1)
            * factorial(2 * l1 + 1)
        )
        * 2^float(-l1 - l2 - l3 - l4 - 4)
        * factorial(l1 + l2 + l3 + l4 + 3)
    )
end

function square_brackets_equation(args)
    n1, n2, n3, n4, l1_, l2_, l3_, l4_, xi1, xi2, xi3, xi4, target = args
    l1 = l1_; l2 = l2_; l3 = l3_; l4 = l4_; 
    if target == :l4
        l4 += 1
    elseif target == :l3
        l3 += 1
        l4, l3 = l3, l4
        n4, n3 = n3, n4
        xi4, xi3 = xi3, xi4
    elseif target == :l2
        l2 += 1
        l4, l2 = l2, l4
        n4, n2 = n2, n4
        xi4, xi2 = xi2, xi4
    elseif target == :l1
        l1 += 1
        l4, l1 = l1, l4
        n4, n1 = n1, n4
        xi4, xi1 = xi1, xi4
    end

    c0 = -xi4*sqrt(n4^2 - l4^2)/l4 * (-1 + (l4 - l3 - l2 - l1 - 2) / (2*l4 + 1)) 
    c1 = xi3*n3/(l3 + 1) + xi2*n2/(l2 + 1) + xi1*n1/(l1 + 1) - xi4*n4/l4 + xi4*n4*(l4 - l3 - l2 - l1 - 2)/(l4*(l4+1))
    c2 = xi4*(l4 - l3 - l2 - l1 - 2)/(2*l4 + 1)*sqrt(n4^2 - (l4 + 1)^2)/(l4 + 1)
    c3 = xi3*sqrt(n3^2 - (l3 + 1)^2)/(l3 + 1)
    c4 = -xi2*sqrt(n2^2 - (l2 + 1)^2)/(l2 + 1)
    c5 = xi1*sqrt(n1^2 - (l1 + 1)^2)/(l1 + 1)

    ret = c1 * square_brackets([n1, n2, n3, n4, l1, l2, l3, l4, xi1, xi2, xi3, xi4]) +
        c2 * square_brackets([n1, n2, n3, n4, l1, l2, l3, l4+1, xi1, xi2, xi3, xi4]) + 
        c3 * square_brackets([n1, n2, n3, n4, l1, l2, l3+1, l4, xi1, xi2, xi3, xi4]) + 
        c4 * square_brackets([n1, n2, n3, n4, l1, l2+1, l3, l4, xi1, xi2, xi3, xi4]) + 
        c5 * square_brackets([n1, n2, n3, n4, l1+1, l2, l3, l4, xi1, xi2, xi3, xi4]) 
    return ret / c0
end

function square_brackets(args)
    n1, n2, n3, n4, l1, l2, l3, l4, xi1, xi2, xi3, xi4 = args
    rtol=1e-8
    if (n1-1 < l1) | (n2-1 < l2) | (n3-1 < l3) | (n4-1 < l4) 
        return 0.
    end

    if (0 > l1) | (0 > l2) | (0 > l3) | (0 > l4) 
        return 0.
    end

    if (n1-1 == l1) && (n2-1 == l2) && (n3-1 == l3) && (n4-1 == l4) 
        # return 2^(n1+n2+n3+n4)/sqrt(
        #     factorial(2n4-1)*factorial(2n3-1)*factorial(2n2-1)*factorial(2n1-1)
        # ) * SpecialFunctions.gamma(n1+n2+n3+n4)
        #=
        integral, err = quadgk(
            r -> r^3*R(n1, l1, xi1, r)*R(n2, l2, xi2, r)*R(n3, l3, xi3, r)*R(n4, l4, xi4, r), 
            0, Infinity.∞, rtol=rtol
        )
        =#
        integral = square_brakcet_integral_root([l1, l2, l3, l4])
        return integral
    end
    
    if n4-1 != l4
        ret = square_brackets_equation([n1, n2, n3, n4, l1, l2, l3, l4, xi1, xi2, xi3, xi4, :l4])
    end

    if n3-1 != l3
        ret = square_brackets_equation([n1, n2, n3, n4, l1, l2, l3, l4, xi1, xi2, xi3, xi4, :l3])
    end

    if n2-1 != l2
        ret = square_brackets_equation([n1, n2, n3, n4, l1, l2, l3, l4, xi1, xi2, xi3, xi4, :l2])
    end

    if n1-1 != l1
        ret = square_brackets_equation([n1, n2, n3, n4, l1, l2, l3, l4, xi1, xi2, xi3, xi4, :l1])
    end

    return ret
end



function round_brackets_1_equation(args)
    n1, n2, n3, n4, l1_, l2_, l3, l4, xi1, xi2, xi3, xi4, l, sign, target = args
    l1 = l1_; l2 = l2_; 
    l_plus = -l - 1
    l_minus = l
    l_plus_minus = sign == +1 ? l_plus : l_minus
    if target == :l2
        l2 += 1
    elseif target == :l1
        l1 += 1
        l2, l1 = l1, l2
        n2, n1 = n1, n2
        xi2, xi1 = xi1, xi2
    end

    c0 = -xi2*sqrt(n2^2 - l2^2)/l2 * (-1 + (l2 + 1 - l1 - 2 - l_plus_minus) / (2*l2 + 1)) 
    c1 = xi1*n1/(l1 + 1) - xi2*n2/l2 + xi2*n2*(l2 + 1 - l1 - l_plus_minus - 2)/(l2*(l2+1))
    c2 = xi2*(l2 + 1 - l1 - l_plus_minus - 2)/(2*l2 + 1)*sqrt(n2^2 - (l2 + 1)^2)/(l2 + 1)
    c3 = xi1*sqrt(n1^2 - (l1 + 1)^2)/(l1 + 1)
    # println("round_brackets_1_equation c0", c0)
    ret = -sign * square_brackets([n1, n2, n3, n4, l1, l2, l3, l4, xi1, xi2, xi3, xi4]) + 
        c1 * round_brackets([n1, n2, n3, n4, l1, l2, l3, l4, xi1, xi2, xi3, xi4, sign, l]) +
        c2 * round_brackets([n1, n2, n3, n4, l1, l2+1, l3, l4, xi1, xi2, xi3, xi4, sign, l]) + 
        c3 * round_brackets([n1, n2, n3, n4, l1+1, l2, l3, l4, xi1, xi2, xi3, xi4, sign, l]) 
    return ret / c0
end

function round_brackets_2_equation(args)
    n1, n2, n3, n4, l1, l2, l3_, l4_, xi1, xi2, xi3, xi4, l, sign, target = args
    l3 = l3_; l4 = l4_; 
    l_plus = -l - 1
    l_minus = l
    l_minus_plus = sign == +1 ? l_minus : l_plus
    if target == :l4
        l4 += 1
    elseif target == :l3
        l3 += 1
        l4, l3 = l3, l4
        n4, n3 = n3, n4
        xi4, xi3 = xi3, xi4
    end

    c0 = -xi4*sqrt(n4^2 - l4^2)/l4 * (-1 + (l4 + 1 - l3 - 2 - l_minus_plus) / (2*l4 + 1)) 
    c1 = xi3*n3/(l3 + 1) - xi4*n4/l4 + xi4*n4*(l4 + 1 - l3 - l_minus_plus - 2)/(l4*(l4+1))
    c2 = xi4*(l4 + 1 - l3 - l_minus_plus - 2)/(2*l4 + 1)*sqrt(n4^2 - (l4 + 1)^2)/(l4 + 1)
    c3 = xi3*sqrt(n3^2 - (l3 + 1)^2)/(l3 + 1)
    # println("round_brackets_2_equation c0", c0)
    ret = sign * square_brackets([n1, n2, n3, n4, l1, l2, l3, l4, xi1, xi2, xi3, xi4]) + 
        c1 * round_brackets([n1, n2, n3, n4, l1, l2, l3, l4, xi1, xi2, xi3, xi4, sign, l]) +
        c2 * round_brackets([n1, n2, n3, n4, l1, l2, l3, l4+1, xi1, xi2, xi3, xi4, sign, l]) + 
        c3 * round_brackets([n1, n2, n3, n4, l1, l2, l3+1, l4, xi1, xi2, xi3, xi4, sign, l]) 
    return ret / c0
end

# @memoize
function round_brackets(args)
    n1, n2, n3, n4, l1, l2, l3, l4, xi1, xi2, xi3, xi4, sign, l = args
    # println("round_brackets", n1, " ", n2, " ", n3, " ", n4, " | ", l1, " ", l2, " ", l3, " ", l4, " | ", xi1, " ", xi2, " ", xi3, " ", xi4, " | ", sign, " ", l)
    rtol=1e-8
    if (n1-1 < l1) | (n2-1 < l2) | (n3-1 < l3) | (n4-1 < l4) 
        return 0.
    end

    if (0 > l1) | (0 > l2) | (0 > l3) | (0 > l4) 
        return 0.
    end
    # println(n1, n2, n3, n4, "/", l1, l2, l3, l4, "/", xi1, xi2, xi3, xi4, "/", sign, " ", l)
    if (n1-1 == l1) && (n2-1 == l2) && (n3-1 == l3) && (n4-1 == l4) 
        # return 2^(n1+n2+n3+n4)/sqrt(
        #     factorial(2n4-1)*factorial(2n3-1)*factorial(2n2-1)*factorial(2n1-1)
        # ) * SpecialFunctions.gamma(n1+n2+n3+n4)
        if sign == +1
            #=
            integral1, err = quadgk(
                r1 -> begin
                    integral2, err = quadgk(
                        r2 -> r2^(2+l)*R(n2, l2, xi2, r2)*R(n1, l1, xi1, r2), 
                        0, r1, rtol=rtol
                    );
                    r1^(2-l-1)*R(n4, l4, xi4, r1)*R(n3, l3, xi3, r1)*integral2
                end,
                0, Infinity.∞, rtol=rtol
            )
            =#
            integral1 = J_integral([2., 2., n1 + n2, n3 + n4, l])
            integral1 *= 2 / sqrt(factorial(2 * (n1-1) + 1)) * (2 * xi1)^(n1-1)
            integral1 *= 2 / sqrt(factorial(2 * (n2-1) + 1)) * (2 * xi2)^(n2-1)
            integral1 *= 2 / sqrt(factorial(2 * (n3-1) + 1)) * (2 * xi3)^(n3-1)
            integral1 *= 2 / sqrt(factorial(2 * (n4-1) + 1)) * (2 * xi4)^(n4-1)
            return integral1
        elseif sign == -1
            #=
            integral1, err = quadgk(
                r1 -> begin
                    integral2, err = quadgk(
                        r2 -> r2^(2-l-1)*R(n2, l2, xi2, r2)*R(n1, l1, xi1, r2), 
                        r1, Infinity.∞, rtol=rtol
                    );
                    r1^(2+l)*R(n4, l4, xi4, r1)*R(n3, l3, xi3, r1)*integral2
                end,
                0, Infinity.∞, rtol=rtol
            )
            =#
            integral1 = J_integral([2., 2., n3 + n4, n1 + n2, l])
            integral1 *= 2 / sqrt(factorial(2 * (n1-1) + 1)) * (2 * xi1)^(n1-1)
            integral1 *= 2 / sqrt(factorial(2 * (n2-1) + 1)) * (2 * xi2)^(n2-1)
            integral1 *= 2 / sqrt(factorial(2 * (n3-1) + 1)) * (2 * xi3)^(n3-1)
            integral1 *= 2 / sqrt(factorial(2 * (n4-1) + 1)) * (2 * xi4)^(n4-1)
            return integral1
        end
    end

    if n4-1 != l4
        ret = round_brackets_2_equation([n1, n2, n3, n4, l1, l2, l3, l4, xi1, xi2, xi3, xi4, l, sign, :l4])
    end

    if n3-1 != l3
        ret = round_brackets_2_equation([n1, n2, n3, n4, l1, l2, l3, l4, xi1, xi2, xi3, xi4, l, sign, :l3])
    end

    if n2-1 != l2
        ret = round_brackets_1_equation([n1, n2, n3, n4, l1, l2, l3, l4, xi1, xi2, xi3, xi4, l, sign, :l2])
    end

    if n1-1 != l1
        ret = round_brackets_1_equation([n1, n2, n3, n4, l1, l2, l3, l4, xi1, xi2, xi3, xi4, l, sign, :l1])
    end
    # println("ret", ret)
    return ret
end

