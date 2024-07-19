using HomotopyContinuation
using DynamicPolynomials

function extend_system(sys_start, vars, x, x_prime)
    N = length(vars)
    sys = sys_start

    while length(sys) < N + 1
        T = sys[end]
        T = sum(differentiate.(T, x) .* x_prime)
        push!(sys, T)
    end

    sys
end

ExampleSystems = []


#====================================#
#                                    #
#             System 1               #
#                                    #
#====================================#

@polyvar x1 x2 v[1:3]

vars = [x1, x2]

sys = [
	x1^2 + x2^2,
	x1^3 + x2^3,
	x1^4 + x2^4,
]
push!(ExampleSystems, (sys, vars))


#====================================#
#                                    #
#             System 2               #
#                                    #
#====================================#

@polyvar E I S K beta epsilon mu gamma v[1:9]

vars = [E, I, S, K, beta, epsilon, mu, gamma]

N = length(vars)

E_prime = 3 * S * I * beta - E * epsilon - E * mu
I_prime = E * epsilon - I * gamma - I * mu
S_prime = -3 * S * I * beta - S * mu + 2

sys = extend_system([I * K, I_prime * K], vars, (E, I, S), (E_prime, I_prime, S_prime))
push!(ExampleSystems, (sys, vars))


#====================================#
#                                    #
#             System 3               #
#                                    #
#====================================#

@polyvar x1 x2 x3 p1 p3 p4 u v[1:9]

vars = [x1, x2, x3, p1, p3, p4, u]

x1_prime = -x1*p1 + u
x2_prime = -x2*p3 + u*p4
x3_prime = x1*u*p4 + x2*u - x3*p1 - x3*p3
u_prime  = u^2
sys = extend_system([x3+0], vars, (x1, x2, x3, u), (x1_prime, x2_prime, x3_prime, u_prime))
push!(ExampleSystems, (sys, vars))


#====================================#
#                                    #
#             System 4               #
#                                    #
#====================================#

@polyvar x1 x2 x3 k01 k21 k31 k12 k13 v[1:9]

vars = [x1, x2, x3, k01, k21, k31, k12, k13]

x1_prime = -x1*k01 - x1*k21 - x1*k31 + x2*k12 + x3*k13
x2_prime = x1*k21 - x2*k12
x3_prime = x1*k31 - x3*k13
sys = extend_system([x1+0], vars, (x1, x2, x3), (x1_prime, x2_prime, x3_prime))
push!(ExampleSystems, (sys, vars))


#====================================#
#                                    #
#             System 5               #
#                                    #
#====================================#

@polyvar mu1 mu2 s1 s2 a v[1:9]

vars = [mu1, mu2, s1, s2, a]

sys = [
    a*mu1 + (1 − a)*mu2,
    a*(mu1^2 + s1) + (1 − a)*(mu2^2 + s2),
    a*(mu1^3 + 3*mu1*s1) + (1 − a)*(mu2^3 + 3*mu2*s2),
    a*(mu1^4 + 6*mu1^2*s1 + 3*s1^2) + (1 − a)*(mu2^4 + 6*mu2^2*s2 + 3*s2^2),
    a*(mu1^5 + 10*mu1^3*s1 + 15*mu1*s1^2) + (1 − a)*(mu2^5 + 10*mu2^3*s2 + 15*mu2*s2^2),
    a*(mu1^6 + 15*mu1^4*s1 + 45*mu1^2*s1^2 + 15*s1^3) + (1 − a)*(mu2^6 + 15*mu2^4*s2 + 45*mu2^2*s2^2 + 15*s2^3),
]
push!(ExampleSystems, (sys, vars))


#====================================#
#                                    #
#             System 6               #
#                                    #
#====================================#

@polyvar x v[1:9]

vars = [x]

sys = [
    x^6,
    x^11 + x^12
]
push!(ExampleSystems, (sys, vars))