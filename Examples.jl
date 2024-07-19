include("IdentifiabilityHomotopy.jl")

#====================================#
#                                    #
#  Compare 4 methods on one system   #
#                                    #
#====================================#

#=
(sys, vars) = ExampleSystems[5]

println("Homotopy continuation")
run_full_homotopy_continuation(sys, vars)
println()

println("Monodromy around random points")
run_random_monodromy(sys, vars, 100)
println()

println("Monodromy around zeros of Jacobian with gradient direction")
run_all_factors(sys, vars, 300, false, [1e-4], "default")
println()

println("Monodromy around zeros of Jacobian with Ju = 0")
run_all_factors(sys, vars, 300, false, [1e-4], "Jacobian")
println()
=#

#====================================#
#                                    #
# By factor statistics for a system  #
#                                    #
#====================================#

println("System 1")
run_all_factors(ExampleSystems[2]..., 30, true, [1e-2], "default")
println()

# println("System 5")
# run_all_factors(ExampleSystems[5]..., 300, true, [1e-4], "default")
# println()

# println("System 2")
# run_all_factors(ExampleSystems[2]..., 30, true, [1e-5], "Jacobian")
# println()

# (sys, vars) = ExampleSystems[5]

# N = length(vars)
# my_vals = [rand() for i in 1:N]

# sub_vals = [
#     eval(subs(eq, vars => my_vals)) for eq in sys
# ]
# sysv = [
#     sys[i] - v[i] for i in 1:(N+1)
# ]
# sysp = System(sysv[1:N]; parameters = v[1:N]);

# sys_with_vals = [
#     sys[i] - sub_vals[i] for i in 1:(N+1)
# ]

# total = 0
# fail = 0
# same_solution = 0
# bad_solution = 0
# good = 0

# mres = monodromy_solve(sysp, my_vals, sub_vals[1:N])

# for sol in solutions(mres)
# 	if norm(sys_with_vals[end](vars => sol)) < 1e-6
# 		println(sol)
# 	end
# end