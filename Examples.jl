include("IdentifiabilityHomotopy.jl")

#====================================#
#                                    #
#  Compare 4 methods on one system   #
#                                    #
#====================================#

(sys, vars) = ExampleSystems[5]

println("Homotopy continuation")
run_full_homotopy_continuation(sys, vars)
println()

println("Monodromy around random points")
run_random_monodromy(sys, vars, 1000)
println()

println("Monodromy around zeros of Jacobian with gradient direction")
run_all_factors(sys, vars, 100)
println()

println("Monodromy around zeros of Jacobian with Ju = 0")
run_all_factors(sys, vars, 100)
println()

#====================================#
#                                    #
# By factor statistics for a system  #
#                                    #
#====================================#

println("System 1")
run_all_factors(ExampleSystems[1]..., 300, true, [1e-4], "Jacobian")
println()

println("System 5")
run_all_factors(ExampleSystems[5]..., 300, true, [1e-4], "default")
println()