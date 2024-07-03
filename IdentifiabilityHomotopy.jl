using HomotopyContinuation
using Nemo
using DynamicPolynomials
using ProgressMeter

include("ExampleSystems.jl")
include("Utils.jl")

function run_all_factors(sys, vars, n_runs=100, by_factor=false, radii=[1e-1, 1e-2, 1e-3], direction_method="default")
	N = length(vars)
	total = 0

	fail = 0 
	same_solution = 0
	bad_solution = 0
	good = 0

	if direction_method == "Jacobian"
		sys_nemo, R, vars_nemo = sys_to_nemo(sys, vars)
		S = matrix_space(R, N, N)
	    J = S()

	    for i in 1:N
	        for j in 1:N
	            J[i, j] = derivative(sys_nemo[i], vars_nemo[j])
	        end
	    end
	end
	
	sys_nemo, R, vars_nemo = sys_to_nemo(sys, vars)
	factors = get_factors(sys_nemo, R, vars_nemo)

	if by_factor
		for (num, factor) in enumerate(factors)
			print("Factor ", num, ": ")
			if length(string(factor)) >= 20
				println("big")
			else
				println(factor)
			end
		end
		println()
		println("Fail\tSame\tBad\tGood\tTotal")
	else
		println("Total ", length(factors) * length(radii) * n_runs, " runs")
	end
	pbar = Progress(length(factors) * length(radii) * n_runs, desc="Processing")

	(zero_point, delta) = (0, 0)

	for factor in factors
		if by_factor
			total = 0

			fail = 0 
			same_solution = 0
			bad_solution = 0
			good = 0
		end

		for it in 1:n_runs
			success = false
			cnt = 0
			n_max_fails = 10
			
			while cnt < n_max_fails && !success
				cnt += 1
		        try
		        	if direction_method == "default"
		            	(zero_point, delta) = get_root_and_direction(factor, vars)
		            else
		            	(zero_point, delta) = get_root_and_direction_J(J, factor, vars)
		            end
		            success = true
		        catch e
		            # println("Caught an error: $e. Retrying...")
		        end
		    end

		   	if !success
	        	println("fail")
				total += length(radii)
	        	fail += length(radii)
		    	for radius in radii
					next!(pbar)
				end
	        	continue
	        end

		    for radius in radii
				next!(pbar)
				total += 1

				use_gradient = true

				solution_1 = [zero_point[i] + delta[i] * radius for i in 1:length(zero_point)]
				omega = (-1 + 0*im)^(2/3)

				vals_0 = [eq(vars => zero_point) for eq in sys]
				vals_1 = [eq(vars => solution_1) for eq in sys]
				vals_2 = vals_0 + (vals_1 - vals_0) * omega^1
				vals_3 = vals_0 + (vals_1 - vals_0) * omega^2

				sysv = [sys[i] - v[i] for i in 1:(N+1)]

				sysp = System(sysv[1:N]; parameters = v[1:N]);

				sys_with_vals = [sys[i] - vals_1[i] for i in 1:(N+1)]	

				vals_1 = vals_1[1:N]
				vals_2 = vals_2[1:N]
				vals_3 = vals_3[1:N]

				tracker_options = HomotopyContinuation.TrackerOptions(
					automatic_differentiation=3,
					extended_precision=true
				)

				res = HomotopyContinuation.solve(sysp, solution_1; start_parameters=vals_1, target_parameters=vals_2, compile=true, tracker_options=tracker_options)
				if length(solutions(res)) == 0
					fail += 1
					continue
				end
				solution_2 = solutions(res)[1]

				res = HomotopyContinuation.solve(sysp, solution_2; start_parameters=vals_2, target_parameters=vals_3, compile=true, tracker_options=tracker_options)
				if length(solutions(res)) == 0
					fail += 1
					continue
				end
				solution_3 = solutions(res)[1]

				res = HomotopyContinuation.solve(sysp, solution_3; start_parameters=vals_3, target_parameters=vals_1, compile=true, tracker_options=tracker_options)
				if length(solutions(res)) == 0
					fail += 1
					continue
				end
				solution_4 = solutions(res)[1]

				if norm(solution_1 .- solution_4) < 1e-8
					same_solution += 1
				elseif norm(sys_with_vals[end](vars => solution_4)) > 1e-6
					bad_solution += 1
				else
					good += 1
				end
			end
		end

		if by_factor
			println()
			println(fail, "\t", same_solution, "\t", bad_solution, "\t", good, "\t", total)
		end
	end

	if !by_factor
		println()
		println("Fail: ", fail, '/', total)
		println("Same solution: ", same_solution, '/', total)
		println("Bad solution: ", bad_solution, '/', total)
		println("Good: ", good, '/', total)
	end
end