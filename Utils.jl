using HomotopyContinuation
using Nemo
using DynamicPolynomials
using LinearAlgebra

function sys_to_nemo(sys, vars)
	R, vars_nemo = polynomial_ring(QQ, [string(v) for v in vars], internal_ordering=:lex)
	matching = Dict(string(x) => x for x in vars_nemo)

	sys_nemo = []

	for eq in sys
		eq_nemo = zero(R)
		for (monom, coef) in zip(DynamicPolynomials.monomials(eq), DynamicPolynomials.coefficients(eq))
			monom_nemo = coef
			for (x, p) in zip(DynamicPolynomials.variables(monom), DynamicPolynomials.exponents(monom))
				monom_nemo = monom_nemo * matching[string(x)]^p
			end
			eq_nemo += monom_nemo
		end
		push!(sys_nemo, eq_nemo)
	end

	sys_nemo, R, vars_nemo
end

function get_factors(sys, R, vars)
    N = length(vars)
    S = matrix_space(R, N, N)
    J = S()

    for i in 1:N
        for j in 1:N
            J[i, j] = derivative(sys[i], vars[j])
        end
    end

    D = det(J)

    cnt = 0

    factors = [f[1] for f in factor(D)]

    factors
end

function get_root_and_direction(factor, vars)
	N = length(vars)
    CC = AcbField(64)
    R_univar, t = polynomial_ring(CC, "t")
    # random_linear = [(1. + rand() / 100) * t + rand() for i in 1:N]
    random_linear = [(1. + rand()) * t + rand() for i in 1:N]

    all_roots = []
    all_roots_2 = []

    P = factor(random_linear...)
    # println("Total roots ", length(roots(P)))
    rt = roots(P)[1]
    root = [convert(ComplexF64, r(rt)) for r in random_linear]
    root_2 = [r(rt) for r in random_linear]

    dist_min = 1e9
    for rt in roots(P)[2:end]
    	other_root = [convert(ComplexF64, r(rt)) for r in random_linear]
    	dist_new = norm(root - other_root)
    	if dist_new < dist_min
    		dist_min = dist_new
    	end
    end
    # println(dist_min)

    Rc, varsc = polynomial_ring(CC, [string(v) for v in vars], internal_ordering=:lex)
    p = factor(varsc...)
    ps = derivative.(p, varsc)

    dir = [convert(ComplexF64, ps[i](root_2...)) for i in 1:N]
    dir = dir / norm(dir)

    (root, dir)
end

function get_radius(factor, vars, zero_point, dir)
	N = length(vars)
    CC = AcbField(64)
    R_univar, t = polynomial_ring(CC, "t")
    random_linear = [dir[i] * t + zero_point[i] for i in 1:N]

    P = factor(random_linear...)
    # println("Total roots ", length(roots(P)))
    radius = 10.

    try
	    for rt in roots(P)
	    	point_new = [convert(ComplexF64, r(rt)) for r in random_linear]
	    	radius_new = norm(point_new - zero_point)
	    	if 1e-6 < radius_new && radius_new < radius
	    		radius = radius_new
	    	end
	    end
	catch e
	    return radius
	end

    radius
end

function get_root_and_direction_J(J, factor, vars)
	N = length(vars)
	CC = AcbField(64)
	_, t = polynomial_ring(CC, "t")
	# random_linear = [i * t + rand() for i in 1:N]
    random_linear = [(1. + rand()) * t + rand() for i in 1:N]

	all_roots = []
	all_roots_2 = []

	P = factor(random_linear...)
	rt = roots(P)[1]
    root = [convert(ComplexF64, evaluate2(r, rt)[1]) for r in random_linear]
	root_2 = [r(rt) for r in random_linear]

	S = matrix_space(CC, N, N)
	I = S()

	_, vars_c = polynomial_ring(CC, [string(v) for v in vars], internal_ordering=:lex)

	for i in 1:N
	    for j in 1:N
	    	Pc = J[i, j](vars_c...)
	    	if typeof(Pc) == QQFieldElem
	    		I[i, j] = CC(Pc)
	    	else	
	        	I[i, j] = Pc(root_2...)
	        end
	    end
	end

	I = [convert(ComplexF64, a) for a in I]
	U, S, Vt = svd(I)
	dir = Vt[:, end]

    (root, dir)
end

function run_full_homotopy_continuation(sys, vars)
	N = length(vars)
	my_vals = [rand() for i in 1:N]
	sub_vals = [
	    convert(Float32, subs(eq, vars => my_vals)) for eq in sys
	]
	sys_with_vals = [
	    sys[i] - sub_vals[i] for i in 1:(N+1)
	]
	sysv = [
	    sys[i] - v[i] for i in 1:(N+1)
	]

	result = HomotopyContinuation.solve(sys_with_vals[1:N])

	total = length(result.path_results)
	# fail = nfailed(result)
	fail = total - length(solutions(result))
	same_solution = 1
	good = sum([norm(sys_with_vals[end](vars => new_solution)) < 1e-4 for new_solution in solutions(result)]) - 1
	bad_solution = total - fail - same_solution - good

	println()
	println("Fail: ", fail, '/', total)
	println("Same solution: ", same_solution, '/', total)
	println("Bad solution: ", bad_solution, '/', total)
	println("Good: ", good, '/', total)
end

function run_random_monodromy(sys, vars, n_runs=1000)
	N = length(vars)
	my_vals = [rand() for i in 1:N]

	sub_vals = [
	    eval(subs(eq, vars => my_vals)) for eq in sys
	]
	sysv = [
	    sys[i] - v[i] for i in 1:(N+1)
	]
	sysp = System(sysv[1:N]; parameters = v[1:N]);

	sys_with_vals = [
	    sys[i] - sub_vals[i] for i in 1:(N+1)
	]

	pbar = Progress(n_runs)

	total = 0
	fail = 0
	same_solution = 0
	bad_solution = 0
	good = 0

	for _ in 1:n_runs
	    mres = monodromy_solve(sysp, my_vals, sub_vals[1:N], target_solutions_count=2)

	    same_solution += length(mres.loops) + 1 - length(mres.results)
	    total += length(mres.loops)

	    for new_solution in solutions(mres)
	        if norm(my_vals .- new_solution) < 1e-6
	            # println("Same solution")
	            same_solution += 1
	        elseif norm(sys_with_vals[end](vars => new_solution)) > 1e-6
	            # println("Bad solution")
	            bad_solution += 1
	        else
	            # println("Good solution!")
	            println(new_solution)
	            println(my_vals)
	            good += 1
	            break
	        end
	    end
	    same_solution -= 1
	    # println()

	    next!(pbar)
	end

	println()
	println("Fail: ", fail, '/', total)
	println("Same solution: ", same_solution, '/', total)
	println("Bad solution: ", bad_solution, '/', total)
	println("Good: ", good, '/', total)
end