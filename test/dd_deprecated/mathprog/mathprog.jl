# DP formulation for additively separable constraints
# State is Array{Real}: right-hand sides of constraints

include("mathprog_eval.jl")

"""^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Transition function for additively separable constraints

Input:
evals: An array of decomposed terms of a constraint
doms: Array of range representing domain of variables
ords: Function for ordering of the variables

Output: Transition function
"""
function create_constraint_transition_function(evals::SeparableFunctionEvaluation, doms::Array{Range}, ords::Ordering)
	@assert(length(doms) == length(ords)) # number of variables are consistent
	n = length(doms)	# number of variables
	min_ahead = Array{Real}(undef, n)		# For each variable, stores the minimum value that can be obtained starting from the next layer to the terminal
	max_ahead = Array{Real}(undef, n)		# For each variable, stores the maximum value that can be obtained starting from the next layer to the terminal

	#initialize values for the last layer
	min_ahead[get_var(ords,n)] = 0
	max_ahead[get_var(ords,n)] = 0
	for i in n-1:-1:1
		var = get_var(ords,i)		# gets the variable of current layer i
		var_next = get_var(ords,i+1)		# gets the variable in layer i+1
		if isassigned(evals.separate_functions, var_next)		#checks if the function is defined for the variable
			min_v, max_v = min_max_calculator(evals.separate_functions[var_next], doms[var_next])		#computes the min and max of the function over varibale domain
			min_ahead[var] = min_v + min_ahead[var_next]
			max_ahead[var] = max_v + max_ahead[var_next]
		else
			min_ahead[var] = min_ahead[var_next]	#if the function is not defined for the next layer, it does not affect the min max values
			max_ahead[var] = max_ahead[var_next]
		end
	end
	return function (state::Real, variable::Int, value::Int)
		@assert(value in doms[variable])		#ensures that the value belongs to the variable domain

		if isassigned(evals.separate_functions, variable)
			new_state = state + evals.separate_functions[variable](value)		#computes the new state value
		else
			new_state = state 		#if the function is not defined for the variable, the state value remains unchanged
		end

		#To determine feasibility, it is assumed that the constraint form is f(x) <= 0
		if new_state <= -max_ahead[variable]			#new state is feasible no matter what variable assignments are chosen from now on
			new_state = -max_ahead[variable]
		elseif new_state > -min_ahead[variable]			#new state is infeasible no matter what variable assignments are chosen
			new_state = nothing
		end

		return new_state
	end
end


"""^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Initial state for additively separable constraints

Input:
evals: An array of decomposed terms of a constraint

"""
function create_constraint_initial_state(evals::SeparableFunctionEvaluation)::Real
    return eval.constant  	# constant value of the constraint as moved to lhs of <= constraints
end


"""^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Domain function for additively separable constraints

Input:
m: math prog model

Output:
domain function that takes a variable index and returns the domain as range
"""
function create_constraint_domain_function(m::Model)
	nvars = MathProgBase.numvar(m)
	ranges = Array{Range}(undef, nvars)
	for v in 1:nvars
		v_obj = Variable(m, v)
		ranges[v] = ceil(Int, getlowerbound(v_obj)):floor(Int, getupperbound(v_obj))
	end
	return ranges
end
