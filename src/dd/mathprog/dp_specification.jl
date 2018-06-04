# DP formulation for additively separable constraints
# State is Array{Real}: right-hand sides of constraints

include_dependency("decomposition.jl")

"""
    create_constraint_transition_function(evals::SeparableFunctionEvaluation, doms::Array{Range}, ords::Ordering)

Compute transition function for additively separable constraints.

# Input:
- `evals`: An array of decomposed terms of a constraint.
- `doms`: Array of range representing domain of variables.
- `ords`: Function for ordering of the variables.
"""
function create_constraint_transition_function(evals::SeparableFunctionEvaluation, doms::Array{Range{Int}}, ords::Ordering)
	n = length(doms)	# number of variables
	min_ahead = Array{Float64}(n)		# For each variable, stores the minimum value that can be obtained starting from the next layer to the terminal
	max_ahead = Array{Float64}(n)		# For each variable, stores the maximum value that can be obtained starting from the next layer to the terminal

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
	return function (state::Float64, variable::Int, value::Int)
		@assert(value in doms[variable])		#ensures that the value belongs to the variable domain

		new_state::Float64 = state
		if isassigned(evals.separate_functions, variable)
			new_state += evals.separate_functions[variable](value)		#computes the new state value
		end

		#To determine feasibility, it is assumed that the constraint form is f(x) <= 0
		if new_state <= -max_ahead[variable]			#new state is feasible no matter what variable assignments are chosen from now on
			new_state = -max_ahead[variable]
		elseif new_state > -min_ahead[variable]			#new state is infeasible no matter what variable assignments are chosen
			return nothing
		end

		return new_state
	end
end


"""
Compute initial state for additively separable constraints.

# Input:
- `eval`: Decomposed terms of a constraint.

"""
function create_constraint_initial_state(eval::SeparableFunctionEvaluation)::Float64
    return eval.constant  	# constant value of the constraint as moved to lhs of <= constraints
end


"""
Compute domain range for additively separable constraints based on a given math prog model.

# Input:
- `m`: math prog model.

# Output:
- An array whose elements return the domain range of variables.
"""
function create_constraint_domain_function(m::Model)
	nvars = MathProgBase.numvar(m)
	ranges = Array{Range{Int}}(nvars)
	for v in 1:nvars
		v_obj = Variable(m, v)
		ranges[v] = ceil(Int, getlowerbound(v_obj)):floor(Int, getupperbound(v_obj))
	end
	return ranges
end
