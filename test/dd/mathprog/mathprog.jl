# DP formulation for additively separable constraints
# State is Array{Real}: right-hand sides of constraints

include("mathprog_eval.jl")

"""Transition function for additively separable constraints"""
function create_mathprog_transition_function(evals::Array{SeparableFunctionEvaluation})
	return function mathprog_transition_function(state::Array{Real}, variable::Int, value::Int)
		@assert(length(evals) == length(state)) # number of constraints are consistent
		new_state = copy(state)
		for (i, eval) in enumerate(evals)
			@assert(isassigned(eval.values, variable) && haskey(eval.values[variable], value))
			new_state[i] -= eval.values[variable][value]
		end

		# TODO: Assuming constraint is of <= form, if min(LHS) > RHS, state is infeasible

		return new_state
	end
end

"""Initial state for additively separable constraints"""
function create_mathprog_initial_state(evals::Array{SeparableFunctionEvaluation})::Array{Real}
    return [-eval.constant for eval in evals]  # negative constant as we move to RHS
end

"""Domain function for additively separable constraints"""
function create_mathprog_domain_function(m::Model)
	nvars = MathProgBase.numvar(m)
	ranges = Array{Range}(nvars)
	for v in 1:nvars
		v_obj = Variable(m, v)
		ranges[v] = ceil(Int, getlowerbound(v_obj)):floor(Int, getupperbound(v_obj))
	end
	return (var::Int -> ranges[var])
end
