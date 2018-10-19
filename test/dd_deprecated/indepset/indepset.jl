# DP formulation for independent set

# State is IntSet

"""Independent set transition function"""
function create_indepset_transition_function(g::Graph)
	return function indepset_transition_function(state::IntSet, variable::Int, value::Int)
		@assert(value == 0 || value == 1)
		new_state = IntSet(state)
		if value == 0
			delete!(new_state, variable)
		else # value == 1
			if !in(variable, state)
				return nothing # indicates false state
			end
			delete!(new_state, variable)
			setdiff!(new_state, out_neighbors(g, variable))
		end
		return new_state
	end
end

"""Independent set initial state"""
function create_indepset_initial_state(g::Graph)
	return IntSet(1:nv(g))
end
