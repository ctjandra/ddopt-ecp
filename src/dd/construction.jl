include("dd.jl")

"""
Construct a decision diagram.
"""
function construct_DD(n::Int, objective::Array, ordering::Ordering, problem::ProblemSpecs)

	# TODO
	# - This is not a top-down search, but a simple depth-first one
	# - Implement long arcs
	# - Implement mapping between problem variables and DD variables

	@assert(n == length(objective))
	dd = DecisionDiagram(n)

	# TODO For top-down, this should probably only be a single Dict
	layer_states = [Dict{Any, Int}() for i=1:n+1]

	# Add source node
	source = add_node!(dd, 1, problem.initial_state)
	set_objval!(dd, source, 0.0)
	layer_states[1][problem.initial_state] = source

	unexplored_nodes = [source]
	niters = 0

	while !isempty(unexplored_nodes)

		# Read next unexplored node
		node = pop!(unexplored_nodes)
		niters += 1
		layer = get_layer(dd, node)

		if layer > n
			continue # terminal node
		end

		state = get_state(dd, node)
		node_objval = get_value(dd, node, :objval)
		var = get_var(ordering, layer)

		println("Layer $(layer): exploring $(node)")

		# Create children
		for val in problem.domain_function(var)
			new_state = problem.transition_function(state, var, val)

			# Nothing to do if new state is the false node
			if new_state == nothing
				continue
			end

			next_layer = layer + 1
			states_to_nodes = layer_states[next_layer]

			new_value = node_objval + val * objective[var]

			# Equivalence check
			if haskey(states_to_nodes, new_state)
				# Equivalent node found: update it
				t = states_to_nodes[new_state]
				t_objval = get_objval(dd, t)
				if t_objval < new_value
					set_objval!(dd, t, new_value)
				end
			else
				# No equivalent node: create new child
				t = add_node!(dd, next_layer, new_state)
				set_objval!(dd, t, new_value)
				states_to_nodes[new_state] = t
			end

			# Insert new node into list of unexplored nodes
			push!(unexplored_nodes, t)

			add_arc!(dd, node, t, val)
			println("Added $val-arc from $node to $t")
		end
	end

	return dd
end
