include("dd.jl")

"""
Construct a decision diagram.
"""
function construct_DD(n::Int, problem::ProblemSpecs; ordering::Ordering = NoOrdering())

	# TODO
	# - This is not a top-down search, but a simple depth-first one
	# - Implement long arcs
	# - Implement mapping between problem variables and DD variables
	# - Implement reduction

	@assert(n == length(objective))
	dd = DecisionDiagram(n)

	# TODO For top-down, this should probably only be a single Dict
	layer_states = [Dict{Any, Int}() for i=1:n+1]

	# Add source node
	source = add_node!(dd, 1, problem.initial_state)
	layer_states[1][problem.initial_state] = source

	unexplored_nodes = [source]
	niters = 0

	while !isempty(unexplored_nodes)

		# Read next unexplored node
		node = pop!(unexplored_nodes)
		niters += 1
		layer = get_node_layer(dd, node)

		if layer > n
			continue # terminal node
		end

		state = get_node_state(dd, node)
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


			# Equivalence check
			if haskey(states_to_nodes, new_state)
				# Equivalent node found: update it
				t = states_to_nodes[new_state]
			else
				# No equivalent node: create new child
				t = add_node!(dd, next_layer, new_state)
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
