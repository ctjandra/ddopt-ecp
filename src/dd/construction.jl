include_dependency("ordering.jl")

using LightGraphs, MetaGraphs



"""Decision diagram structure"""
mutable struct DecisionDiagram
	graph::MetaDiGraph				#represents the directed graph of the DD
    layers::Array{Array{Int}}		#stores index of nodes at each node layer
end

DecisionDiagram(nvars::Int) = DecisionDiagram(MetaDiGraph(), [[] for i=1:nvars+1])

"""Return the number of node layers in the decision diagram"""
nlayers(dd::DecisionDiagram) = size(layers)

"""Return the number of variables (or number of arc layers) in the decision diagram"""
nvars(dd::DecisionDiagram) = size(layers) - 1


"""Problem specifications required to construct decision diagram"""
struct ProblemSpecs
	initial_state::Any               # Initial state at the root
	transition_function::Function    # Transition function of DP formulation
	domain_range::Array{Range}       # domain range for each variable
end



"""
Construct a decision diagram.

# Input
- `n`: Number of variables (arc layers).
- `problem`: Problem specification (initial state, transition function, domain function) based on DP formulation of the problem.
- `ordering`: Ordering of the variables in DD.
"""
function construct_DD(n::Int, problem::ProblemSpecs; ordering::Ordering = NoOrdering())

	# TODO
	# - This is not a top-down search, but a simple depth-first one
	# - Implement long arcs
	# - Implement mapping between problem variables and DD variables
	# - Implement reduction

	dd = DecisionDiagram(n)

	# TODO For top-down, this should probably only be a single Dict
	layer_states = [Dict{Any, Int}() for i=1:n+1]

	# Add source node
	source = add_node!(dd, 1, problem.initial_state)	#the index of the created node
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

		# println("Layer $(layer): exploring $(node) with state $(state)")

		# Create children
		for val in problem.domain_range[var]
			new_state = problem.transition_function(state, var, val)

			# Nothing to do if new state is the false node
			if new_state == nothing
				continue
			end

			#new_state = floor(new_state)
			new_state = floor(10*new_state)/10		# rounds down the fractional state values: ONLY use to construct RELAXED DD when state values are fractional, e.g., for the pricing problem

			next_layer = layer + 1
			states_to_nodes = layer_states[next_layer]

			# Equivalence check
			if haskey(states_to_nodes, new_state)
				# Equivalent node found: update it
				t = states_to_nodes[new_state]		#the index of the node with that state
			else
				# No equivalent node: create new child
				t = add_node!(dd, next_layer, new_state)
				states_to_nodes[new_state] = t
				# Insert new node into list of unexplored nodes
				push!(unexplored_nodes, t)
			end

			add_arc!(dd, node, t, val)
			# println("Added $val-arc from $node to $t")
		end
	end

	return dd
end
