include_dependency("ordering.jl")
include_dependency("core.jl")


"""Problem specifications required to construct decision diagram"""
struct ProblemSpecs{S<:AbstractFloat}
	initial_state::S                 # Initial state at the root
	transition_function::Function    # Transition function of DP formulation
	domain_range::Array{Range{Int}, 1}  # domain range for each variable
end


"""
Construct a decision diagram.

# Input
- `n`: Number of variables (arc layers).
- `problem`: Problem specification (initial state, transition function, domain function) based on DP formulation of the problem.
- `ordering`: Ordering of the variables in DD.
- `reduced_arc`: Determines whether parallel arcs are merged to their lower and upper bound label values
"""
function construct_DD(n::Int, problem::ProblemSpecs{S}; ordering::Ordering = NoOrdering(), reduced_arc::Bool = false)::DecisionDiagram where S<:AbstractFloat

	# TODO
	# - This is not a top-down search, but a simple depth-first one
	# - Implement long arcs
	# - Implement mapping between problem variables and DD variables
	# - Implement reduction

	dd::DecisionDiagram{S} = DecisionDiagram(n)

	# TODO For top-down, this should probably only be a single Dict
	layer_states::Array{Dict{S, Int}, 1} = [Dict{S, Int}() for i=1:n+1]

	# Add source node
	source::Int = add_node!(dd, 1, problem.initial_state)	#the index of the created node
	layer_states[1][problem.initial_state] = source

	unexplored_nodes::Array{Tuple{Int, Int}, 1} = Array{Tuple{Int, Int}, 1}(0)
	push!(unexplored_nodes, (1, source))
	niters::Int = 0

	while !isempty(unexplored_nodes)

		# Read next unexplored node
		node::Tuple{Int, Int} = pop!(unexplored_nodes)
		niters += 1
		layer::Int = node[1]
		id::Int = node[2]

		if layer > n
			continue # terminal node
		end

		state::S = get_node_state(dd, layer, id)
		var::Int = get_var(ordering, layer)

		# println("Layer $(layer): exploring $(node) with state $(state)")

		# Create children
		for val::Int in problem.domain_range[var]
			new_state::Union{S,Void} = problem.transition_function(state, var, val)

			# Nothing to do if new state is the false node
			if new_state == nothing
				continue
			end

			#new_state = floor(new_state/10)
			new_state = floor(new_state)
			#new_state = floor(10*new_state)/10		# rounds down the fractional state values: ONLY use to construct RELAXED DD when state values are fractional, e.g., for the pricing problem

			next_layer::Int = layer + 1
			states_to_nodes::Dict{S, Int} = layer_states[next_layer]

			# Equivalence check
			if haskey(states_to_nodes, new_state)
				# Equivalent node found: update it
				t::Int = states_to_nodes[new_state]		#the index of the node with that state value
			else
				# No equivalent node: create new child
				t = add_node!(dd, next_layer, new_state)
				states_to_nodes[new_state] = t
				# Insert new node into list of unexplored nodes
				push!(unexplored_nodes, (next_layer,t))
			end

			#if reduced_arc
			#	add_reduced_arc!(dd, node, t, val)
			#else
				add_arc!(dd, layer, id, next_layer, t, val)
			#end
			# println("Added $val-arc from $node to $t")
		end
	end

	return dd
end
