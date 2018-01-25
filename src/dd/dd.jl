using LightGraphs, MetaGraphs

"""Problem specifications required to construct decision diagram"""
struct ProblemSpecs
	initial_state::Any               # Initial state at the root
	transition_function::Function    # Transition function of DP formulation
	domain_function::Function        # Function providing domain for each variable
end

"""Decision diagram structure"""
mutable struct DecisionDiagram
	graph::MetaGraph
    layers::Array{Array{Int}}
end

DecisionDiagram(nvars::Int) = DecisionDiagram(MetaGraph(), [[] for i=1:nvars+1])

"""Attach a value to a decision diagram node"""
set_value!(dd::DecisionDiagram, node::Int, key::Symbol, value::Any) = set_prop!(dd.graph, node, key, value)

# Helper functions to obtain properties of node
get_value(dd::DecisionDiagram, node::Int, key::Symbol) = get_prop(dd.graph, node, key)
get_state(dd::DecisionDiagram, node::Int) = get_prop(dd.graph, node, :state)
get_layer(dd::DecisionDiagram, node::Int) = get_prop(dd.graph, node, :layer)
get_idlayer(dd::DecisionDiagram, node::Int) = get_prop(dd.graph, node, :idlayer)

# Non-essential properties that are common
set_objval!(dd::DecisionDiagram, node::Int, value::Real) = set_prop!(dd.graph, node, :objval, value)
get_objval(dd::DecisionDiagram, node::Int) = get_prop(dd.graph, node, :objval)

"""Add a new node to the decision diagram with a given state and return the node's id"""
function add_node!(dd::DecisionDiagram, layer::Int, state::Any)
	g = dd.graph
	add_vertex!(g)
	id::Int = nv(g)
	push!(dd.layers[layer], id)
	# state, layer, and idlayer are essential properties; others are optional and can be set with set_value
	set_props!(g, id, Dict(:state => state, :layer => layer, :idlayer => length(dd.layers[layer])))
	return id
end

"""Add a new arc to the decision diagram between two nodes"""
function add_arc!(dd::DecisionDiagram, node1::Int, node2::Int, value::Int)
	@assert(get_layer(dd, node1) < get_layer(dd, node2))
	g = dd.graph
	e = Edge(node1, node2)
	# We use a list for values because LightGraphs does not support parallel edges
	if !has_edge(g, e)
		add_edge!(g, e)
		set_prop!(g, e, :value, [value])
	else
		push!(get_prop(g, e, :value), value)
	end
end

"""Binary domain for all variables for DD construction"""
domain_function_binary(var::Int) = [0, 1]
