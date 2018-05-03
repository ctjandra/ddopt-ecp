using LightGraphs, MetaGraphs


"""Attach a value to a decision diagram node"""
set_node_value!(dd::DecisionDiagram, node::Int, key::Symbol, value::Any) = set_prop!(dd.graph, node, key, value)
get_node_value(dd::DecisionDiagram, node::Int, key::Symbol) = get_prop(dd.graph, node, key)

# Helper functions to obtain properties of node
set_node_state!(dd::DecisionDiagram, node::Int, value::Any) = set_prop!(dd.graph, node, :state, value)
get_node_state(dd::DecisionDiagram, node::Int) = get_prop(dd.graph, node, :state)
set_node_layer!(dd::DecisionDiagram, node::Int, value::Int) = set_prop!(dd.graph, node, :layer, value)
get_node_layer(dd::DecisionDiagram, node::Int) = get_prop(dd.graph, node, :layer)
set_node_idlayer!(dd::DecisionDiagram, node::Int, value::Int) = set_prop!(dd.graph, node, :idlayer, value)	#idlayer represents the position of the node in its layer
get_node_idlayer(dd::DecisionDiagram, node::Int) = get_prop(dd.graph, node, :idlayer)

# Non-essential properties that are common
set_node_obj_val!(dd::DecisionDiagram, node::Int, value::Real) = set_prop!(dd.graph, node, :obj_val, value)	#obj_val represents the longest path from source to the node
get_node_obj_val(dd::DecisionDiagram, node::Int) = get_prop(dd.graph, node, :obj_val)
set_node_obj_path!(dd::DecisionDiagram, node::Int, value::Real) = set_prop!(dd.graph, node, :obj_path, value)	#obj_path represents the parent node on the longest path to the node
get_node_obj_path(dd::DecisionDiagram, node::Int) = get_prop(dd.graph, node, :obj_path)
set_node_obj_label!(dd::DecisionDiagram, node::Int, value::Real) = set_prop!(dd.graph, node, :obj_label, value)	#obj_label represents the label of the longest path to the node
get_node_obj_label(dd::DecisionDiagram, node::Int) = get_prop(dd.graph, node, :obj_label)


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


"""Attach a value to a decision diagram arc"""
#method to operate on the tail and head nodes of the edge
set_arc_value!(dd::DecisionDiagram, node1::Int, node2::Int, key::Symbol, value::Any) = set_prop!(dd.graph, Edge(node1, node2), key, value)
get_arc_value(dd::DecisionDiagram, node1::Int, node2::Int, key::Symbol) = get_prop(dd.graph, Edge(node1, node2), key)
#method to operate on edge directly
set_arc_value!(dd::DecisionDiagram, arc::LightGraphs.SimpleGraphs.SimpleEdge{Int64}, key::Symbol, value::Any) = set_prop!(dd.graph, arc, key, value)
get_arc_value(dd::DecisionDiagram, arc::LightGraphs.SimpleGraphs.SimpleEdge{Int64}, key::Symbol) = get_prop(dd.graph, arc, key)


# Helper functions to obtain properties of arcs
#NOTE: arc labels are stored in arrays
set_arc_label!(dd::DecisionDiagram, node1::Int, node2::Int, value::Array{T} where T<:Number) = set_prop!(dd.graph, Edge(node1, node2), :label, value)
get_arc_label(dd::DecisionDiagram, node1::Int, node2::Int) = get_prop(dd.graph, Edge(node1, node2), :label)
set_arc_layer!(dd::DecisionDiagram, node1::Int, node2::Int, value::Int) = set_prop!(dd.graph, Edge(node1, node2), :layer, value)	#layer of an arc is computed wrt the head node (this definition accounts for long arcs)
get_arc_layer(dd::DecisionDiagram, node1::Int, node2::Int) = get_prop(dd.graph, Edge(node1, node2), :layer)
#method to operate on edge directly
set_arc_label!(dd::DecisionDiagram, arc::LightGraphs.SimpleGraphs.SimpleEdge{Int64}, value::Array{T} where T<:Number) = set_prop!(dd.graph, arc, :label, value)
get_arc_label(dd::DecisionDiagram, arc::LightGraphs.SimpleGraphs.SimpleEdge{Int64}) = get_prop(dd.graph, arc, :label)
set_arc_layer!(dd::DecisionDiagram, arc::LightGraphs.SimpleGraphs.SimpleEdge{Int64}, value::Int) = set_prop!(dd.graph, arc, :layer, value)	#layer of an arc is computed wrt the head node (this definition accounts for long arcs)
get_arc_layer(dd::DecisionDiagram, arc::LightGraphs.SimpleGraphs.SimpleEdge{Int64}) = get_prop(dd.graph, arc, :layer)


# Non-essential properties that are common
#NOTE: arc labels are stored in arrays, one weight for each arc label
set_arc_weight!(dd::DecisionDiagram, node1::Int, node2::Int, value::Array{T} where T<:Number) = set_prop!(dd.graph, Edge(node1, node2), :weight, value)
get_arc_weight(dd::DecisionDiagram, node1::Int, node2::Int) = get_prop(dd.graph, Edge(node1, node2), :weight)
#method to operate on edge directly
set_arc_weight!(dd::DecisionDiagram, arc::LightGraphs.SimpleGraphs.SimpleEdge{Int64}, value::Array{T} where T<:Number) = set_prop!(dd.graph, arc, :weight, value)
get_arc_weight(dd::DecisionDiagram, arc::LightGraphs.SimpleGraphs.SimpleEdge{Int64}) = get_prop(dd.graph, arc, :weight)



"""Add a new arc to the decision diagram between two nodes"""
function add_arc!(dd::DecisionDiagram, node1::Int, node2::Int, label::Number)
	@assert(get_node_layer(dd, node1) < get_node_layer(dd, node2))
	g = dd.graph
	e = Edge(node1, node2)
	# We use a list for labels because LightGraphs does not support parallel edges
	if !has_edge(g, e)
		add_edge!(g, e)
		set_prop!(g, e, :layer, get_node_layer(dd, node2) - 1)	#if the head node has layer k, the arc layer is k-1 (account for long arcs too)
		set_prop!(g, e, :label, [label])
	else
		push!(get_prop(g, e, :label), label)
	end
end

"""Binary domain for all variables for DD construction"""
domain_function_binary(var::Int) = [0, 1]
