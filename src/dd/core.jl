# Core structure and functions for decision diagrams

# All higher level functions should use these functions instead of accessing
# the internal structure of a decision diagram directly

# Functions in this file must be careful to keep graph and layers synchronized

using LightGraphs
using ResumableFunctions
import LightGraphs.SimpleGraphs.SimpleEdge

"""Decision diagram structure"""
mutable struct DecisionDiagram{S}
	graph::SimpleDiGraph{Int}            # represents the directed graph of the DD
    layers::Array{Array{Int,1},1}        # stores index of nodes at each node layer

	# Implementation ensures that node_layers, node_layerids, states always have size nv(graph)

	# Node properties
	node_layers::Array{Int,1}            # layer of node (i in layers[node_layers[i]])
	node_layerids::Array{Int,1}          # layer id of node (i == layers[node_layers[i]][node_layerids[i]])
	states::Array{S,1}                   # node states

	# Arc properties
	arc_labels::Dict{SimpleEdge{Int},Array{Int,1}}  # labels of arcs (due to no parallel arcs in LightGraphs, store list of labels)
end


DecisionDiagram(nvars::Int) = DecisionDiagram(SimpleDiGraph{Int}(), [Array{Int,1}(0) for i=1:nvars+1],
	Array{Int,1}(0), Array{Int,1}(0), [], Dict{SimpleEdge{Int},Array{Int,1}}())

"""Return the number of node layers in the decision diagram"""
nlayers(dd::DecisionDiagram) = size(dd.layers)

"""Return the number of variables (or number of arc layers) in the decision diagram"""
nvars(dd::DecisionDiagram) = size(dd.layers) - 1

nnodes(dd::DecisionDiagram) = nv(dd.graph)
narcs(dd::DecisionDiagram) = ne(dd.graph)

root(dd::DecisionDiagram) = dd.layers[1][1]
terminal(dd::DecisionDiagram) = dd.layers[end][1]


"""Set state of an existing node"""
function set_node_state!(dd::DecisionDiagram, node::Int, state::S) where S
	dd.states[node] = state
end

"""Get state of a node"""
function get_node_state(dd::DecisionDiagram, node::Int)
	return dd.states[node]
end

"""Get layer of a node"""
function get_node_layer(dd::DecisionDiagram, node::Int)::Int
	return dd.node_layers[node]
end

"""Get layer id of a node"""
function get_node_layerid(dd::DecisionDiagram, node::Int)::Int
	return dd.node_layerids[node]
end

"""Add a new node to the decision diagram with a given state and return the node's id"""
function add_node!(dd::DecisionDiagram, layer::Int, state::S) where S
	g = dd.graph
	add_vertex!(g)
	id::Int = nv(g)
	push!(dd.layers[layer], id)
	push!(dd.states, state)
	push!(dd.node_layers, layer)
	push!(dd.node_layerids, length(dd.layers[layer]))
	@assert(length(dd.states) == id)
	@assert(length(dd.node_layers) == id)
	@assert(length(dd.node_layerids) == id)
	return id
end

# TODO: rem_node! is untested
"""Remove a node from the decision diagram"""
function rem_node!(dd::DecisionDiagram, node::Int)
	v in vertices(g) || return false

	# Remove arc labels
	for outnode in outneighbors(dd.graph, node)
		delete!(dd.arc_labels, Edge(node, outnode))
	end
	for innode in inneighbors(dd.graph, node)
		delete!(dd.arc_labels, Edge(innode, node))
	end

	# Move end node in the layer to node's position within the layer
	layer = dd.node_layers[node]
	layerid = dd.node_layerids[node]
	dd.layers[layer][layerid] = pop!(dd.layers[layer])
	dd.layerids[dd.layers[layer][layerid]] = layerid

	# Renumber last node to node according to LightGraphs implementation
	dd.states[node] = pop!(dd.states)
	dd.node_layers[node] = pop!(dd.node_layers)
	dd.node_layerids[node] = pop!(dd.node_layerids)
	dd.layers[dd.node_layers[node]][dd.node_layerids[node]] = node

	# Remove vertex from graph structure
	removed = rem_vertex!(dd.graph, node)
	@assert(removed)
end

"""Remove a node from the decision diagram given its layer and layerid"""
function rem_node!(dd::DecisionDiagram, layer::Int, layerid::Int)
	rem_node!(dd, dd.layers[layer][layerid])
end

"""Add a new arc to the decision diagram between two nodes"""
function add_arc!(dd::DecisionDiagram, node1::Int, node2::Int, label::Int)
	@assert(dd.node_layers[node1] < dd.node_layers[node2])
	g = dd.graph
	e = Edge(node1, node2)
	# We use a list for labels because LightGraphs does not support parallel edges
	if !has_edge(g, e)
		add_edge!(g, e)
		dd.arc_labels[e] = [label]
	else
		push!(dd.arc_labels[e], label)
	end
end

# TODO rem_arc! is untested
"""Remove a labeled arc from a decision diagram"""
function rem_arc!(dd::DecisionDiagram, node1::Int, node2::Int, label::Int)
	e = Edge(node1, node2)
	deleteat!(dd.arc_labels[e], findfirst(dd.arc_labels[e] .== label))
	if isempty(dd.arc_labels[e])
		rem_edge!(dd.graph, e)
		delete!(dd.arc_labels, e)
	end
end

# TODO rem_arcs! is untested
"""Remove all (parallel) arcs between two nodes from a decision diagram"""
function rem_arcs!(dd::DecisionDiagram, node1::Int, node2::Int)
	e = Edge(node1, node2)
	rem_edge!(dd.graph, e)
	delete!(dd.arc_labels, e)
end

"""Return all arc labels between two nodes"""
get_arc_labels(dd::DecisionDiagram, node1::Int, node2::Int) = dd.arc_labels[Edge(node1, node2)]

# if the head node has layer k, the arc layer is k-1 (account for long arcs too)
get_arc_layer(dd::DecisionDiagram, arc::LightGraphs.SimpleGraphs.SimpleEdge{Int64})::Int = dd.node_layers[dst(arc)] - 1

"""Return the in neighbors of a node"""
function inneighbors(dd::DecisionDiagram, node::Int)
	result = []
	for v in inneighbors(dd.graph, node)
		for label in dd.arc_labels[Edge(v, node)]
			push!(result, (v, label))
		end
	end
	return result
end

"""Return the out neighbors of a node"""
function outneighbors(dd::DecisionDiagram, node::Int)
	result = []
	for v in outneighbors(dd.graph, node)
		for label in dd.arc_labels[Edge(node, v)]
			push!(result, (v, label))
		end
	end
	return result
end
