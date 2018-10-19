# Core structure and functions for decision diagrams
# This structure stores both node and arc layers in arrays as fields of DD
# It is fast to compute longest path, but (possibly?!) slow in removing nodes and arcs!
# Commented out sections are not adjusted...

# All higher level functions should use these functions instead of accessing
# the internal structure of a decision diagram directly

# Functions in this file must be careful to keep graph and layers synchronized


using LightGraphs
#using ResumableFunctions
#import LightGraphs.SimpleGraphs.SimpleEdge
import LightGraphs.inneighbors


"""Node structure"""
mutable struct Node{S<:AbstractFloat}
	layer::Int				# node layer
	state::S				# node state
#	out_neighbor::Array{Int,1}		# outgoing arcs of the node, stores arc id in the arc layer
#	in_neighbor::Array{Int,1}		# outgoing arcs of the node, stores arc id in the arc layer
	out_neighbor::Dict{Int, Array{Int, 1}}		# key: the id of node neighbor in the next layer, value: set of arcs (id) connected to that node
	in_neighbor::Dict{Int, Array{Int, 1}}		# key: the id of node neighbor in the previous layer, value: set of arcs (id) connected to that node
	aux_val::Float64		# auxiliary value at the node e.g., the longest path value
end


"""Decision diagram structure"""
struct DecisionDiagram{S<:AbstractFloat}
	graph::SimpleDiGraph{Int}            # represents the directed graph of the DD
    node_layers::Array{Array{Node{S},1},1}        # stores instances of the Node class in each node layer
	arc_layers::Array{Array{Tuple{Int, Int, Float64},1},1}		# stores arc info as tuples (tail, head, label) in arc layers
end


"""
A constrcter for type `DecisionDiagram`

# Input
- `t::Type`: Data type of the state values of nodes
- `nvars::Int`: Number of variables (arc layers)
"""

DecisionDiagram(nvars::Int; t::DataType = Float64)::DecisionDiagram{t} = DecisionDiagram{t}(SimpleDiGraph{Int}(), [Array{Node{t},1}(  0) for i=1:nvars+1], [Array{Tuple{Int, Int, Float64},1}(  0) for i=1:nvars])






"""Return the number of variables (or number of arc layers) in the decision diagram"""
nvars(dd::DecisionDiagram)::Int = length(dd.node_layers) - 1

#nnodes(dd::DecisionDiagram) = nv(dd.graph)
nnodes(dd::DecisionDiagram)::Int = sum(length(dd.node_layers[i]) for i = 1:length(dd.node_layers))

#narcs(dd::DecisionDiagram) = ne(dd.graph)
narcs(dd::DecisionDiagram)::Int = sum(length(dd.arc_layers[i]) for i = 1:length(dd.arc_layers))


"""Set state of an existing node"""
function set_node_state!(dd::DecisionDiagram, layer::Int, id::Int, state::S) where S<:AbstractFloat
	dd.node_layers[layer][id].state = state
end

"""Get state of a node"""
get_node_state(dd::DecisionDiagram, layer::Int, id::Int)::S where S<:AbstractFloat = dd.node_layers[layer][id].state


"""Return arc label"""
get_arc_label(dd::DecisionDiagram, layer::Int, id::Int)::Float64 = dd.arc_layers[layer][id][3]
#get_arc_labels(dd::DecisionDiagram, arc::LightGraphs.SimpleGraphs.SimpleEdge{Int64}) = collect(dd.nodes[src(e)].out_neighbor[dst(e)])


"""Return arc tail and head"""
get_arc_tail(dd::DecisionDiagram, layer::Int, id::Int)::Int = dd.arc_layers[layer][id][1]
get_arc_head(dd::DecisionDiagram, layer::Int, id::Int)::Int = dd.arc_layers[layer][id][2]

get_node_layer_size(dd::DecisionDiagram, layer::Int)::Int = length(dd.node_layers[layer])
get_arc_layer_size(dd::DecisionDiagram, layer::Int)::Int = length(dd.arc_layers[layer])

get_width(dd::DecisionDiagram)::Int = maximum(length(dd.node_layers[i]) for i=1:length(dd.node_layers))

#=
function get_arc_labels(dd::DecisionDiagram, node1::Int, node2::Int)
	i = findfirst(dd.nodes[node1].out_neighbor_nodes, node2)
	#if i == 0
	#	error("No arc exists between node ", node1, " and node ", node2)
	#else
		return dd.nodes[node1].out_neighbor_nodes[i]
	#end
end
=#


"""Return the array of incoming arcs of a node"""
inneighbors(dd::DecisionDiagram, layer::Int, id::Int)::Dict{Int, Array{Int, 1}} = dd.node_layers[layer][id].in_neighbor
#inneighbors(dd::DecisionDiagram, layer::Int, id::Int)::Array{Int, 1} = dd.node_layers[layer][id].in_neighbor

"""Return the array of outgoing nodes of a node"""
outneighbors(dd::DecisionDiagram, layer::Int, id::Int)::Dict{Int, Array{Int, 1}} = dd.node_layers[layer][id].out_neighbor
#outneighbors(dd::DecisionDiagram, layer::Int, id::Int)::Array{Int, 1} = dd.node_layers[layer][id].out_neighbor






"""Add a new node to the decision diagram with a given state and return the node's id at its layer"""
function add_node!(dd::DecisionDiagram, layer::Int, state::S)::Int where S<:AbstractFloat
	#g = dd.graph
	#add_vertex!(g)
	#aux_node::Node{S} = Node(layer, state, Array{Int, 1}(0), Array{Int,1}(0), 0.0)		#This is for the case where neighbor field stores array of arcs only
	aux_node::Node{S} = Node(layer, state, Dict{Int, Array{Int, 1}}(), Dict{Int, Array{Int, 1}}(), 0.0)
	push!(dd.node_layers[layer], aux_node)
	id::Int = length(dd.node_layers[layer])
	return id
end




#=
# TODO: rem_node! is untested
"""Remove a node from the decision diagram"""
function rem_node!(dd::DecisionDiagram, node::Int)

	# ensure the node belongs to the vertex set
	# node in vertices(g) || return false
	node <= length(dd.nodes) || return false


	# Remove incoming and outgoing arcs of the node
	for outnode in keys(dd.nodes[node].out_neighbor)
		delete!(dd.nodes[outnode].in_neighbor, node)
	end
	for innode in keys(dd.nodes[node].in_neighbor)
		delete!(dd.nodes[innode].out_neighbor, node)
	end

	# Move end node in the layer to node's position within the layer
	layer = dd.nodes[node].layer
	layerid = dd.nodes[node].layerid
	dd.layers[layer][layerid] = pop!(dd.layers[layer])
	dd.nodes[dd.layers[layer][layerid]].layerid = layerid

	# Renumber the last node of the graph to the removing node according to LightGraphs implementation
	lastid = length(dd.nodes)
	dd.nodes[node] = pop!(dd.nodes)			# swap the last node with the removed node
	dd.layers[dd.nodes[node].layer][dd.nodes[node].layerid] = node			# update the last node number in dd.layers
	# Update the head/tail number in the incoming/outgoing arcs of the last node
	for outnode in keys(dd.nodes[node].out_neighbor)
		dd.nodes[outnode].in_neighbor[node] = dd.nodes[outnode].in_neighbor[lastid]
		delete!(dd.nodes[outnode].in_neighbor, lastid)
	end
	for innode in keys(dd.nodes[node].in_neighbor)
		dd.nodes[innode].out_neighbor[node] = dd.nodes[innode].out_neighbor[lastid]
		delete!(dd.nodes[innode].out_neighbor, lastid)
	end


	# Remove vertex from graph structure
	#removed = rem_vertex!(dd.graph, node)
	#@assert(removed)
end
=#


#=
"""Remove a node from the decision diagram given its layer and layerid"""
function rem_node!(dd::DecisionDiagram, layer::Int, layerid::Int)
	rem_node!(dd, dd.layers[layer][layerid])
end
=#



"""Add a new arc to the decision diagram between two nodes and retunrs arc id in its layer"""
function add_arc!(dd::DecisionDiagram, layer1::Int, node1_id::Int, layer2::Int, node2_id::Int, label::Float64)::Int
	#@assert((dd.nodes[node1].layer < dd.nodes[node2].layer) && (node1 <= length(dd.nodes)) && (node2 <= length(dd.nodes)))
	#g = dd.graph
	#e = Edge(node1, node2)
	@assert(layer1 == layer2 - 1)
#	if isa(label, Int)					# converts Int label type to floating point
#		lbl = float(label)
#	else
#		lbl = label
#	end
	push!(dd.arc_layers[layer1], (node1_id, node2_id, label))
	id::Int = length(dd.arc_layers[layer1])
#	push!(dd.node_layers[layer1][node1_id].out_neighbor, id)	#for case where neighbors are array of arcs only
#	push!(dd.node_layers[layer2][node2_id].in_neighbor, id)	#for case where neighbors are array of arcs only
	if haskey(dd.node_layers[layer1][node1_id].out_neighbor, node2_id)
		push!(dd.node_layers[layer1][node1_id].out_neighbor[node2_id], id)
		push!(dd.node_layers[layer2][node2_id].in_neighbor[node1_id], id)
	else
		dd.node_layers[layer1][node1_id].out_neighbor[node2_id] = [id]
		dd.node_layers[layer2][node2_id].in_neighbor[node1_id] = [id]
	end
	return id
end



"""
	Add a new arc to the decision diagram between two nodes, where the arc_merging is on. Returns the arc id in the layer.
	If the arc is parallel, it does not add the arc if its label is not lower or upper bound among other arcs.
"""
#NOTE: this code considers the case where inneighbor/outneighbor field is a Dict (NOT array)
function add_reduced_arc!(dd::DecisionDiagram, layer1::Int, node1_id::Int, layer2::Int, node2_id::Int, label::Float64)::Int
	#g = dd.graph
	#e = Edge(node1, node2)
	id::Int = 0			# if the arc is not added, the returning id will be zero
	# check if an arc between two nodes already exists
	if !haskey(dd.node_layers[layer1][node1_id].out_neighbor, node2_id)
		#add_edge!(g, e)
		push!(dd.arc_layers[layer1], (node1_id, node2_id, label))
		id = length(dd.arc_layers[layer1])
		dd.node_layers[layer1][node1_id].out_neighbor[node2_id] = [id]
		dd.node_layers[layer2][node2_id].in_neighbor[node1_id] = [id]
	# if there is only a single arc, the arc is added to in/outneighbor in an increasing order wrt label value
	elseif length(dd.node_layers[layer1][node1_id].out_neighbor[node2_id]) == 1
		push!(dd.arc_layers[layer1], (node1_id, node2_id, label))
		id = length(dd.arc_layers[layer1])
		if dd.arc_layers[layer1][dd.node_layers[layer1][node1_id].out_neighbor[node2_id][1]][3] <= label
			 push!(dd.node_layers[layer1][node1_id].out_neighbor[node2_id], id)	# adds the new label to the parallel arcs
			 push!(dd.node_layers[layer2][node2_id].in_neighbor[node1_id], id)
		 else
			 pushfirst!(dd.node_layers[layer1][node1_id].out_neighbor[node2_id], id)	# adds the new label to the parallel arcs
 			 pushfirst!(dd.node_layers[layer2][node2_id].in_neighbor[node1_id], id)
 		 end
	# if there are two arcs between the nodes, the new arc is added only if it does not belong to the label interval of the current arcs
	# To add the new arc, we only need to change the label of the arc that is being replaced
	else
		@assert(length(dd.node_layers[layer1][node1_id].out_neighbor[node2_id]) == 2)
		id1 = dd.node_layers[layer1][node1_id].out_neighbor[node2_id][1]
		id2 = dd.node_layers[layer1][node1_id].out_neighbor[node2_id][2]
		if label < dd.arc_layers[layer1][id1][3]
			id = id1
			dd.arc_layers[layer1][id] = (node1_id, node2_id, label)
		elseif label > dd.arc_layers[layer1][id2][3]
			id = id2
			dd.arc_layers[layer1][id] = (node1_id, node2_id, label)
		end
	end
	return id
end





#=
# TODO rem_arc! is untested
"""Remove a labeled arc from a decision diagram"""
function rem_arc!(dd::DecisionDiagram, node1::Int, node2::Int, label::Int)
	# checks if the arc exists
	if haskey(dd.nodes[node1].out_neighbor, node2)
		setdiff!(dd.nodes[node1].out_neighbor[node2], label)
		setdiff!(dd.nodes[node2].in_neighbor[node1], label)
		if isempty(dd.nodes[node1].out_neighbor[node2])
			#e = Edge(node1, node2)
			#rem_edge!(dd.graph, e)
			delete!(dd.nodes[node1].out_neighbor, node2)
			delete!(dd.nodes[node2].in_neighbor, node1)
		end
	end
end
=#

#=
# TODO rem_arcs! is untested
"""Remove all (parallel) arcs between two nodes from a decision diagram"""
function rem_arcs!(dd::DecisionDiagram, node1::Int, node2::Int)
	# checks if the arc exists
	if haskey(dd.nodes[node1].out_neighbor, node2)
		#e = Edge(node1, node2)
		#rem_edge!(dd.graph, e)
		delete!(dd.nodes[node1].out_neighbor, node2)
		delete!(dd.nodes[node2].in_neighbor, node1)
	end
end
=#
