include_dependency("core.jl")

#=
# NOTE: This function needs to be adjusted.
"""
Find a longest path on DD.

# Input
- `dd::DecisionDiagram`: The DD.
- `func::Array{Function}`: An array of weight functions (each component represents the separable function for that variable)

# Output
- The longest path (composed of arc labels)
- The length of the longest path
"""
function longest_path(dd::DecisionDiagram, func::Array{Function})

    n = nvars(dd)
    @assert(n == length(func))       #makes sure the size of the fractional point matches the dimension of variables represented by dd

    node_num = nnodes(dd)
    node_obj_vals = Array{Float64}(node_num)   #stores the maximum objective value computed from the source to each node
    node_obj_paths = Array{Int}(node_num)   #at each node, stores the index of the previous node of longest path from source to that node
    node_obj_labels = Array{Int}(node_num)  #stores the label value of the last arc of the longest path from source to each node

    #setting the root objective value equal to zero
    node_obj_vals[root(dd)] = 0

    #computing the longest path from source to a node, by searching through all nodes and their incoming arcs, layer by layer
    for i=2:n+1, node in dd.layers[i]
        temp_v = -Inf
        temp_p = 0
        temp_l = 0
        for parent in inneighbors(dd, node), label in get_arc_labels(dd, parent, node)
            val = node_obj_vals[parent] + func[i-1](label)    #computes the head node value based on the tail and arc label
            if val > temp_v
                temp_v = val
                temp_p = parent
                temp_l = label
            end
        end
        node_obj_vals[node] = temp_v
        node_obj_paths[node] = temp_p
        node_obj_labels[node] = temp_l
    end

    terminal_node = terminal(dd)
    lpv = node_obj_vals[terminal_node]   #the length of the longest path

    lp = Array{Int}(0)
    p = terminal_node
    #computing the arc labels on the longest path through a backtracking
    for i=n:-1:1
        # push!(lp, get_node_obj_label(dd, p))
        push!(lp, node_obj_labels[p])
        # p = get_node_obj_path(dd, p)
        p = node_obj_paths[p]
    end

    return reverse(lp), lpv
end
=#


"""
Second method for function Longest_Path, where weight functions are linear for all variables.

# Input
- `dd::DecisionDiagram`: The DD.
- `cf::Array{T} where T<:Float64`: An array of coefficients of the variables in the linear objective function.

function longest_path(dd::DecisionDiagram, cf::Array{T} where T<:Float64)
    func = Array{Function}(length(cf))
    for (i, c) in enumerate(cf)                 #creates an array of separable functions for each variable
        func[i] = (x::Int64) -> c*x
    end
    return longest_path(dd, func)
end
"""
#This version is for DD structure where we store arc layers, and use them directly to compute the longest path
function longest_path(dd::DecisionDiagram, cf::Array{Float64})::Tuple{Array{Int, 1}, Float64}

    n::Int = nvars(dd)
    @assert(n == length(cf))       #makes sure the size of the fractional point matches the dimension of variables represented by dd


    aux_val::Array{Array{Float64,1}, 1} = [Array{Float64,1}(0) for i=1:n+1]

    #setting the root objective value equal to zero
    push!(aux_val[1], 0)

    #computing the longest path from source to a node, by searching through all nodes and their incoming arcs, layer by layer
    for i::Int=1:n
        aux_val[i+1] = fill(-Inf, get_node_layer_size(dd, i+1))
        for j::Int=1:get_arc_layer_size(dd, i)
            val::Float64 = aux_val[i][get_arc_tail(dd, i, j)] + cf[i]*get_arc_label(dd, i, j)    #computes the head node value based on the tail and arc label
            if val > aux_val[i+1][get_arc_head(dd, i, j)]
                aux_val[i+1][get_arc_head(dd, i, j)] = val
            end
        end
    end

    # Computing the max obj for the terminal node
    lpv::Float64 = aux_val[n+1][1]   #the length of the longest path

    lp::Array{Int, 1} = Array{Int, 1}(0)
    p::Int = 1                   #terminal node id at the last layer
    # Computing the arc labels on the longest path through a backtracking
    for i::Int=n+1:-1:2
        for in_arc::Int in inneighbors(dd, i, p)
            val::Float64 = aux_val[i-1][get_arc_tail(dd, i-1, in_arc)] + cf[i-1]*get_arc_label(dd, i-1, in_arc)    #computes the head node value based on the tail and arc label
            if val == aux_val[i][p]
                push!(lp, get_arc_label(dd, i-1, in_arc))
                p = get_arc_tail(dd, i-1, in_arc)
                break
            end
        end
    end

    return reverse(lp), lpv
end
